import ppl
import threading
import subprocess

from ast import literal_eval
from warnings import warn
from queue import Queue

from sage.misc.verbose import verbose
from sage.all import vector, ZZ, matrix, identity_matrix, zero_matrix, block_matrix

try:
    from ortools.sat.python import cp_model as ort
except ImportError:
    ort = None


ORTOOLS = 'ortools'
PPL = 'ppl'

FPLLL_PATH = 'fplll'


def _sage_to_fplll(M):
    return '[' + '\n'.join('[' + ' '.join(map(str, row)) + ']' for row in M) + ']'


def _fplll_to_sage(s):
    return matrix(ZZ, literal_eval(s.strip()))


def BKZ(block_size=20):
    assert block_size >= 1
    def wrap(M):
        if M.nrows() > M.ncols():
            raise ValueError('BKZ is broken for matrices where nrows > ncols, '
                             'see https://github.com/fplll/fplll/issues/525')
        res = subprocess.check_output(f'{FPLLL_PATH} -a bkz -b {int(block_size)} -of uk'.split(),
            input=_sage_to_fplll(M).encode()).decode()
        U = _fplll_to_sage(res)
        return U*M, U
    return wrap


def LLL(**kwargs):
    def wrap(M):
        return M.LLL(transformation=True, **kwargs)
    return wrap


def cvp_coords(B, t, is_reduced=False, reduce=LLL()):
    '''
    Returns both the (approximate) closest vector
    to t and its coordinates in the lattice B
    '''

    t = vector(ZZ, t)

    if not is_reduced: B, R = reduce(B)
    else: R = identity_matrix(ZZ, B.nrows())

    # an LLL reduced basis is ordered 
    # by increasing norm
    S = B[-1].norm().round()+1

    L = block_matrix([
        [B,         0],
        [matrix(t), S]
    ])
    
    L, U = reduce(L)
    for u, v in zip(U, L):
        if abs(u[-1]) == 1:
            # *u[-1] cancels the sign to be positive
            # just in case
            return t - v[:-1]*u[-1], -u[:-1]*u[-1]*R
    raise ValueError("babai failed? plz msg @blupper on discord (unless you didn't reduce?)")


def cvp(B, t, is_reduced=False, reduce=LLL()):
    return cvp_coords(B, t, is_reduced, reduce)[0]


def _cp_model(problem, lp_bound=100):
    M, b = problem

    model = ort.CpModel()
    X = [model.NewIntVar(-lp_bound, lp_bound, f'x{i}') for i in range(M.nrows())]

    # X*M >= b
    Y = [sum([int(c)*x for c, x in zip(col, X)]) for col in M.columns()]
    for yi, bi in zip(Y, b):
        model.Add(yi >= int(bi))
    return model, X


def _validate_model(model):
    err = model.Validate()
    if err == '': return
    if 'Possible integer overflow' in err:
        raise ValueError('Possible overflow when trying to solve, try decreasing lp_bound')
    raise ValueError(f'Model rejected by ortools: {model.Validate()!r}')


def _solve_ppl(problem, lp_bound=100):
    M, b = problem

    cs = ppl.Constraint_System()
    X = [ppl.Variable(i) for i in range(M.nrows())]

    # X*M >= b
    Y = [sum([int(c)*x for c, x in zip(col, X)]) for col in M.columns()]
    for yi, bi in zip(Y, b):
        cs.insert(yi >= int(bi))

    # TODO: check if lp_bound would improve/decrease performance of ppl
    # in my experience it performs worse hence it's not used atm
    #for i in range(len(X)):
    #    cs.insert(X[i] <= lp_bound)
    #    cs.insert(-lp_bound<=X[i])

    prob = ppl.MIP_Problem(len(X), cs, 0)
    prob.add_to_integer_space_dimensions(ppl.Variables_Set(X[0], X[-1]))
    try:
        gen = prob.optimizing_point()
    except ValueError:
        raise ValueError('no solution found')

    assert gen.is_point() and gen.divisor() == 1
    return tuple(ZZ(c) for c in gen.coefficients())


def gen_solutions(problem, solver=None, lp_bound=100, **_):
    '''
    Return a generator which yields all solutions to a
    problem instance, this will be slower at finding a single
    solution because ortools can't parallelize the search
    '''

    # solver parameter is a little redundant here since
    # only ortools is supported, kept here for future possibilities

    if solver == PPL:
        raise ValueError('ppl does not support enumeration')
    elif solver is not None and solver != ORTOOLS:
        raise ValueError(f'unknown solver {solver!r}')
    
    if ort is None:
        raise ImportError('ortools is needed for enumeration but is not installed,'
                         ' install with `pip install ortools`')

    try:
        model, X = _cp_model(problem, lp_bound)
    except OverflowError:
        raise ValueError('problem too large for ortools, enumeration not possible')

    _validate_model(model)

    queue = Queue()
    search_event = threading.Event()
    stop_event = threading.Event()

    slvr = ort.CpSolver()
    slvr.parameters.enumerate_all_solutions = True

    t = threading.Thread(target=_solver_thread,
        args=(model, X, queue, search_event, stop_event))
    t.start()

    try:
        while True:
            search_event.set()
            sol = queue.get()
            if sol is None:
                break
            yield sol
    finally:
        stop_event.set()
        search_event.set()
        t.join()


def _solver_thread(model, X, queue, search_event, stop_event):
    slvr = ort.CpSolver()
    slvr.parameters.enumerate_all_solutions = True
    solution_collector = _SolutionCollector(X, queue, search_event, stop_event)
    search_event.wait()
    search_event.clear()
    slvr.Solve(model, solution_collector)
    queue.put(None)


if ort is not None:
    class _SolutionCollector(ort.CpSolverSolutionCallback):
        def __init__(self, X, queue, search_event, stop_event):
            super().__init__()
            self.X = X
            self.queue = queue
            self.search_event = search_event
            self.stop_event = stop_event

        def on_solution_callback(self):
            self.queue.put(tuple(self.Value(x) for x in self.X))
            self.search_event.wait()
            self.search_event.clear()
            if self.stop_event.is_set():
                self.StopSearch()
                return


def find_solution(problem, solver=None, lp_bound=100, **_):
    '''
    Finds a single solution to a problem instance
    '''

    # checks if 0 >= b, in that case 0 is a solution
    # since it's always part of the lattice
    if all(x <= 0 for x in problem[1]):
        verbose('trivial solution found, no linear programming used', level=1)
        return tuple([0]*problem[0].nrows())

    if solver == PPL or (solver is None and ort is None):
        return _solve_ppl(problem, lp_bound)
    
    if solver is not None and solver != ORTOOLS:
        raise ValueError(f'unknown solver {solver!r}')
    
    if ort is None:
        raise ImportError('ortools is not installed, install with `pip install ortools`')

    try:
        model, X = _cp_model(problem, lp_bound)
    except OverflowError:
        # if the user explicitly requested ortools we throw
        if solver == ORTOOLS:
            raise ValueError('problem too large for ortools')
        # otherwise we fall back to ppl
        verbose('instance too large for ortools, falling back to ppl', level=1)
        return _solve_ppl(problem, lp_bound)

    _validate_model(model)

    slvr = ort.CpSolver()
    status = slvr.Solve(model)
    if status == ort.MODEL_INVALID:
        print(model.Validate())
    if status not in [ort.OPTIMAL, ort.FEASIBLE]:
        raise ValueError('no solution found')
    return tuple(slvr.Value(v) for v in X)


# https://library.wolfram.com/infocenter/Books/8502/AdvancedAlgebra.pdf page 80
def _build_system(M, Mineq, b, bineq, reduce=LLL(), **_):
    '''
    Returns a tuple (problem, f) where problem is a tuple of the form (M, b)
    where a solution to x*M >= b is sought, and f is a function that will
    transform a solution back to the original space of the query
    '''

    assert Mineq.ncols() == M.ncols()
    assert Mineq.nrows() == len(bineq)
    assert M.nrows() == len(b)

    # find unbounded solution
    D, U, V = M.smith_form()
    try:
        s = V*D.solve_right(U*vector(ZZ, b))
        s = s.change_ring(ZZ)
    except (TypeError, ValueError):
        raise ValueError('no solution (even without bounds)')

    # TODO faster kernel finding
    ker = M.right_kernel_matrix().change_ring(ZZ)

    # switch to left multiplication
    Mker = ker*Mineq.T

    if Mker.rank() < Mker.nrows():
        warn('underdetermined inequalities, beware of many solutions')

    # degenerate case
    if Mker.ncols() == 0: 
        Mred = matrix([[0]*ker.nrows()]).T
        bred = vector([0])

        def f(sol):
            verbose(f'solution paramaters: {sol}', level=1)
            return vector(ZZ, sol)*ker + s
    
        return (Mred, bred), f

    if Mker.nrows() == 0:
        warn('fully determined equalities, no enumeration will occur')
        Mred = matrix([[1, -1]])

        if all(a>=b for a, b in zip(Mineq*s, bineq)):
            bred = vector([0, 0]) # one solution
        else:
            bred = vector([0, 1]) # no solutions

        return (Mred, bred), lambda _: s

    Mred, R = reduce(Mker)

    bineq = vector(ZZ, bineq) - Mineq*s

    verbose('running cvp', level=1)
    bineq_cv, bineq_coord = cvp_coords(Mred, bineq, is_reduced=True, reduce=reduce)
    bineq -= bineq_cv

    # we then let a solver find an integer solution to
    # x*Mred >= bineq
    
    # precompute the operation (x+bineq_coord)*R*ker + s
    # as x*T + c
    T = R*ker
    c = bineq_coord*T + s
    
    def f(sol):
        verbose(f'solution paramaters: {sol}', level=1)
        return vector(ZZ, sol)*T + c

    verbose('model processing done', level=1)
    return (Mred, bineq), f


class LinineqSolver:
    '''
    Helper class for setting up systems of
    linear equalities/inequalities and solving them

    See example usage in example_solver.py
    '''
    def __init__(self, polyring):
        assert polyring.base_ring() is ZZ
        self.polyring = polyring
        self.vars = polyring.gens()

        self.eqs_lhs = [] # matrix of integers
        self.eqs_rhs = [] # vector of integers

        self.ineqs_lhs = [] # matrix of integers
        self.ineqs_rhs = [] # vector of integers
    
    def eq(self, lhs, rhs):
        poly = self.polyring(lhs - rhs)
        assert poly.degree() <= 1
        self.eqs_lhs.append([poly.coefficient(v) for v in self.vars])
        self.eqs_rhs.append(-poly.constant_coefficient())

    def ge(self, lhs, rhs):
        poly = self.polyring(lhs - rhs)
        assert poly.degree() <= 1
        self.ineqs_lhs.append([poly.coefficient(v) for v in self.vars])
        self.ineqs_rhs.append(-poly.constant_coefficient())

    def gt(self, lhs, rhs):
        self.ge(lhs, rhs+1)
    
    def le(self, lhs, rhs):
        self.ge(rhs, lhs)

    def lt(self, lhs, rhs):
        self.ge(rhs, lhs+1)

    def _to_system(self, **kwargs):
        dim = len(self.vars)

        M = matrix(ZZ, len(self.eqs_lhs), dim, self.eqs_lhs)
        Mineq = matrix(ZZ, len(self.ineqs_rhs), dim, self.ineqs_lhs)

        problem, f = _build_system(M, Mineq, self.eqs_rhs, self.ineqs_rhs, **kwargs)
        return problem, lambda sol: {v: x for v, x in zip(self.vars, f(sol))}
    
    def solve(self, **kwargs):
        problem, f = self._to_system(**kwargs)
        return f(find_solution(problem, **kwargs))
    
    def solve_gen(self, **kwargs):
        problem, f = self._to_system(**kwargs)
        yield from map(f, gen_solutions(problem, **kwargs))


def solve_eq_ineq_gen(M, Mineq, b, bineq, **kwargs):
    '''
    Returns a generetor yielding all* solutions to:
    M*x = b where Mineq*x >= bineq

    *depending on the lp_bound parameter
    '''

    problem, f = _build_system(M, Mineq, b, bineq, **kwargs)
    yield from map(f, gen_solutions(problem, **kwargs))


def solve_eq_ineq(M, Mineq, b, bineq, **kwargs):
    '''
    Finds a solution to:
    M*x = b where Mineq*x >= bineq
    '''

    problem, f = _build_system(M, Mineq, b, bineq, **kwargs)
    return f(find_solution(problem, **kwargs))


def solve_ineq_gen(Mineq, bineq, **kwargs):
    '''
    Returns a generator yielding all* solutions to:
    Mineq*x >= bineq

    *depending on the lp_bound parameter
    '''

    # 0xn matrix, the kernel will be the nxn identity
    M = matrix(ZZ, 0, Mineq.ncols())
    yield from solve_eq_ineq_gen(M, Mineq, [], bineq, **kwargs)


def solve_ineq(Mineq, bineq, **kwargs):
    '''
    Finds a solution to:
    Mineq*x >= bineq
    '''

    # 0xn matrix, the kernel will be the nxn identity
    M = matrix(ZZ, 0, Mineq.ncols())
    return solve_eq_ineq(M, Mineq, [], bineq, **kwargs)


def _build_bounded_system(M, b, lb, ub, **kwargs):
    assert len(lb) == len(ub) == M.ncols()

    Mineq = identity_matrix(M.ncols())
    Mineq = Mineq.stack(-Mineq)
    bineq = [*lb] + [-x for x in ub]

    return _build_system(M, Mineq, b, bineq, **kwargs)


def solve_bounded_gen(M, b, lb, ub, **kwargs):
    '''
    Returns a generator yielding all* solutions to:
    M*x = b where lb <= x <= ub

    *depending on the lp_bound parameter
    '''

    problem, f = _build_bounded_system(M, b, lb, ub, **kwargs)
    yield from map(f, gen_solutions(problem, **kwargs))


def solve_bounded(M, b, lb, ub, **kwargs):
    '''
    Finds a solution to:
    M*x = b where lb <= x <= ub
    '''

    problem, f = _build_bounded_system(M, b, lb, ub, **kwargs)
    return f(find_solution(problem, **kwargs))


def _build_mod_system(M, b, lb, ub, N, **kwargs):
    neqs = M.nrows()
    nvars = M.ncols()

    M = M.augment(identity_matrix(neqs)*N)
    Mineq = identity_matrix(nvars).augment(zero_matrix(nvars, neqs))
    Mineq = Mineq.stack(-Mineq)
    bineq = [*lb] + [-x for x in ub]

    problem, f = _build_system(M, Mineq, b, bineq, **kwargs)
    return problem, lambda sol: f(sol)[:nvars]


def solve_bounded_mod_gen(M, b, lb, ub, N, **kwargs):
    '''
    Returns a generator yielding all* solutions to:
    M*x = b (mod N) where lb <= x <= ub

    *depending on the lp_bound parameter
    '''

    problem, f = _build_mod_system(M, b, lb, ub, N, **kwargs)
    yield from map(f, gen_solutions(problem, **kwargs))


def solve_bounded_mod(M, b, lb, ub, N, **kwargs):
    '''
    Finds a solution to:
    M*x = b (mod N) where lb <= x <= ub
    '''

    problem, f = _build_mod_system(M, b, lb, ub, N, **kwargs)
    return f(find_solution(problem, **kwargs))


def _build_bounded_lcg_system(a, b, m, lb, ub, **kwargs):
    assert len(lb) == len(ub)
    n = len(lb)
    B = vector(ZZ, [(b*(a**i-1)//(a-1))%m for i in range(n)])

    lb = vector(ZZ, lb) - B
    ub = vector(ZZ, ub) - B
    bineq = [*lb] + [-x for x in ub]

    L = identity_matrix(n)*m
    L.set_column(0, [(a**i)%m for i in range(n)])
    Mineq = L.stack(-L)

    # no equalities need to hold
    M = matrix(ZZ, 0, L.ncols())

    problem, f = _build_system(M, Mineq, [], bineq, **kwargs)
    return problem, lambda sol: L*f(sol)+B


def solve_bounded_lcg(a, b, m, lb, ub, **kwargs):
    '''
    Solves for consecutive outputs of the LCG s[i+1]=(a*s[i]+b)%m
    where lb[i] <= s[i] <= ub[i]
    '''
    problem, f = _build_bounded_lcg_system(a, b, m, lb, ub, **kwargs)
    return f(find_solution(problem, **kwargs))


def solve_bounded_lcg_gen(a, b, m, lb, ub, **kwargs):
    '''
    Returns a generator yielding all* possible consecutive
    outputs of the LCG s[i+1]=(a*s[i]+b)%m where
    lb[i] <= s[i] <= ub[i]

    *depending on the lp_bound parameter
    '''
    problem, f = _build_bounded_lcg_system(a, b, m, lb, ub, **kwargs)
    yield from map(f, gen_solutions(problem, **kwargs))


def solve_truncated_lcg(a, b, m, ys, ntrunc, **kwargs):
    '''
    Solve for consecutive outputs of the LCG s[i+1]=(a*s[i]+b)%m
    where s[i] >> ntrunc = ys[i]
    '''
    lb = [y << ntrunc for y in ys]
    ub = [((y+1) << ntrunc)-1 for y in ys]
    return solve_bounded_lcg(a, b, m, lb, ub, **kwargs)


def solve_truncated_lcg_gen(a, b, m, ys, ntrunc, **kwargs):
    '''
    Returns a generator yielding all* possible consecutive
    outputs of the LCG s[i+1]=(a*s[i]+b)%m where
    s[i] >> ntrunc = ys[i]

    *depending on the lp_bound parameter
    '''
    lb = [y << ntrunc for y in ys]
    ub = [((y+1) << ntrunc)-1 for y in ys]
    yield from solve_bounded_lcg_gen(a, b, m, lb, ub, **kwargs)