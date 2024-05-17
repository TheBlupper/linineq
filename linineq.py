import ppl
import threading
from warnings import warn
from sage.misc.verbose import verbose
from sage.all import vector, ZZ, matrix, identity_matrix, zero_matrix, block_matrix
from ortools.sat.python import cp_model as ort
from queue import Queue


ORTOOLS = 'ortools'
PPL = 'ppl'


def babai(B, t):
    '''
    Returns both the (approximate) closest vector
    to tgt and its coordinates in the already
    reduced lattice M

    (this doesn't really do Babai's algorithm it's just LLL)
    '''

    # fails if B is non-invertible or smth
    # so we skip this and trust the user
    #assert B.is_LLL_reduced()
    
    t = vector(ZZ, t)

    # an LLL reduced basis is ordered 
    # by increasing norm
    S = B[-1].norm().round()+1

    L = block_matrix([
        [B,         0],
        [matrix(t), S]
    ])
    
    L, U = L.LLL(transformation=True)
    for i, u in enumerate(U):
        if abs(u[-1]) == 1:
            # *u[-1] cancels the sign to be positive
            # just in case
            return t - (L[i]*u[-1])[:-1], -u[:-1]
    else:
        raise ValueError('babai failed? plz msg @blupper on discord')

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
    gen = prob.optimizing_point()

    assert gen.is_point() and gen.divisor() == 1
    return tuple(ZZ(c) for c in gen.coefficients())


def gen_solutions(problem, solver=ORTOOLS, lp_bound=100, **_):
    '''
    Return a generator which yields all solutions to a
    problem instance, this will be slower at finding a single
    solution because ortools can't parallelize the search
    '''

    if solver == PPL:
        warn("using ppl which doesn't support enumeration, only one solution "
             "will be returned")
        yield _solve_ppl(problem, lp_bound)
        return
    elif solver != ORTOOLS:
        raise ValueError(f'unknown solver {solver!r}')

    try:
        model, X = _cp_model(problem, lp_bound)
    except OverflowError:
        warn("instance too large for ortools, falling back to ppl which "
             "doesn't support enumeration, only one solution will be returned")
        yield find_solution(problem, solver, lp_bound)
        return

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


def find_solution(problem, solver=ORTOOLS, lp_bound=100, **_):
    '''
    Finds a single solution to a problem instance
    '''

    if solver == PPL:
        return _solve_ppl(problem, lp_bound)
    elif solver != ORTOOLS:
        raise ValueError(f'unknown solver {solver!r}')

    try:
        model, X = _cp_model(problem, lp_bound)
    except OverflowError:
        verbose('instance too large for ortools, falling back to ppl', level=1)
        return _solve_ppl(problem, lp_bound)

    _validate_model(model)

    slvr = ort.CpSolver()
    status = slvr.Solve(model)
    if status == ort.MODEL_INVALID:
        print(model.Validate())
    if status not in [ort.OPTIMAL, ort.FEASIBLE]:
        raise ValueError('No solution found')
    return tuple(slvr.Value(v) for v in X)


# https://library.wolfram.com/infocenter/Books/8502/AdvancedAlgebra.pdf page 80
def _build_system(M, Mineq, b, bineq, reduction='LLL', bkz_block_size=20, **_):
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
    s = V*D.solve_right(U*vector(ZZ, b))
    try:
        s = s.change_ring(ZZ)
    except TypeError:
        raise ValueError('no solution (even without bounds)')

    ker = M.right_kernel_matrix().change_ring(ZZ)

    # switch to left multiplication
    Mker = ker*Mineq.T
    
    if Mker.rank() < Mker.nrows():
        warn('underdetermined inequalities, beware of many solutions')

    # using BKZ instead might help in some cases
    if reduction == 'LLL':
        Mred, R = Mker.LLL(transformation=True)
    elif reduction == 'BKZ':
        Mred = Mker.BKZ(block_size=bkz_block_size)

        # BKZ doesnt provide a transformation matrix
        # Mker is not always invertible hence we use solve_left
        R = Mker.solve_left(Mred).change_ring(ZZ)
    else: raise ValueError(f"reduction must be 'LLL' or 'BKZ', not {reduction!r}")

    bineq = vector(ZZ, bineq) - Mineq*s

    verbose('running babai', level=1)
    bineq_cv, bineq_coord = babai(Mred, bineq)
    bineq -= bineq_cv

    # we then let a solver find an integer solution to
    # x*Mred >= bineq
    
    # precompute the operation (x+v)*R*ker + s
    # as x*T + c
    T = R*ker
    c = bineq_coord*T + s
    
    def f(sol):
        verbose(f'solution paramaters: {sol}', level=1)
        return vector(ZZ, sol)*T + c

    verbose('model processing done', level=1)
    return (Mred, bineq), f


def solve_bounded_gen(M, Mineq, b, bineq, **kwargs):
    '''
    Returns a generetor yielding all* solutions to:
    M*x = b where Mineq*x >= bineq

    *depending on the lp_bound parameter
    '''

    problem, f = _build_system(M, Mineq, b, bineq, **kwargs)
    yield from map(f, gen_solutions(problem, **kwargs))


def solve_bounded(M, Mineq, b, bineq, **kwargs):
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
    yield from solve_bounded_gen(M, Mineq, [], bineq, **kwargs)


def solve_ineq(Mineq, bineq, **kwargs):
    '''
    Finds a solution to:
    Mineq*x >= bineq
    '''

    # 0xn matrix, the kernel will be the nxn identity
    M = matrix(ZZ, 0, Mineq.ncols())
    return solve_bounded(M, Mineq, [], bineq, **kwargs)


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