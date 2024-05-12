import threading
from warnings import warn
from sage.misc.verbose import verbose
from sage.all import vector, ZZ, matrix, identity_matrix, zero_matrix
from ortools.sat.python import cp_model as ort
from queue import Queue
from fpylll import IntegerMatrix, GSO, FPLLL


def babai_fplll(M, tgt, prec=4096):
    '''
    Returns both the (approximate) closest vector
    to tgt and its coordinates in the already
    reduced lattice M
    '''

    prev_prec = FPLLL.get_precision()
    FPLLL.set_precision(prec)

    Mf = IntegerMatrix.from_matrix(M)
    G = GSO.Mat(Mf, float_type='mpfr')
    G.update_gso()
    w = vector(ZZ, G.babai(list(tgt)))

    FPLLL.set_precision(prev_prec)
    return w*M, w


def _validate_model(model):
    err = model.Validate()
    if err == '': return
    if 'Possible integer overflow' in err:
        raise ValueError('Overflow when trying to solve, try decreasing lp_bound')
    raise ValueError(f'Model rejected by ortools: {model.Validate()!r}')


def gen_solutions(model, variables):
    '''
    Return a generator which yields all solutions to an
    ortools model, this will be slower at finding a single
    solution because ortools can't parallelize the search
    '''

    _validate_model(model)

    queue = Queue()
    search_event = threading.Event()
    stop_event = threading.Event()

    slvr = ort.CpSolver()
    slvr.parameters.enumerate_all_solutions = True

    t = threading.Thread(target=_solver_thread,
        args=(model, variables, queue, search_event, stop_event))
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


def _solver_thread(model, variables, queue, search_event, stop_event):
    slvr = ort.CpSolver()
    slvr.parameters.enumerate_all_solutions = True
    solution_collector = _SolutionCollector(variables, queue, search_event, stop_event)
    search_event.wait()
    search_event.clear()
    slvr.Solve(model, solution_collector)
    queue.put(None)


class _SolutionCollector(ort.CpSolverSolutionCallback):
    def __init__(self, vars, queue, search_event, stop_event):
        super().__init__()
        self.vars = vars
        self.queue = queue
        self.search_event = search_event
        self.stop_event = stop_event

    def on_solution_callback(self):
        self.queue.put(tuple(self.Value(v) for v in self.vars))
        self.search_event.wait()
        self.search_event.clear()
        if self.stop_event.is_set():
            self.StopSearch()
            return


def find_solution(model, variables):
    '''
    Finds a single solution to an ortools model
    '''

    _validate_model(model)

    slvr = ort.CpSolver()
    status = slvr.Solve(model)
    if status == ort.MODEL_INVALID:
        print(model.Validate())
    if status not in [ort.OPTIMAL, ort.FEASIBLE]:
        raise ValueError('No solution found')
    return tuple(slvr.Value(v) for v in variables)


# https://library.wolfram.com/infocenter/Books/8502/AdvancedAlgebra.pdf page 80
def _build_system(M, Mineq, b, lb, ub, lp_bound=100, reduction='LLL', bkz_block_size=20, babai_prec=None):
    '''
    Returns a tuple (model, X, f) where model is an ortools model,
    X is a list of variables we want the solution for, and f is a
    function that will transform a solution back to the original space
    '''

    assert Mineq.ncols() == M.ncols()
    assert Mineq.nrows() == len(lb) == len(ub)
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

    # using BKZ instead might help in some cases
    if reduction == 'LLL':
        Mred, R = Mker.LLL(transformation=True)
    elif reduction == 'BKZ':
        Mred = Mker.BKZ(block_size=bkz_block_size)

        # BKZ doesnt provide a transformation matrix
        # Mker is not always invertible hence we use solve_left
        R = Mker.solve_left(Mred).change_ring(ZZ)
    else: raise ValueError(f"reduction must be 'LLL' or 'BKZ', not {reduction!r}")

    Mineq_s = Mineq*s
    lb = vector(ZZ, lb) - Mineq_s
    ub = vector(ZZ, ub) - Mineq_s

    if babai_prec is None:
        # very heuristic, lmk if this causes issues
        babai_prec = max(4096, 2*Mred.norm().round().nbits())

    verbose(f'running babai with {babai_prec} bits of precision', level=1)
    mid = (lb+ub).apply_map(lambda x: x>>1)
    mid_cv, mid_coord = babai_fplll(Mred, mid, prec=int(babai_prec))
    lb_red = lb - mid_cv
    ub_red = ub - mid_cv

    model = ort.CpModel()
    X = [model.NewIntVar(-lp_bound, lp_bound, f'x{i}') for i in range(Mred.nrows())]

    # lb_red <= Mred*X <= ub_red
    Y = [sum([int(c)*x for c, x in zip(col, X)]) for col in Mred.columns()]
    for yi, lbi, ubi in zip(Y, lb_red, ub_red):
        model.Add(int(lbi) <= yi)
        model.Add(yi <= int(ubi))
    
    if Mker.rank() < Mker.nrows():
        warn('underdetermined inequalities, beware of many solutions')
    
    # precompute the operation (x+v)*R*ker + s
    # as x*T + c
    T = R*ker
    c = mid_coord*T + s
    
    def f(sol):
        verbose(f'solution paramaters: {sol}', level=1)
        return vector(ZZ, sol)*T + c

    verbose('model processing done', level=1)
    return model, X, f


def solve_bounded_gen(M, Mineq, b, lb, ub, **kwargs):
    '''
    Returns a generetor yielding all* solutions to:
    M*x = b where lb <= Mineq*x <= ub

    *depending on the lp_bound parameter
    '''

    model, X, f = _build_system(M, Mineq, b, lb, ub, **kwargs)
    yield from map(f, gen_solutions(model, X))


def solve_bounded(M, Mineq, b, lb, ub, **kwargs):
    '''
    Finds a solution to:
    M*x = b where lb <= Mineq*x <= ub
    '''

    model, X, f = _build_system(M, Mineq, b, lb, ub, **kwargs)
    return f(find_solution(model, X))


def solve_ineq_gen(Mineq, lb, ub, **kwargs):
    '''
    Returns a generator yielding all* solutions to:
    lb <= Mineq*x <= ub

    *depending on the lp_bound parameter
    '''

    # 0xn matrix, the kernel will be the nxn identity
    M = matrix(ZZ, 0, Mineq.ncols())
    yield from solve_bounded_gen(M, Mineq, [], lb, ub, **kwargs)


def solve_ineq(Mineq, lb, ub, **kwargs):
    '''
    Finds a solution to:
    lb <= Mineq*x <= ub
    '''

    # 0xn matrix, the kernel will be the nxn identity
    M = matrix(ZZ, 0, Mineq.ncols())
    return solve_bounded(M, Mineq, [], lb, ub, **kwargs)


def _build_mod_system(M, b, lb, ub, N, **kwargs):
    '''
    Returns a tuple (model, X, f) where model is an ortools model,
    X is a list of variables we want the solution for, and f is a
    function that will transform the solution back to the original space
    '''

    neqs = M.nrows()
    nvars = M.ncols()

    M = M.augment(identity_matrix(neqs)*N)
    Mineq = identity_matrix(nvars).augment(zero_matrix(nvars, neqs))

    model, X, f = _build_system(M, Mineq, b, lb, ub, **kwargs)
    return model, X, lambda sol: f(sol)[:nvars]


def solve_bounded_mod_gen(M, b, lb, ub, N, **kwargs):
    '''
    Returns a generator yielding all* solutions to:
    M*x = b (mod N) where lb <= x <= ub

    *depending on the lp_bound parameter
    '''

    model, X, f = _build_mod_system(M, b, lb, ub, N, **kwargs)
    yield from map(f, gen_solutions(model, X))


def solve_bounded_mod(M, b, lb, ub, N, **kwargs):
    '''
    Finds a solution to:
    M*x = b (mod N) where lb <= x <= ub
    '''

    model, X, f = _build_mod_system(M, b, lb, ub, N, **kwargs)
    return f(find_solution(model, X))


def _build_bounded_lcg_system(a, b, m, lb, ub, **kwargs):
    assert len(lb) == len(ub)
    n = len(lb)
    B = vector(ZZ, [(b*(a**i-1)//(a-1))%m for i in range(n)])

    lb = vector(ZZ, lb) - B
    ub = vector(ZZ, ub) - B

    L = identity_matrix(n)*m
    L.set_column(0, [(a**i)%m for i in range(n)])

    # no equalities need to hold
    M = matrix(ZZ, 0, L.ncols())

    model, X, f = _build_system(M, L, [], lb, ub, **kwargs)
    return model, X, lambda sol: L*f(sol)+B


def solve_bounded_lcg(a, b, m, lb, ub, **kwargs):
    '''
    Solves for consecutive outputs of the LCG s[i+1]=(a*s[i]+b)%m
    where lb[i] <= s[i] <= ub[i]
    '''
    model, X, f = _build_bounded_lcg_system(a, b, m, lb, ub, **kwargs)
    return f(find_solution(model, X))


def solve_bounded_lcg_gen(a, b, m, lb, ub, **kwargs):
    '''
    Returns a generator yielding all* possible consecutive
    outputs of the LCG s[i+1]=(a*s[i]+b)%m where
    lb[i] <= s[i] <= ub[i]

    *depending on the lp_bound parameter
    '''
    model, X, f = _build_bounded_lcg_system(a, b, m, lb, ub, **kwargs)
    yield from map(f, gen_solutions(model, X))


def solve_truncated_lcg(a, b, m, ys, ntrunc):
    '''
    Solve for consecutive outputs of the LCG s[i+1]=(a*s[i]+b)%m
    where s[i] >> ntrunc = ys[i]
    '''
    lb = [y << ntrunc for y in ys]
    ub = [((y+1) << ntrunc)-1 for y in ys]
    return solve_bounded_lcg(a, b, m, lb, ub)


def solve_truncated_lcg_gen(a, b, m, ys, ntrunc):
    '''
    Returns a generator yielding all* possible consecutive
    outputs of the LCG s[i+1]=(a*s[i]+b)%m where
    s[i] >> ntrunc = ys[i]

    *depending on the lp_bound parameter
    '''
    lb = [y << ntrunc for y in ys]
    ub = [((y+1) << ntrunc)-1 for y in ys]
    yield from solve_bounded_lcg_gen(a, b, m, lb, ub)