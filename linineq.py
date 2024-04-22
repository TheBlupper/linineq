import threading
from sage.all import vector, ZZ, matrix, identity_matrix, zero_matrix
from ortools.sat.python import cp_model as ort
from queue import Queue


def babai_coords(M, tgt):
    '''
    Returns both the (approximate) closest vector
    to tgt and its coordinates in the already
    reduced lattice M
    '''

    G = M.gram_schmidt()[0]
    diff = tgt

    coord = []
    for i in reversed(range(G.nrows())):
        c = ((diff * G[i]) / (G[i] * G[i])).round()
        coord.append(c)
        diff -= c*M[i]
    return tgt - diff, vector(coord[::-1])


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
def _build_system(M, Mineq, b, bineq, lp_bound=1000):
    '''
    Returns a tuple (model, X, f) where model is an ortools model,
    X is a list of variables we want the solution for, and f is a
    function that will transform a solution back to the original space
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

    ker = M.right_kernel().basis_matrix().change_ring(ZZ)

    Mker = Mineq*ker.T

    # using BKZ instead might help in some cases
    Mred = Mker.T.BKZ().T

    # matrix magic that will transform our
    # solution back to the original space
    R = ((Mker.T*Mker)**-1 * (Mker.T*Mred)).change_ring(ZZ)

    bineq = bineq - Mineq*s
    bineq_cv, v = babai_coords(Mred.T, bineq)
    bineq_red = bineq - bineq_cv

    model = ort.CpModel()
    X = [model.NewIntVar(-lp_bound, lp_bound, f'x{i}') for i in range(Mred.ncols())]

    # Mred*X >= bineq_red
    Y = [sum([c*x for c, x in zip(row, X)]) for row in Mred]
    for i, yi in enumerate(Y):
        model.Add(yi >= bineq_red[i])
    
    # precompute the operation R*(x+v)*ker + s
    # as T*x + c
    T = ker.T*R
    c = T*v + s
    f = lambda sol: T*vector(ZZ, sol) + c

    return model, X, f


def solve_bounded_gen(M, Mineq, b, bineq, **kwargs):
    '''
    Returns a generetor yielding all* solutions to:
    M*x = b where Mineq*x >= bineq

    *depending on the lp_bound parameter
    '''

    model, X, f = _build_system(M, Mineq, b, bineq, **kwargs)
    yield from map(f, gen_solutions(model, X))


def solve_bounded(M, Mineq, b, bineq, **kwargs):
    '''
    Finds a solution to:
    M*x = b where Mineq*x >= bineq
    '''

    model, X, f = _build_system(M, Mineq, b, bineq, **kwargs)
    return f(find_solution(model, X))


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
    '''
    Returns a tuple (model, X, f) where model is an ortools model,
    X is a list of variables we want the solution for, and f is a
    function that will transform the solution back to the original space
    '''

    neqs = M.nrows()
    nvars = M.ncols()

    M = M.augment(identity_matrix(neqs)*N)
    I = identity_matrix(nvars)
    Mineq = I.stack(-I).augment(zero_matrix(nvars*2, neqs))
    bineq = vector([*lb] + [-x for x in ub])

    model, X, f = _build_system(M, Mineq, b, bineq, **kwargs)
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