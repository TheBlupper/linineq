import ppl
import re
import threading
import subprocess
import itertools

from warnings import warn
from queue import Queue
from typing import Callable, Optional
from functools import partial # often comes in handy

from sage.misc.verbose import verbose
from sage.all import ZZ, QQ, vector, matrix, identity_matrix, zero_matrix, block_matrix, xsrange, zero_vector, lcm
from fpylll import IntegerMatrix, GSO, FPLLL

try:
    from ortools.sat.python import cp_model as ort
except ImportError:
    ort = None


ORTOOLS = 'ortools'
PPL = 'ppl'

_DEFAULT_FLATTER_PATH = 'flatter'
_DEFAULT_FPLLL_PATH = 'fplll'

_PROBLEM_LP = 'lp' # base case where we use linear programming
_PROBLEM_UNRESTRICTED = 'unrestricted' # infinite number of solutions
_PROBLEM_ONE_SOLUTION = 'one_solution' # only one solution
_PROBLEM_NO_SOLUTION = 'no_solution'

# LATTICE REDUCTION FUNCTIONS

def _sage_to_fplll(M):
    return '[' + '\n'.join('[' + ' '.join(map(str, row)) + ']' for row in M) + ']'


def _fplll_to_sage(s, nrows, ncols):
    return matrix(ZZ, nrows, ncols, map(ZZ, re.findall(r'-?\d+', s)))


def BKZ(
    M,
    transformation: bool = False,
    no_cli: bool = False,
    block_size: int = 20,
    fplll_path: str = _DEFAULT_FPLLL_PATH,
    auto_abort: bool = True
):
    '''
    Computes the BKZ reduction of the lattice M.
    
    If no_cli is True it will just use M.BKZ(), otherwise it will use the
    fplll CLI which skips having to postcompute the transformation matrix.

    Args:
        M: The input lattice basis.
        transformation (optional): If True, returns the tuple (L, R) where
            L is the reduced basis and R*M = L
        block_size (optional): The block size to use for BKZ.

    Returns:
        The matrix L or the tuple of matrices (L, R)
    '''
    assert block_size >= 1

    M, d = M._clear_denom()

    if M.nrows() > M.ncols() or (not transformation) or no_cli:
        L = M.BKZ()
        if d != 1: L /= d
        if transformation:
            verbose("Retroactively (slowly) computing transformation matrix for BKZ "
                 "because fplll currently doesn't support this for matrices where "
                 "nrows > ncols, see https://github.com/fplll/fplll/issues/525", level=1)
            return L, transformation_matrix(M, L)
        return L
        
    cmd = f'{fplll_path} -a bkz -b {int(block_size)} -of u'.split()
    if auto_abort: cmd.append('-bkzautoabort')

    res = subprocess.check_output(cmd,
        input=_sage_to_fplll(M).encode()).decode()
    R = _fplll_to_sage(res, M.nrows(), M.nrows())
    L = R*M
    return L if d == 1 else L / d, R


def LLL(M, **kwargs):
    '''
    Wrapper around `M.LLL()` for consistency, serves no real purpose.

    Args:
        M: The input lattice basis.
        **kwargs: Passed onto `M.LLL()`

    Returns:
        The result of `M.LLL(**kwargs)`
    '''
    M, d = M._clear_denom()
    M = M.LLL(**kwargs)
    if isinstance(M, tuple): # transformation
        return M[0] if d == 1 else M[0] / d, M[1]
    return M if d == 1 else M / d


_flatter_supports_transformation = None
def flatter(M, transformation: bool=False, path: str=_DEFAULT_FLATTER_PATH):
    '''
    Runs flatter on the lattice basis M using the flatter CLI located
    at `path`.

    Args:
        M: The input lattice basis.
        transformation (optional): If True, returns the tuple (L, R) where
            L is the reduced basis and R*M = L
        path (optional): The path to the flatter CLI.

    Returns:
        The matrix L or the tuple of matrices (L, R)
    '''
    global _flatter_supports_transformation

    M, d = M._clear_denom()

    if M.is_zero():
        if transformation:
            return M, identity_matrix(ZZ, M.nrows())
        return M

    M_rank = M.rank()
    kerdim = M.nrows() - M_rank
    if kerdim > 0:
        # find a basis for the lattice generated by the rows of M
        # this is done by computing the HNF (hermite normal form)
        verbose(f'computing hnf of a {M.nrows()} x {M.ncols()} matrix', level=1)

        M, R = M.echelon_form(transformation=True)
        assert M[M_rank:].is_zero()
        # we put the zero-rows at the top for convenience
        M = M[:M_rank][::-1]
        R = R[::-1]
    else:
        R = identity_matrix(ZZ, M.nrows())
    # R is the transformation into HNF form

    if M.nrows() > 1:
        if transformation and _flatter_supports_transformation is None:
            res = subprocess.check_output([path, '-h'])
            _flatter_supports_transformation = '-of' in res.decode()

            # only warn once
            if not _flatter_supports_transformation:
                warn('the installed version of flatter doesn\'t support providing'
                     ' transformation matrices so it will be calculated after the'
                     ' fact (slowly). consider building the feature branch'
                     ' output_unimodular from the flatter repo')
        verbose(f'running flatter on a {M.nrows()} x {M.ncols()} matrix', level=1)

        if transformation and _flatter_supports_transformation:
            res = subprocess.check_output([path, '-of', 'b', '-of', 'u'],
                input=_sage_to_fplll(M).encode())
            res_lines = res.decode().splitlines()
            L = _fplll_to_sage('\n'.join(res_lines[:M.nrows()+1]), M.nrows(), M.ncols())
            T = _fplll_to_sage('\n'.join(res_lines[M.nrows()+1:]), M.nrows(), M.nrows())
            # T is the transformation into LLL form
        else:
            res = subprocess.check_output([path], input=_sage_to_fplll(M).encode())
            L = _fplll_to_sage(res.decode(), M.nrows(), M.ncols())
    else:
        T = identity_matrix(ZZ, M.nrows())
        L = M

    L = zero_matrix(ZZ, kerdim, L.ncols()).stack(L)
    if d != 1: L /= d

    if transformation:
        if _flatter_supports_transformation:
            return L, R[:kerdim].stack(T*R[kerdim:])
        return L, R[:kerdim].stack(transformation_matrix(M, L[kerdim:])*R[kerdim:])
    return L


def solve_right_int(A, B):
    '''
    Solves the linear integer system of equations
    A*X = B for a matrix A and a vector or matrix B.

    Args:
        A: The matrix A.
        B: The vector or matrix B.

    Returns:
        The vector or matrix X.

    Raises:
        ValueError: If either no rational solution or no integer solution exists.
    '''
    verbose(f'computing smith normal form of a {A.nrows()} x {A.ncols()} matrix', level=1)
    D, U, V = A.smith_form()
    try:
        return (V*D.solve_right(U*B)).change_ring(ZZ)
    except TypeError:
        raise ValueError('no integer solution')
    except ValueError:
        raise ValueError('no solution')


def transformation_matrix(M, L):
    '''
    Finds an integer matrix R s.t R*M = L
    (assuming L is LLL reduced)

    In my experience this produces a unimodular matrix
    but I have no proof or assurance of this.

    Args:
        M: The matrix M.
        L: The matrix L.

    Returns:
        The matrix R.

    Raises:
        ValueError: If no valid transformation matrix can be found, this
                    should not happen if L is an LLL reduction of M.
    '''
    verbose(f'computing transformation matrix of a {M.nrows()} x {M.ncols()} lattice', level=1)
    # this assumes L is LLL reduced...
    kerdim = sum(r.is_zero() for r in L.rows())
    if kerdim != 0:
        R = M.left_kernel_matrix(algorithm='pari')
        assert R.nrows() == kerdim, 'kernel dimension mismatch, is L LLL reduced?'
    else:
        R = matrix(ZZ, 0, M.nrows())
    
    try:
        # ... and so does L[kerdim:]
        return R.stack(solve_right_int(M.T, L[kerdim:].T).T)
    except ValueError:
        raise ValueError('failed to calculate transformation, message @blupper on discord plz')


# CVP FUNCTIONS


def kannan_cvp(B, t, is_reduced: bool=False, reduce: Callable=LLL, coords: bool=False):
    '''
    Computes the (approximate) closest vector to t in the lattice spanned by B
    using the Kannan embedding.

    Args:
        B: The lattice basis.
        t: The target vector.
        is_reduced (optional): If True, assumes B is already reduced.
            Otherwise it will reduce B first.
        reduce (optional): The lattice reduction function to use.
        coords (optional): If True, returns the coordinates of the closest vector.
    
    Returns:
        The closest vector u to t in the lattice spanned by B, or the tuple
        (u, v) where v*B = u if coords is True.
    '''
    if not is_reduced:
        if coords:
            B, R = reduce(B, transformation=True)
        else:
            B = reduce(B)
    elif coords: R = identity_matrix(ZZ, B.nrows())

    t = vector(t)

    # an LLL reduced basis is ordered 
    # by increasing norm
    S = B[-1].norm().round()+1

    L = block_matrix([
        [B,         0],
        [matrix(t), S]
    ])
    
    L, U = reduce(L, transformation=True)
    for u, l in zip(U, L):
        if abs(u[-1]) == 1:
            # *u[-1] cancels the sign to be positive
            # just in case
            res = t - l[:-1]*u[-1]
            if coords:
                return res, -u[:-1]*u[-1]*R
            return res
    raise ValueError("babai failed? plz msg @blupper on discord (unless you didn't reduce?)")


# handy alias
cvp = kannan_cvp


def babai_cvp(B, t, is_reduced: bool=False, reduce: Callable=LLL, coords: bool=False):
    '''
    Computes the (approximate) closest vector to t in the lattice spanned by B
    using Babai's nearest plane algorithm. This can be very slow for large B.

    Args:
        B: The lattice basis.
        t: The target vector.
        is_reduced (optional): If True, assumes B is already reduced.
            Otherwise it will reduce B first.
        reduce (optional): The lattice reduction function to use.
        coords (optional): If True, returns the coordinates of the closest vector.
    
    Returns:
        The closest vector u to t in the lattice spanned by B, or the tuple
        (u, v) where v*B = u if coords is True.
    '''
    if not is_reduced:
        if coords:
            B, R = reduce(B, transformation=True)
        else:
            B = reduce(B)
    elif coords: R = identity_matrix(ZZ, B.nrows())

    if (rank := B.rank()) < B.nrows():
        assert B[:-rank].is_zero(), 'B is not LLL reduced'

    G = B[-rank:].gram_schmidt()[0]
    diff = t

    v = []
    for i in reversed(range(G.nrows())):
        c = ((diff * G[i]) / (G[i] * G[i])).round('even')
        if coords: v.append(c)
        diff -= c*B[-rank+i]
    res = t - diff

    if coords:
        if rank < B.nrows():
            v += [0]*(B.nrows()-rank)
        return res, vector(ZZ, v[::-1])*R
    return res


def fplll_cvp(B, t, prec: int=4096, is_reduced: bool=False, reduce: Callable=LLL, coords: bool=False):
    '''
    Computes the (approximate) closest vector to t in the lattice spanned by B
    using fplll's MatGSO.babai() function. Beware of the precision used.

    Args:
        B: The lattice basis.
        t: The target vector.
        prec (optional): The precision to use for the floating point calculations in fplll.
        is_reduced (optional): If True, assumes B is already reduced.
            Otherwise it will reduce B first.
        reduce (optional): The lattice reduction function to use.
        coords (optional): If True, returns the coordinates of the closest vector.
    
    Returns:
        The closest vector u to t in the lattice spanned by B, or the tuple
        (u, v) where v*B = u if coords is True.
    '''
    if not is_reduced:
        if coords:
            B, R = reduce(B, transformation=True)
        else:
            B = reduce(B)
    elif coords: R = identity_matrix(ZZ, B.nrows())

    prev_prec = FPLLL.get_precision()
    FPLLL.set_precision(prec)

    BZZ, dB = B._clear_denom()
    dT = 1 if t.base_ring() == ZZ else t.denominator()
    d = lcm(dB, dT)
    BZZ *= ZZ(d / dB)
    t = (t*d).change_ring(ZZ)

    Mf = IntegerMatrix.from_matrix(BZZ)
    G = GSO.Mat(Mf, float_type='mpfr')
    G.update_gso()
    v = vector(ZZ, G.babai(list(t)))

    FPLLL.set_precision(prev_prec)
    if coords:
        return v*B, v*R
    return v*B


def rounding_cvp(B, t, is_reduced: bool=False, reduce: Callable=LLL, coords: bool=False):
    '''
    Computes the (approximate) closest vector to t in the lattice spanned by B
    using Babai's rounding-off algorithm. This is the fastest cvp algorithm provided.

    Args:
        B: The lattice basis
        t: The target vector
        is_reduced (optional): If True, assumes B is already reduced.
            Otherwise it will reduce B first.
        reduce (optional): The lattice reduction function to use.
        coords (optional): If True, returns the coordinates of the closest vector.
    
    Returns:
        The closest vector u to t in the lattice spanned by B, or the tuple
        (u, v) where v*B = u if coords is True.
    '''
    if not is_reduced:
        if coords:
            B, R = reduce(B, transformation=True)
        else:
            B = reduce(B)
    elif coords: R = identity_matrix(ZZ, B.nrows())

    if B.is_square() and B.det() != 0:
        exact = B.solve_left(t)
    else:
        # we could also do t*B.pseudoinverse() but it's slower
        exact = (B*B.T).solve_left(t*B.T)

    v = vector(ZZ, [QQ(x).round('even') for x in exact])
    if coords:
        return v*B, v*R
    return v*B


# LINEAR PROGRAMMING SOLVERS


def _cp_model(M, b, lp_bound: int=100):
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


def _solve_ppl(M, b, lp_bound: int=100):
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


def _gen_solutions_ppl(M, b):
    '''
    Uses a branch-and-bound algorithm together with
    PPL as a LP subroutine
    '''
    X = [ppl.Variable(i) for i in range(M.nrows())]

    # pre-constructed system of equations with
    # one more variable eliminated each iteration
    # probably a negligible speedup
    expr_systems = []
    for i in range(len(X)):
        es = []
        for col in M.columns():
            es.append(sum([int(c)*x for c, x in zip(col[i:], X[:len(X)-i])]))
        expr_systems.append(es)

    def recurse(tgt, assignment):
        i = len(assignment)
        if i == len(X):
            yield assignment
            return

        cs = ppl.Constraint_System()
        for expr, lb in zip(expr_systems[i], tgt):
            cs.insert(expr >= int(lb))

        prob = ppl.MIP_Problem(len(X)-i, cs, 0)
        prob.set_objective_function(X[0])

        try:
            prob.set_optimization_mode('minimization')
            lb = QQ(prob.optimal_value()).ceil()
            prob.set_optimization_mode('maximization')
            ub = QQ(prob.optimal_value()).floor()
        except ValueError:
            return # no solution

        for v in range(lb, ub+1):
            yield from recurse(tgt - v*M[i], assignment + (v,))
    yield from recurse(b, ())


def _gen_solutions_ortools(M, b, lp_bound: int=100):
    model, X = _cp_model(M, b, lp_bound) # raises OverflowError

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


def gen_solutions(problem, solver: Optional[str]=None, lp_bound: int=100, **_):
    '''
    Generate solutions to a problem instance. This will be slower at
    finding a single solution because ortools can't parallelize the search.

    Args:
        problem: A problem instance from a _build_xxx function.
        solver (optional): The solver to use, either `ORTOOLS` or `PPL`, or `None` for automatic selection.
        lp_bound (optional): The bounds on the unknown variables in ortools.
    
    Returns:
        A generator yielding solutions to the problem instance.
    '''
    problem_type, params = problem
    if problem_type == _PROBLEM_NO_SOLUTION:
        return
    if problem_type == _PROBLEM_ONE_SOLUTION:
        yield params[0]
        return
    if problem_type == _PROBLEM_UNRESTRICTED:
        ker, s = params
        for v in itertools.product(xsrange(-lp_bound, lp_bound+1), repeat=ker.nrows()):
            yield vector(v)*ker + s
        return
    assert problem_type == _PROBLEM_LP, f'unknown problem type {problem_type!r}'
    M, b, f = params

    if solver == PPL:
        yield from map(f, _gen_solutions_ppl(M, b))
        return
    elif solver != ORTOOLS and solver is not None:
        raise ValueError(f'unknown solver {solver!r}')
    
    if ort is None:
        if solver == ORTOOLS: # explicitly requested ortools
            raise ImportError('ortools is not installed, install with `pip install ortools`')
        verbose('ortools is not installed, falling back to ppl', level=1)
        yield from map(f, _gen_solutions_ppl(M, b))
        return
    
    try:
        yield from map(f, _gen_solutions_ortools(M, b, lp_bound))
    except OverflowError:
        if solver == ORTOOLS: # explicitly requested ortools
            raise ValueError('problem too large for ortools, try using ppl')
        verbose('instance too large for ortools, falling back to ppl', level=1)
        yield from map(f, _gen_solutions_ppl(M, b))


def find_solution(problem, solver: Optional[str]=None, lp_bound: int=100, **_):
    '''
    Find a single solution to a problem instance.

    Args:
        problem: A problem instance from a _build_xxx function.
        solver (optional): The solver to use, either `ORTOOLS`, `PPL` or `None`.
        lp_bound (optional): The bounds on the unknown variables in ortools, not used by ppl.
    
    Returns:
        A tuple of integers representing a solution to the problem instance.
    '''
    problem_type, params = problem

    if problem_type == _PROBLEM_NO_SOLUTION:
        raise ValueError('no solution found')
    if problem_type == _PROBLEM_ONE_SOLUTION:
        return params[0]
    if problem_type == _PROBLEM_UNRESTRICTED:
        ker, s = params
        return s
    assert problem_type == _PROBLEM_LP
    M, b, f = params

    # checks if 0 >= b, in that case 0 is a solution
    # since it's always part of the lattice
    if all(0 >= x for x in b):
        verbose('trivial solution found, no linear programming used', level=1)
        return f([0]*M.nrows())

    if solver == PPL or (solver is None and ort is None):
        return f(_solve_ppl(M, b, lp_bound))
    
    if solver is not None and solver != ORTOOLS:
        raise ValueError(f'unknown solver {solver!r}')
    
    if ort is None:
        raise ImportError('ortools is not installed, install with `pip install ortools`')

    try:
        model, X = _cp_model(M, b, lp_bound)
    except OverflowError:
        # if the user explicitly requested ortools we throw
        if solver == ORTOOLS:
            raise ValueError('problem too large for ortools')
        # otherwise we fall back to ppl
        verbose('instance too large for ortools, falling back to ppl', level=1)
        return f(_solve_ppl(M, b, lp_bound))

    _validate_model(model)

    slvr = ort.CpSolver()
    status = slvr.Solve(model)
    if status == ort.MODEL_INVALID:
        print(model.Validate())
    if status not in [ort.OPTIMAL, ort.FEASIBLE]:
        raise ValueError('no solution found')
    return f([slvr.Value(v) for v in X])


# This is where the magic happens
# based on https://library.wolfram.com/infocenter/Books/8502/AdvancedAlgebra.pdf page 80
def _build_problem(M, Mineq, b, bineq, reduce: Callable = LLL, cvp: Callable = kannan_cvp, kernel_algo: Optional[str] = None, **_):
    '''
    Accepts a system of linear (in)equalities:
    M*x = b where Mineq*x >= bineq

    And reduces the problem as much as possible to make it
    easy for a linear programming solver to solve.

    Args:
        M: The matrix of equalities.
        Mineq: The matrix of inequalities.
        b: The target vector for the equalities.
        bineq: The target vector for the inequalities.
        reduce (optional): The lattice reduction function to use.
        cvp (optional): The CVP function to use.
        kernel_algo (optional): The algorithm to use for finding the kernel of an internal matrix,
            this is passed to the `algorithm` parameter of `Matrix_integer_dense.right_kernel_matrix()`.
            If None it is automatically chosen heuristically.
    
    Returns:
        A tuple containing the problem type and parameters needed to solve it. This should only be passed to
        `find_solution` or `gen_solutions`.
    '''
    assert Mineq.ncols() == M.ncols()
    assert Mineq.nrows() == len(bineq)
    assert M.nrows() == len(b)

    # find unbounded solution
    if M.nrows() == 0:
        s = zero_vector(ZZ, M.ncols())
        ker = identity_matrix(M.ncols())
    else:
        try:
            s = solve_right_int(M, vector(ZZ, b))
        except ValueError:
            return (_PROBLEM_NO_SOLUTION, ())

        if kernel_algo is None:
            # TODO: improve this heuristic
            kernel_algo = 'pari' if M.nrows()/M.ncols() < 0.25 else 'flint'
        ker = M.right_kernel_matrix(algorithm=kernel_algo).change_ring(ZZ)

    # switch to left multiplication
    Mker = ker*Mineq.T

    if Mker.rank() < Mker.nrows():
        warn('underdetermined inequalities, beware of many solutions')

    # degenerate cases
    if Mker.ncols() == 0:
        return (_PROBLEM_UNRESTRICTED, (ker, s))

    if Mker.nrows() == 0:
        verbose('fully determined system', level=1)
        if all(a>=b for a, b in zip(Mineq*s, bineq)):
            return (_PROBLEM_ONE_SOLUTION, (s,))
        return (_PROBLEM_NO_SOLUTION, ())

    Mred, R = reduce(Mker, transformation=True)

    bineq = vector(ZZ, bineq) - Mineq*s

    verbose('running cvp', level=1)
    bineq_cv, bineq_coord = cvp(Mred, bineq, is_reduced=True, coords=True)
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
    return (_PROBLEM_LP, (Mred, bineq, f))


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

    def _to_problem(self, **kwargs):
        dim = len(self.vars)

        M = matrix(ZZ, len(self.eqs_lhs), dim, self.eqs_lhs)
        Mineq = matrix(ZZ, len(self.ineqs_rhs), dim, self.ineqs_lhs)

        problem = _build_problem(M, Mineq, self.eqs_rhs, self.ineqs_rhs, **kwargs)
        return problem, lambda sol: {v: x for v, x in zip(self.vars, sol)}
    
    def solve(self, **kwargs):
        problem, to_dict = self._to_problem(**kwargs)
        return to_dict(find_solution(problem, **kwargs))
    
    def solve_gen(self, **kwargs):
        problem, to_dict = self._to_problem(**kwargs)
        yield from map(to_dict, gen_solutions(problem, **kwargs))


def solve_eq_ineq_gen(M, Mineq, b, bineq, **kwargs):
    '''
    Returns a generetor yielding solutions to:
    M*x = b where Mineq*x >= bineq
    '''
    problem = _build_problem(M, Mineq, b, bineq, **kwargs)
    yield from gen_solutions(problem, **kwargs)


def solve_eq_ineq(M, Mineq, b, bineq, **kwargs):
    '''
    Finds a solution to:
    M*x = b where Mineq*x >= bineq
    '''
    problem = _build_problem(M, Mineq, b, bineq, **kwargs)
    return find_solution(problem, **kwargs)


def solve_ineq_gen(Mineq, bineq, **kwargs):
    '''
    Returns a generator yielding solutions to:
    Mineq*x >= bineq
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


def _build_bounded_problem(M, b, lb, ub, **kwargs):
    assert len(lb) == len(ub) == M.ncols()

    Mineq = identity_matrix(M.ncols())
    Mineq = Mineq.stack(-Mineq)
    bineq = [*lb] + [-x for x in ub]

    return _build_problem(M, Mineq, b, bineq, **kwargs)


def solve_bounded_gen(M, b, lb, ub, **kwargs):
    '''
    Returns a generator yielding solutions to:
    M*x = b where lb <= x <= ub
    '''
    problem = _build_bounded_problem(M, b, lb, ub, **kwargs)
    yield from gen_solutions(problem, **kwargs)


def solve_bounded(M, b, lb, ub, **kwargs):
    '''
    Finds a solution to:
    M*x = b where lb <= x <= ub
    '''
    problem = _build_bounded_problem(M, b, lb, ub, **kwargs)
    return find_solution(problem, **kwargs)


def _build_mod_problem(M, b, lb, ub, N, **kwargs):
    neqs = M.nrows()
    nvars = M.ncols()

    M = M.augment(identity_matrix(neqs)*N)
    Mineq = identity_matrix(nvars).augment(zero_matrix(nvars, neqs))
    Mineq = Mineq.stack(-Mineq)
    bineq = [*lb] + [-x for x in ub]

    problem = _build_problem(M, Mineq, b, bineq, **kwargs)
    return problem, lambda sol: sol[:nvars]


def solve_bounded_mod_gen(M, b, lb, ub, N, **kwargs):
    '''
    Returns a generator yielding solutions to:
    M*x = b (mod N) where lb <= x <= ub
    '''
    problem, f = _build_mod_problem(M, b, lb, ub, N, **kwargs)
    yield from map(f, gen_solutions(problem, **kwargs))


def solve_bounded_mod(M, b, lb, ub, N, **kwargs):
    '''
    Finds a solution to:
    M*x = b (mod N) where lb <= x <= ub
    '''
    problem, f = _build_mod_problem(M, b, lb, ub, N, **kwargs)
    return f(find_solution(problem, **kwargs))


def _build_bounded_lcg_problem(a: int, b: int, m: int, lb, ub, **kwargs):
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

    problem = _build_problem(M, Mineq, [], bineq, **kwargs)

    # also return a function which transforms the solution
    return problem, lambda sol: L*sol+B


def solve_bounded_lcg(a, b, m, lb, ub, **kwargs):
    '''
    Solves for consecutive outputs of the LCG s[i+1]=(a*s[i]+b)%m
    where lb[i] <= s[i] <= ub[i]
    '''
    problem, f = _build_bounded_lcg_problem(a, b, m, lb, ub, **kwargs)
    return f(find_solution(problem, **kwargs))


def solve_bounded_lcg_gen(a, b, m, lb, ub, **kwargs):
    '''
    Returns a generator yielding possible consecutive
    outputs of the LCG s[i+1]=(a*s[i]+b)%m where
    lb[i] <= s[i] <= ub[i]
    '''
    problem, f = _build_bounded_lcg_problem(a, b, m, lb, ub, **kwargs)
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
    Returns a generator yielding possible consecutive
    outputs of the LCG s[i+1]=(a*s[i]+b)%m where
    s[i] >> ntrunc = ys[i]
    '''
    lb = [y << ntrunc for y in ys]
    ub = [((y+1) << ntrunc)-1 for y in ys]
    yield from solve_bounded_lcg_gen(a, b, m, lb, ub, **kwargs)