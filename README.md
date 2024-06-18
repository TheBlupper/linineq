# Linear inequality solver

`linineq.py` contains functions for solving linear inequalities in integers. It utilizes methods outlined [here](https://library.wolfram.com/infocenter/Books/8502/AdvancedAlgebra.pdf) on pages 80-81. It combines lattice reduction and linear programming to first reduce the problem to a simpler one and then solve it.

## Solver functions
All $\le$ denote component-wise comparison.

 - `solve_eq_ineq(M, Mineq, b, bineq)` solves $\mathbf{Mx} = \mathbf{b}$ and $\mathbf{M_{ineq}x} \ge \mathbf{b_{ineq}}$. This is the most general form of the problem.

 - `solve_ineq(Mineq, bineq)` solves $\mathbf{M_{ineq}x} \ge \mathbf{b_{ineq}}$

 - `solve_bounded(M, b, lb, ub)` solves $\mathbf{Mx} = \mathbf{b}$ and $\mathbf{lb} \le \mathbf{x} \le \mathbf{ub}$.

 - `solve_bounded_mod(M, b, lb, ub, N)` solves $\mathbf{Mx} \equiv \mathbf{b}\ (\bmod{\ N})$ and $\mathbf{lb} \le \mathbf{x} \le \mathbf{ub}$.

 - `solve_bounded_lcg(a, b, m, lb, ub)` solves for $n$ consecutive outputs $\mathbf{s}=(s_0, s_1, ..., s_{n-1})$ of the LCG given by $s_{i+1} \equiv a s_i + b \pmod{m}$ where $\mathbf{lb} \le \mathbf{s} \le \mathbf{ub}$.

 - `solve_truncated_lcg(a, b, m, ys, ntrunc)` solves for $n$ consecutive outputs $\mathbf{s}=(s_0, s_1, ..., s_{n-1})$ of the LCG given by $s_{i+1} \equiv a s_i + b \pmod{m}$ where $(\mathbf{s}>>\mathrm{ntrunc}) = \mathbf{ys}$

Each of these methods has a `_gen` variant (e.g `solve_bounded_gen`) which returns a generator instead of a single solution. This generator will yield all solutions which are findable for a given `lp_bound` parameter.

Keyword arguments:

 - `solver` (default `ORTOOLS`) is one of `ORTOOLS` and `PPL` and decides which integer programming solver to use. `ortools` is faster and usually preferred, but as it only supports 64-bit numbers it is sometimes insufficient. In those cases `ppl` is used automatically. Only `ortools` supports iterating solutions.

 - `lp_bound` (default `100`) is an internal restriction on the size of the unknowns for the linear programming step, and is currently only used by `ortools`. This does *not* correspond to the size of the solutions you will find, rather the size of the coefficients a lattice-reduced matrix is multiplied by. You can often get by with surprisingly small values for `lp_bound` (e.g. `5`).<br><br>
 If you're unsure what to set `lp_bound` to you can call `set_verbose(1)` before solving an instance. It will then log what internal coefficients where used to find that solution. These will be `<= lp_bound` and `>= -lp_bound`.

 - `reduce` (default `wLLL()`) is a function which will be used to reduce a lattice basis. The wrapper functions `wLLL`, `wBKZ` and `wflatter` are provided for convenience, use them like `solve_bounded(..., reduce=wBKZ(block_size=20))`. The passed function should behave similarly to the [predefined ones](#lattice-reduction-algorithms).

 - `cvp` (default `wkannan_cvp()`) is function which will be used to solve the (approximate) closest vector problem. The wrapper functions `wkannan_cvp`, `wbabai_cvp` and `wfplll_cvp` are provided for convenience, use them like `solve_bounded(..., cvp=wfplll_cvp(prec=4096))`. The passed function should accept the same parameters as the [these](#cvp-solvers) and behave similarly.

## Lattice reduction algorithms
 - `BKZ(M, transformation=False)` returns the BKZ reduction of $\mathbf{M}$.

 - `flatter(M, path='flatter', transformation=False)` returns a the result of running [flatter](https://github.com/keeganryan/flatter) on $\mathbf{M}$. This requires `flatter` to be installed and in `PATH`, or you can specify the path to the executable with the `path` argument.

 - `LLL(M, transformation=False)` returns the LLL reduction of $\mathbf{M}$. This is just a wrapper around `M.LLL()` and is only here for consistency.

Each of these has a wrapper variant (e.g `wBKZ`) which makes it easy to specify reduction parameters when passing the `reduce` argument to solvers, see the `reduce` parameter above.

The `transformation` parameter indicates if the function should instead return the tuple $(\mathbf{L}, \mathbf{R})$ where $\mathbf{L}$ is the reduced lattice basis and $\mathbf{R M} = \mathbf{L}$.

> [!WARNING]  
> Neither BKZ nor flatter provide a transformation matrix themselves, it is calculated after the fact and this can be slow for large matrices. Use `set_verbose(1)` and look for if it freezes on `computing smith normal form...`

## CVP solvers
 - `kannan_cvp(B, t)` (alias `cvp()`) finds an approximate closest vector to $\mathbf{t}$ in the lattice $\mathbf{B}$ using the Kannan embedding. This uses lattice reduction so the `reduce` argument is relevant even if `is_reduced=True`.

 - `babai_cvp(B, t)` finds an approximate closest vector to $\mathbf{t}$ in the lattice $\mathbf{B}$ using Babai's closest plane algorithm.

 - `rounding_cvp(B, t)` finds an approximate closest vector to $\mathbf{t}$ in the lattice $\mathbf{B}$ using Babai's rounding-off algorithm.

 - `fplll_cvp(B, t, prec=4096)` finds an approximate closest vector to $\mathbf{t}$ in the lattice $\mathbf{B}$ using `fplll`'s `MatGSO.babai()` method.

Each of these accept the following keyword arguments:
 - `reduce` (default `wLLL()`) is a function which will be used to reduce the lattice basis, see `reduce` above.

 - `is_reduced` (default `False`) indicates if the given lattice basis is already reduced, else it will be done internally using `reduce`.

 - `coords` (default `False`) if `True` a tuple $(\mathbf{u}, \mathbf{v})$ will instead be returned where $\mathbf{u} = \mathbf{v B}$ and $\mathbf{u}$ is the closest (approximate) vector to $\mathbf{t}$.

Each of these has a wrapper variant (e.g `wbabai_cvp`) which makes it easy to specify parameters when passing the `cvp` argument to solvers, see the `cvp` parameter above.


## Installation
You will need [Sage](https://doc.sagemath.org/html/en/installation/index.html), and you can optionally `pip install ortools` which allows for enumerating solutions and is often faster than `ppl` (which comes bundled with Sage).

Then either
 - `pip install git+https://github.com/TheBlupper/linineq.git`

 - Type `load('https://raw.githubusercontent.com/TheBlupper/linineq/main/linineq.py')` in Sage

 - Or manually download [linineq.py](./linineq.py)