# Linear inequality solver

`linineq.py` contains functions for solving linear inequalities in integers. It utilizes methods outlined [here](https://library.wolfram.com/infocenter/Books/8502/AdvancedAlgebra.pdf) on pages 80-81. It combines lattice reduction and linear programming to first reduce the problem to a simpler one and then solve it.

## Functions
All $\le$ denote component-wise comparison.

 - `solve_eq_ineq(M, Mineq, b, bineq)` solves $\mathbf{Mx} = \mathbf{b}$ and $\mathbf{M_{ineq}x} \ge \mathbf{b_{ineq}}$. This is the most general form of the problem.

 - `solve_ineq(Mineq, bineq)` solves $\mathbf{M_{ineq}x} \ge \mathbf{b_{ineq}}$

 - `solve_bounded(M, b, lb, ub)` solves $\mathbf{Mx} = \mathbf{b}$ and $\mathbf{lb} \le \mathbf{x} \le \mathbf{ub}$.

 - `solve_bounded_mod(M, b, lb, ub, N)` solves $\mathbf{Mx} \equiv \mathbf{b} \pmod{N}$ and $\mathbf{lb} \le \mathbf{x} \le \mathbf{ub}$.

 - `solve_bounded_lcg(a, b, m, lb, ub)` solves for $n$ consecutive outputs $\mathbf{s}=(s_0, s_1, ..., s_{n-1})$ of the LCG given by $s_{i+1} \equiv a s_i + b \pmod{m}$ where $\mathbf{lb} \le \mathbf{s} \le \mathbf{ub}$.

 - `solve_truncated_lcg(a, b, m, ys, ntrunc)` solves for $n$ consecutive outputs $\mathbf{s}=(s_0, s_1, ..., s_{n-1})$ of the LCG given by $s_{i+1} \equiv a s_i + b \pmod{m}$ where $(\mathbf{s}>>\mathrm{ntrunc}) = \mathbf{ys}$

Each of these methods has a `_gen` variant (e.g `solve_bounded_gen`) which returns a generator instead of a single solution. This generator will yield all solutions which are findable for a given `lp_bound` parameter.

Keyword arguments:

 - `solver` (default `'ortools'`) is one of `'ortools'` and `'ppl'` and decides which integer programming solver to use. `ortools` is faster and usually preferred, but as it only supports 64-bit numbers it is sometimes insufficient. In those cases `ppl`is used automatically and a warning is emitted. Only `ortools` supports iterating solutions, if you use `ppl` the `_gen` variants will only yield a single solution.

 - `lp_bound` (default `100`) is an internal restriction on the size of the unknowns for the linear programming step, and is currently only used by `ortools`. This does *not* correspond to the size of the solutions you will find, rather the size of the coefficients a lattice-reduced matrix is multiplied by. You can often get by with surprisingly small values for `lp_bound` (e.g. `5`).<br><br>
 If you're unsure what to set `lp_bound` to you can call `set_verbose(1)` before solving an instance. It will then log what internal coefficients where used to find that solution. These will be `<= lp_bound` and `>= -lp_bound`.

 - `reduction` (default `'LLL'`) is one of `'LLL'` and `'BKZ'` and denotes which lattice reduction algorithm to use.

 - `bkz_block_size` (default `20`) is only applicable if `reduction='BKZ'`, and denotes the block size the BKZ algorithm should use.

## Installation
You will need `ortools` (`pip install ortools`) and [Sage](https://doc.sagemath.org/html/en/installation/index.html).

Then either manually download [linineq.py](./linineq.py), or if you're feeling trusting just type `load('https://raw.githubusercontent.com/TheBlupper/linineq/main/linineq.py')` into Sage (no guarantees that functionality won't change in the future though).