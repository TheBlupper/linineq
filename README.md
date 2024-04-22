# Linear inequality solver

`linineq.py` contains functions for solving linear inequalities in integers. It utilizes methods outlined [here](https://library.wolfram.com/infocenter/Books/8502/AdvancedAlgebra.pdf) on pages 80-81. It combines lattice reduction and linear programming to first reduce the problem to a simpler one and then solve it.

All $\ge$ denote component-wise comparison.

 - `solve_bounded(M, Mineq, b, bineq)` solves $Mx = b$ and $M_{ineq}x \ge b_{ineq}$. This is the most general form of the problem.
 - `solve_ineq(Mineq, bineq)` solves $M_{ineq}x \ge b_{ineq}$
 - `solve_bounded_mod(M, b, lb, ub, N)` solves $Mx \equiv b \pmod{N}$ and $lb \le x \le ub$.

Each of these methods has a `_gen` variant (e.g `solve_bounded_gen`) which returns a generator instead of a single solution. This generator will yield all solutions which are findable for the given `lp_bound` parameter.

`lp_bound` is an internal restriction on the size of the unknowns for the linear programming step. This does *not* correspond to the size of the solutions you will find, rather the size of the coefficients a lattice-reduced matrix is multiplied by. You can usually get by with surprisingly small values for `lp_bound` (e.g. 10).

You will need `ortools` (`pip install ortools`) and [Sage](https://doc.sagemath.org/html/en/installation/index.html).