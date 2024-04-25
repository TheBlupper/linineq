# Linear inequality solver

`linineq.py` contains functions for solving linear inequalities in integers. It utilizes methods outlined [here](https://library.wolfram.com/infocenter/Books/8502/AdvancedAlgebra.pdf) on pages 80-81. It combines lattice reduction and linear programming to first reduce the problem to a simpler one and then solve it.

## Functions
All $\ge$ denote component-wise comparison.

 - `solve_bounded(M, Mineq, b, bineq)` solves $Mx = b$ and $M_{ineq}x \ge b_{ineq}$. This is the most general form of the problem.
 - `solve_ineq(Mineq, bineq)` solves $M_{ineq}x \ge b_{ineq}$
 - `solve_bounded_mod(M, b, lb, ub, N)` solves $Mx \equiv b \pmod{N}$ and $lb \le x \le ub$.

Each of these methods has a `_gen` variant (e.g `solve_bounded_gen`) which returns a generator instead of a single solution. This generator will yield all solutions which are findable for a given `lp_bound` parameter.

Keyword arguments:

 - `lp_bound` (default `100`) is an internal restriction on the size of the unknowns for the linear programming step. This does *not* correspond to the size of the solutions you will find, rather the size of the coefficients a lattice-reduced matrix is multiplied by. You can often get by with surprisingly small values for `lp_bound` (e.g. `10`).<br><br>
 If you're unsure what to set `lp_bound` to you can call `set_verbose(1)` before solving an instance. It will then log what internal coefficients where used to find that solution. These will be `<= lp_bound` and `>= -lp_bound`.

 - `reduction` (default `'LLL'`) is one of `'LLL'` and `'BKZ'` and denotes which lattice reduction algorithm to use.

 - `bkz_block_size` (default `10`) is only applicable if `reduction='BKZ'`, and denotes the block size the BKZ algorithm should use.

 - `babai_prec` (default `4096`) is the float precision FPLLL should use when running Babai's algorithm, in number of bits. If this is `-1` then a much slower but more precise implementation is used.


## Installation
You will need `ortools` (`pip install ortools`) and [Sage](https://doc.sagemath.org/html/en/installation/index.html).

Then either manually download `linineq.py`, or if you're feeling trusting just type `load('https://raw.githubusercontent.com/TheBlupper/linineq/main/linineq.py')` into Sage.