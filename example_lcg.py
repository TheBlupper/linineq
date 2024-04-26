from sage.all import *
from linineq import solve_bounded_lcg_gen

# finds all possible solutions for any% ssg
# from LA CTF 2024

# the goal was to find a seed for which a given LCG
# outputed 16 consecutive numbers which where all
# above a certain bound

# this can be solved using a normal lattice reduction
# followed by Babai CVP, but finding all possible
# solutions is normally a bit more tricky. Linear
# programming trivializes this

mask = 3473400794307473
m = 1 << 52
b = 4164880461924199
a = 2760624790958533
n = 16

lb = [9 * (1 << 52)//10]*n
ub = [m-1]*n

for sol in solve_bounded_lcg_gen(a, b, m, lb, ub):
    seed = sol[0]
    for _ in range(27):
        seed = pow(a, -1, m)*(seed-b) % m
    print(seed^mask)