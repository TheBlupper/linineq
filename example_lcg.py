from sage.all import *
from linineq import solve_ineq_gen

# finds all possible solutions for any% ssg
# from LA CTF 2024

# the goal was to find a seed for which a given LCG
# outputed 16 consecutive numbers which where all
# above a certain bound

# this can be solved using a normal lattice reduction
# followed by Babai CVP, but finding all possible
# solutions is normally a bit more tricky

mask = 3473400794307473
m = 1 << 52
b = 4164880461924199
a = 2760624790958533
n = 16

# constant term of each seed state
# (not dependant on s0)
B = vector([(b*(a**i-1)//(a-1)) for i in range(n)])

lb = vector([9 * (1 << 52)//10]*n) - B
ub = vector([m-1]*n)               - B

# modulo m
L = identity_matrix(n)*m

# coefficient in front of s0
L.set_column(0, [(a**i) for i in range(n)])

# one for lower bound one for upper bound
# x <= y  <=>  -x >= -y
Mineq = L.stack(-L)
bineq = vector(ZZ, [*lb] + [-x for x in ub])

for sol in solve_ineq_gen(Mineq, bineq):
    seed = (L*sol+B)[0]
    for _ in range(27):
        seed = pow(a, -1, m)*(seed-b) % m
    print(seed^mask)