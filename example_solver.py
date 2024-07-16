from sage.all import *
from linineq import LinineqSolver

# solve the same problem as example_lcg.py
# but via the LinineqSolver interface

mask = 3473400794307473
m = 1 << 52
b = 4164880461924199
a = 2760624790958533
n = 16

P, (s, *ks) = PolynomialRing(ZZ, 'x', n).objgens()
slvr = LinineqSolver(P)

# we intentionally don't add a +k*m
# term here since it would mean there are
# many identical solutions where s is just
# offset by a multiple of m
out = [s]

# now we can just emulate the LCG symbolically
for i in range(n-1):
    out.append(a*out[-1] + b + ks[i]*m)

for o in out:
    slvr.ge(o, 9*m//10)
    slvr.lt(o, m)

for sol in slvr.solve_gen():
    seed = sol[s]

    for _ in range(27):
        seed = pow(a, -1, m)*(seed-b) % m
    print(seed^mask)