from linineq import *
from sage.all import *

m = 1 << 52
b = 4164880461924199
a = 2760624790958533
n = 7

s = randrange(m)
ys = []
for _ in range(n):
    s = (a*s+b)%m
    ys.append(s)

ys = vector(ys)
trunc = 45
for sol in solve_truncated_lcg_gen(a, b, m, [y>>trunc for y in ys], trunc):
    print(sol, '<-- correct' if sol == ys else '')