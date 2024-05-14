from linineq import *
from sage.all import *

# this problem is way too large for ortools so
# it will automatically fall back to ppl instead

m = 1 << 512
b = randrange(m)
a = randrange(m)
n = 7

s = randrange(m)
ys = []
for _ in range(n):
    s = (a*s+b)%m
    ys.append(s)

ys = vector(ys)
trunc = 100

set_verbose(1)
sol = solve_truncated_lcg(a, b, m, [y>>trunc for y in ys], trunc)
print('recovered:', sol)
print('correct:  ', ys)
print(ys == sol)