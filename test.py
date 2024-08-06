from sage.all import *
from linineq import *
from itertools import product

reds = [BKZ, LLL, flatter]
cvps = [rounding_cvp, babai_cvp, fplll_cvp] + [partial(kannan_cvp, reduce=red) for red in reds]

print('testing edge cases...')

# fully determined equations
Meq = random_matrix(ZZ, 10, 10, algorithm='unimodular')
Mineq = identity_matrix(ZZ, Meq.ncols())
x = vector(ZZ, [1337]*Meq.ncols())
a = Meq*x

# should always succeed
assert solve_eq_ineq(Meq, Mineq, a, [0]*Meq.nrows()) == x

# should never succeed
try:
    solve_eq_ineq(Meq, -Mineq, a, [0]*Meq.nrows())
except ValueError: pass
else: assert False, 'wait how'

# unbounded
Meq = random_matrix(ZZ, 1, 10)
Mineq = matrix(ZZ, 0, Meq.ncols())
x = vector(ZZ, [1337]*Meq.ncols())
a = Meq*x
assert len(list(solve_eq_ineq_gen(Meq, Mineq, a, [], lp_bound=1))) == 3**9

print('testing reduction transformation...')
for M in [random_matrix(ZZ, r, c) for r, c in [(50, 50), (50, 60), (60, 50)]]:
    for red in reds:
        L, R = red(M, transformation=True)
        assert R*M == L
        assert abs(R.det()) == 1 # unimodular


print('testing cvp...')
B = random_matrix(ZZ, 50, 100)
t = random_vector(100, 2000)
for cvp in cvps:
    u, v = cvp(B, t, coords=True)
    assert u == v*B

print('testing bounded lcg...')
mask = 3473400794307473
m = 1 << 52
b = 4164880461924199
a = 2760624790958533
n = 16

lb = [9 * (1 << 52)//10]*n
ub = [m-1]*n

for cvp, red in product(cvps, reds):
    assert {s[0] for s in solve_bounded_lcg_gen(a, b, m, lb, ub, reduce=red, cvp=cvp)} == {4337090406850605, 4359404241011232}

print('all good :)')