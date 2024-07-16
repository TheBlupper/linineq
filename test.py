from sage.all import *
from linineq import *
from itertools import product

reds = [BKZ, LLL, flatter]
cvps = [rounding_cvp, babai_cvp, fplll_cvp] + [partial(kannan_cvp, reduce=red) for red in reds]

B = random_matrix(ZZ, 50, 100)
t = random_vector(100, 2000)
for cvp in cvps:
    u, v = cvp(B, t, coords=True)
    assert u == v*B


mask = 3473400794307473
m = 1 << 52
b = 4164880461924199
a = 2760624790958533
n = 16

lb = [9 * (1 << 52)//10]*n
ub = [m-1]*n

for cvp, red in product(cvps, reds):
    assert {s[0] for s in solve_bounded_lcg_gen(a, b, m, lb, ub, reduce=red, cvp=cvp)} == {4337090406850605, 4359404241011232}
