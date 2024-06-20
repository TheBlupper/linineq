from sage.all import *
from linineq import *
from itertools import product

reds = [wBKZ(), wLLL(), wflatter()]
cvps = [wrounding_cvp(), wbabai_cvp(), wfplll_cvp()] + [wkannan_cvp(reduce=red) for red in reds]

mask = 3473400794307473
m = 1 << 52
b = 4164880461924199
a = 2760624790958533
n = 16

lb = [9 * (1 << 52)//10]*n
ub = [m-1]*n

B = random_matrix(ZZ, 50, 100)
t = random_vector(100, 2000)
for cvp in cvps:
    u, v = cvp(B, t, coords=True)
    assert u == v*B


for cvp, red in product(cvps, reds):
    assert {s[0] for s in solve_bounded_lcg_gen(a, b, m, lb, ub, reduce=red, cvp=cvp)} == {4337090406850605, 4359404241011232}
