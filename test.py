from sage.all import *
from linineq import *
from itertools import product
from traceback import print_exc

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


# test cases from https://github.com/kionactf/coppersmith/blob/1726e06b9bbaac03d8d075e6ba7417d1b360d5ae/lll.py#L452
lattice_tests = [
    ("zerodim", [[0,0,0]]),
    ("onedim", [[1,2,3]]),
    ("twodim_indep", [[1,2,3],[4,5,6]]),
    ("twodim_indep_qq", [[QQ(1)/2,QQ(4)/5,QQ(3)/4],[QQ(1)/4,QQ(5)/3,QQ(3)/5]]),
    ("twodim_dep", [[1,2,3],[2,4,6]]),
    ("threedim_indep", [[1,2,3],[4,5,6],[7,8,8]]),
    ("threedim_one_dep", [[1,2,3],[2,4,6],[8,9,10]]),
    ("threedim_two_dep", [[1,2,3],[2,4,6],[3,6,9]]),
    ("overdim", [[1,2,3],[4,5,6],[7,8,9],[10,11,12]]),
    ("overdim_onedep", [[1,2,3],[4,5,6],[3,6,9],[5,6,7]]),
    ("multiple_2_ker", [[-2,-4,-6],[1,2,3],[3,6,9]]),
]

print('testing reduction transformation and cvp...')
for testname, M in lattice_tests:
    M = matrix(M)
    for red in reds:
        try:
            L, R = red(M, transformation=True)
            assert R*M == L
            assert abs(R.det()) == 1 # unimodular
            assert L.dimensions() == M.dimensions()
        except:
            print_exc()
            print(f'{red.__name__} failed reduction on {testname}')

    for cvp in cvps:
        try:
            c = random_vector(ZZ, M.nrows(), x=-100, y=100)
            u = c*M + random_vector(ZZ, M.ncols(), x=-1, y=1)
            t, cc = cvp(M, u, coords=True)
            assert cc*M == t
        except:
            print_exc()
            print(f'{cvp} failed cvp on {testname}')


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