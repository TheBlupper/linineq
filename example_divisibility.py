from sage.all import *
from Crypto.Util.number import bytes_to_long
from linineq import solve_ineq

# problem posed by @neobeo on discord:
# what's the smallest integer i such that
# bytes_to_long(str(i).encode()) % 13**37 == 0

# this does not find *the* smallest i, that is unsolved
# afaik. this only finds a single solution which is
# 55 digits long, but you can go further by just throwing
# more computing power at it

# the best i've managed is 42 digits (with lp_bound=2):
# 375100768340560904583476889463350224435184

ndig = 55

L = identity_matrix(ndig)
for i in range(ndig-1):
    L[i, i+1] = -256
L[0,0] = 13**37
L = L[::-1, :]

zero = ord('0')
lb = [zero+1] + [zero]*(ndig-1)
ub = [zero+9]*ndig

Mineq = L.stack(-L)
bineq = lb + [-x for x in ub]
sol = solve_ineq(Mineq, bineq, lp_bound=1, solver='ortools')
n = int(bytes(L*sol))
print(n)
assert bytes_to_long(str(n).encode()) % 13**37 == 0