from sage.all import *
from Crypto.Cipher import AES
from Crypto.Util.Padding import unpad
from linineq import *

# "Lucky Roll" from osu!gaming CTF 2024
# this problem is too large for ortools
# so it will fall back to ppl

p = 4420073644184861649599
a = 1144993629389611207194
b = 3504184699413397958941

# these are outputs of an lcg then taken modulo 100 
trunc = [39, 47, 95, 1, 77, 89, 77, 70, 99, 23, 44, 38, 87, 34, 99, 42, 10, 67, 24, 3, 2, 80, 26, 87, 91, 86, 1, 71, 59, 97, 69, 31, 17, 91, 73, 78, 43, 18, 15, 46, 22, 68, 98, 60, 98, 17, 53, 13, 6, 13, 19, 50, 73, 44, 7, 44, 3, 5, 80, 26, 10, 55, 27, 47, 72, 80, 53, 2, 40, 64, 55, 6]
n = len(trunc)

P, (sym_s, *ks) = PolynomialRing(ZZ, 'x', 2*n).objgens()
k100 = ks[:n]
kp = ks[n:]

slvr = LinineqSolver(P)
out = [sym_s]

# now we can just emulate the LCG symbolically
for i in range(1, n):
    out.append(a*out[-1] + b)

for i in range(n):
    # mod 100
    out[i] -= k100[i]*100

    # mod p, but not on the first one to avoid
    # multiple solutions
    if i != 0: out[i] -= kp[i-1]*p

for i, o in enumerate(out):
    slvr.gt(k100[i], -p//100)
    slvr.lt(k100[i], p//100)
    slvr.eq(o, trunc[i])

s = slvr.solve()[sym_s]

def lcg():
    global s
    s = (a*s + b) % p
    return s % 100

for _ in range(72-1): lcg()

ct = bytes.fromhex('34daaa9f7773d7ea4d5f96ef3dab1bbf5584ecec9f0542bbee0c92130721d925f40b175e50587196874e14332460257b')
key = bytes([lcg() for _ in range(16)])
iv = bytes([lcg() for _ in range(16)])
cipher = AES.new(key, AES.MODE_CBC, iv)
print(unpad(cipher.decrypt(ct), 16).decode())