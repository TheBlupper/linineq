from sage.all import *
from Crypto.Cipher import AES
from Crypto.Util.Padding import unpad
from linineq import solve_bounded

# "Lucky Roll" from osu!gaming CTF 2024
# this problem is too large for ortools
# so it will fall back to ppl

p = 4420073644184861649599
a = 1144993629389611207194
b = 3504184699413397958941

# these are outputs of an lcg then taken modulo 100 
trunc = [39, 47, 95, 1, 77, 89, 77, 70, 99, 23, 44, 38, 87, 34, 99, 42, 10, 67, 24, 3]
n = len(trunc)

As = [(a**i) % p for i in range(n)]
Bs = [(b*(a**i-1)//(a-1)) % p for i in range(n)]

M = block_matrix([
    [column_matrix(As),
     identity_matrix(n)*p,
     identity_matrix(n)*100],
])

lb = [0] + [-p]*n + [-p//100]*n
ub = [p] + [p]*n + [p//100]*n

sol = solve_bounded(M, vector([(t-b)%p for t, b in zip(trunc, Bs)]), lb, ub)
s = sol[0]

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