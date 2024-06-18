from sage.all import matrix
from linineq import solve_bounded_mod_gen

# adapted from https://connor-mccartney.github.io/cryptography/other/Trying-to-crack-COD-FNV-hashes
FNV_INIT = 0xCBF29CE484222325
p = 0x100000001B3

def fnv64(s):
    hsh = FNV_INIT
    for c in s.lower().replace(b"\\", b"/"):
        hsh = ((hsh^c)*p)%2**64
    return hsh

def rev(sol):
    ret = []
    h = FNV_INIT
    for s in sol:
        ch = (h+s)^h
        if ch not in range(32, 128):
            return None
        ret.append(ch)
        h += s
        h *= p
    return bytes(ret)

def solve(target, n):
    hsh = FNV_INIT
    rets = []
    M = matrix([[p**(n - i) for i in range(n)]])

    # lp_bound could maybe be increased slightly if the
    # solution isn't found
    for sol in solve_bounded_mod_gen(
        M, [target - hsh*p**n], [-128]*n, [128]*n, 2**64, lp_bound=10
    ):
        ret = rev(sol)
        if ret is None: continue
        print(ret)
        rets.append(ret)
    return rets


print(solve(fnv64(b'abcdefghi'), 9))

# this takes a while
print(solve(fnv64(b'abcdefghij'), 10))