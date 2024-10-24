import re
from sage.all import *
from linineq import *

# solves "onelinecrypto" from SEETF 2023 by neobeo

n = 23
off = int.from_bytes(b'SEE{' + b'\0'*n + b'}', 'big')

for sol in solve_bounded_mod_gen(
    matrix([256**(n-i) for i in range(n)]),
    vector([0 - off]), [48]*n, [122]*n, 13**37
):
    flag = bytes(sol).decode()
    if re.match(f'\w{{{n}}}', flag):
        break
    print(flag)
else:
    exit('No solution found')

print('found! SEE{' + bytes(sol).decode() + '}')
