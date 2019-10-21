"""Test matlab command form python."""
import matlab.engine

eng = matlab.engine.start_matlab()
# tf = eng.isprime(37)
# print(tf)

net = eng.network_swing_simple()

print(net.N)
