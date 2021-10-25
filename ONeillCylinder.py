#!/usr/bin/env python3
import math
import numpy as np

a = 9.8  # m/s^2

in2m = 25.4/1000   # 1in = 25.4mm
mi2m = in2m*5280*12

mph2mps = mi2m / 3600 # miles per hour to meters per second

print('mph2mps = ',mph2mps)


Rmi = 56  # miles
Rin = 56*mi2m
R   = Rin*in2m

print('R = ',R,' meters')
print('a = ',a,' acceleration')

v = math.sqrt(R*a)
print('v = ',v,' m/s')
print('v = ',v/mph2mps,' mph')
