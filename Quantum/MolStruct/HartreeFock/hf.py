#!/usr/bin/env python3
'''
Implementation of the Hartree-Fock, Self-Consistent Field method for
computing molecular energies.

D. A. Brue
'''
import numpy as np
import math
import scipy as sp
import matplotlib.pyplot as plt


if (__name__ == '__main__'):
    import sys
    rtn = hartreefock()
