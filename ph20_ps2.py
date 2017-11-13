#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Ph20 Problem Set #2
Morgaine Mandigo-Stoba

This is the edit I am making.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.integrate import romberg

# Problem #2
def ext_trap(func, a, b, N):
    '''
    Inputs:
        func - a function to be integrated
        a - starting value of the integral
        b - ending value of the integral
        N - number of the steps to be used in the trapezoid rule approximation
    
    Output:
        total - the approximate value of the integral according to the 
        extended trapezoid rule     
    '''
    hn = (float((b - a)) / N)
    steps = np.linspace(a, b, N)
    new_func = np.vectorize(func)
    vals = new_func(steps)
    vals[0] = (vals[0] / 2.0)
    vals[-1] = (vals[-1] / 2.0)
    total = hn * (np.sum(vals))
    return total

# Problem #3
def ext_simp(func, a, b, N):
    '''
    Inputs:
        func - a functio to be integrated
        a - starting value of the integral
        b - ending value of the integral
        N - number of the steps to be used in the simpson's rule approximation
        
    Output: 
        total - the approximate value of the integral according to the 
        extended simpson's rule
    '''
    hn = (float((b - a)) / N)
    steps = np.linspace(a, b, ((2 * N) +1))
    new_func = np.vectorize(func)
    vals = new_func(steps)
    vals = (vals / 3.0)
    for i in range(1, len(vals) - 1):
        if (i % 2) == 1:
            vals[i] = (vals[i] * 2.0)
    vals[0] = (vals[0] / 2.0)
    vals[-1] = (vals[-1] / 2.0)
    total = hn * (np.sum(vals))
    return total

# Problem #4
def ex(x):
    '''
    Inputs:
        x - an x value
    
    Output:
        y - the corresponding y value of the equation e^x
    '''
    y = np.exp(x)
    return y

real_val = np.e - 1
N_vals = np.arange(1, 5000, 100)
trap_vals = np.array([])
simp_vals = np.array([])
for val in N_vals:
    trap_vals = np.append(trap_vals, (ext_trap(ex, 0, 1, val)))
    simp_vals = np.append(simp_vals, (ext_simp(ex, 0, 1, val)))
trap_err = np.absolute(trap_vals - real_val)
simp_err = np.absolute(simp_vals - real_val)

plt.plot(N_vals, trap_err)
plt.loglog()
plt.xlabel('N')
plt.ylabel('Error')
plt.title('Trapezoid Rule')
plt.show()

plt.plot(N_vals, simp_err)
plt.loglog()
plt.xlabel('N')
plt.ylabel('Error')
plt.title("Simpson's Rule")
plt.show()

# Problem #6
def comp_rel_error(val1, val2):
    '''
    Inputs:
        val1 - the first value to compare
        val2 - the second value to compare
    Output:
        error - the difference in the two values in reference to val1
    '''
    error = (np.absolute(val1 - val2)) / val1
    return error
    
def rel_error(func, a, b, des_err):
    '''
    Inputs:
        func - the function to be integrated
        des_err - the desired relative error, at which the process will stop
    
    Output:
        total - the approximate value of the integral once the specified
        error tolerance has been met
    '''
    N = 10
    val1 = ext_simp(func, a, b, N)
    val2 = ext_simp(func, a, b, (2 * N))
    while comp_rel_error(val1, val2) > des_err:
        N *= 2
        val1 = ext_simp(func, a, b, N)
        val2 = ext_simp(func, a, b, (2 * N))
    total = ext_simp(func, a, b, N)
    return total

print(rel_error(ex, 0, 1, 0.000000001))

def squared(x):
    '''
    Inputs:
        x - an x value
    
    Output:
        y - the corresponding y value of the function y = x^2
    '''
    y = x**2
    return y

print(rel_error(squared, 0, 10, 0.00000000001))

# Problem #7
real_val = np.e - 1
quad_val = quad(ex, 0, 1)
print(quad_val)

rom_val = romberg(ex, 0, 1)
rom_err = (rom_val - real_val)
print(rom_err)

    

    
