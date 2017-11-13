# -*- coding: utf-8 -*-
"""
Morgaine Mandigo-Stoba
Ph20 Problem Set #1
"""

import numpy as np
import matplotlib.pyplot as plt

'''
Problem #2
'''

def comp_x(f_x, a_x, t):
    return (a_x * (np.cos(2 * (np.pi) * f_x * t)))

def comp_y(f_y, a_y, phi, t):
    return (a_y * (np.sin((2 * (np.pi) * f_y * t) + phi)))

def trig_list(f_x, f_y, a_x, a_y, phi, del_t, n):
    x = []
    y = []
    z = []
    for i in range(n):
        t = del_t * (i + 1)
        temp_x = comp_x(f_x, a_x, t)
        temp_y = comp_y(f_y, a_y, phi, t)
        temp_z = (temp_x + temp_y)
        x.append(temp_x)
        y.append(temp_y)
        z.append(temp_z)
    np.savetxt('liss.txt', [x, y, z])
    return [x, y, z]
    
def trig_array(f_x, f_y, a_x, a_y, phi, del_t, n):
    x = np.array([])
    y = np.array([])
    z = np.array([])
    for i in range(n):
        t = del_t * (i + 1)
        temp_x = comp_x(f_x, a_x, t)
        temp_y = comp_y(f_y, a_y, phi, t)
        temp_z = (temp_x + temp_y)
        x.append([temp_x])
        y.append([temp_y])
        z.append([temp_z])
    all_vals = np.array(x, y, z)
    np.savetxt('liss.txt', all_vals)
    return all_vals

'''
Problem 3
'''

# part a
vals = trig_list(2, 4, 1.5, 1.5, 0, 0.005, 200)
plt.plot(vals[0], vals[1])
plt.title('Rational f_x/f_y')
plt.show()

vals = trig_list(1, 4, 2, 3, 0.7, 0.005, 200)
plt.plot(vals[0], vals[1])
plt.title('Rational f_x/f_y')
plt.show()

vals = trig_list(5, 2, 7, 9, 2, 0.005, 200)
plt.plot(vals[0], vals[1])
plt.title('Rational f_x/f_y')
plt.show()


# part b
vals = trig_list(100, 1, 1.5, 1.5, 0, 0.005, 200)
plt.plot(vals[0], vals[1])
plt.title('f_x/f_y >> 1')
plt.show()

vals = trig_list(9, 2, 1.5, 1.5, 0, 0.005, 200)
plt.plot(vals[0], vals[1])
plt.title('f_x/f_y > 1')
plt.show()

vals = trig_list(0, 4, 1.5, 1.5, 0, 0.005, 200)
plt.plot(vals[0], vals[1])
plt.title('f_x/f_y = 0')
plt.show()

vals = trig_list(4, 0, 1.5, 1.5, 0, 0.005, 200)
plt.plot(vals[0], vals[1])
plt.title('f_x/f_y = 0')
plt.show()

vals = trig_list(2, 7, 1.5, 1.5, 0, 0.005, 200)
plt.plot(vals[0], vals[1])
plt.title('f_x/f_y < 1')
plt.show()

vals = trig_list(1, 25, 1.5, 1.5, 0, 0.005, 200)
plt.plot(vals[0], vals[1])
plt.title('f_x/f_y << 1')
plt.show()

# part c
vals = trig_list(1.0, 1.0, 1.0, 1.0, 0, 0.005, 200)
plt.plot(vals[0], vals[1])
plt.title('f_x = f_y = 1, phi = 0')
plt.show()

vals = trig_list(1.0, 1.0, 1.0, 1.0, ((np.pi) / 6.0), 0.005, 200)
plt.plot(vals[0], vals[1])
plt.title('f_x = f_y = 1, phi = pi / 6')
plt.show()

vals = trig_list(1.0, 1.0, 1.0, 1.0, ((np.pi) / 3.0), 0.005, 200)
plt.plot(vals[0], vals[1])
plt.title('f_x = f_y = 1, phi = pi / 3')
plt.show()

vals = trig_list(1.0, 1.0, 1.0, 1.0, ((np.pi) / 2.0), 0.005, 200)
plt.plot(vals[0], vals[1])
plt.title('f_x = f_y = 5, phi = pi / 2')
plt.show()

vals = trig_list(1.0, 1.0, 1.0, 1.0, (2* (np.pi) / 3), 0.005, 200)
plt.plot(vals[0], vals[1])
plt.title('f_x = f_y = 1, phi = 2pi / 3')
plt.show()

vals = trig_list(1.0, 1.0, 1.0, 1.0, (5 * (np.pi) / 6), 0.005, 200)
plt.plot(vals[0], vals[1])
plt.title('f_x = f_y = 7, phi = 5pi / 6')
plt.show()

vals = trig_list(1.0, 1.0, 1.0, 1.0, (np.pi), 0.005, 200)
plt.plot(vals[0], vals[1])
plt.title('f_x = f_y = 7, phi = pi')
plt.show()

vals = trig_list(1.0, 1.0, 1.0, 1.0, (7*(np.pi) / 6), 0.005, 200)
plt.plot(vals[0], vals[1])
plt.title('f_x = f_y = 7, phi = 7pi / 6')
plt.show()

vals = trig_list(1.0, 1.0, 1.0, 1.0, (4 * (np.pi) / 3), 0.005, 200)
plt.plot(vals[0], vals[1])
plt.title('f_x = f_y = 7, phi = 4pi / 3')
plt.show()

vals = trig_list(1.0, 1.0, 1.0, 1.0, (3 * (np.pi) / 2), 0.005, 200)
plt.plot(vals[0], vals[1])
plt.title('f_x = f_y = 7, phi = 3pi / 2')
plt.show()

vals = trig_list(1.0, 1.0, 1.0, 1.0, (5 * (np.pi) / 3), 0.005, 200)
plt.plot(vals[0], vals[1])
plt.title('f_x = f_y = 7, phi = 5pi / 3')
plt.show()

vals = trig_list(1.0, 1.0, 1.0, 1.0, (11 * (np.pi) / 6), 0.005, 200)
plt.plot(vals[0], vals[1])
plt.title('f_x = f_y = 7, phi = 11 pi / 6')
plt.show()

vals = trig_list(1.0, 1.0, 1.0, 1.0, (2 * (np.pi)), 0.005, 200)
plt.plot(vals[0], vals[1])
plt.title('f_x = f_y = 7, phi = 2 pi')
plt.show()

'''
Problem 4
'''
ts = []
for i in range(500):
    temp_t = 0.005 * (i + 1)
    ts.append(temp_t)
vals = trig_list(4.0, 5.0, 1.0, 1.0, 0.0, 0.005, 500)
plt.plot(ts, vals[2])
plt.title('Beats with f_x = 4, f_y = 5')
plt.show()

vals = trig_list(20.0, 21.0, 1.0, 1.0, 0.0, 0.005, 500)
plt.plot(ts, vals[2])
plt.title('Beats with f_x = 20, f_y = 21')
plt.show()

vals = trig_list(100.0, 110.0, 1.0, 1.0, 0.0, 0.005, 500)
plt.plot(ts, vals[2])
plt.title('Beats with f_x = 100, f_y = 110')
plt.show()

    
    
