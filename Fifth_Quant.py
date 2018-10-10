import numpy as np
from matplotlib import mlab

delta = 0.0001
epsX, epsE = 0.0005, 0.0005
xmin = 2**(1/6)
npts = 40 #узлы интегрирования
gamma = 21.7
#точки поворота
sxin, sxout = [], []

def potential(x):
    return 4*(x**(-12) - x**(-6))

def integrate(func, min_lim, max_lim, n):
    integral = 0.0
    step = (max_lim - min_lim) / n
    x = min_lim + step
    '''
    for x in mlab.frange(min_lim + step / 2, max_lim - step / 2, step):
        integral += step / 6 * ( func(x - step / 2) + 4 * func(x) + func(x + step / 2) )
    return integral*step/3
    '''
    while (x < max_lim):
        integral += 4 * func(x)
        x = x + step
        integral += 2 * func(x)
        x += step

    integral = step / 3 * (integral + func(min_lim) - func(max_lim))
    return integral

def x1_x2(E):
    return (1 - np.sqrt(1 - E)) / 2, (1 - np.sqrt(1 + E)) / 2

def s_effect(E):
    xmin, xmax = x1_x2(E)
    integral = integrate(potential, xmin, xmax, npts)
    nmax = gamma * integral - 0.5
    return nmax, xmin, xmax

def main_calc(En, ar_xin, ar_xout):
    nmax, xin, xout = s_effect(-0.0005)
    nmax = int(nmax)
    print(nmax)
    #start
    En.append(-0.0005)
    ar_xin.append(xin)
    ar_xout.append(xout)
    E = -1
    de = 0.0001
    for n in range(nmax):
        if n == 100:
            break
        E += de
        de = 2*epsE
        E, xin, xout = calc_energy_level(E, n, de)
        En.append(E)
        ar_xin.append(xin)
        ar_xout.append(xout)

def calc_energy_level(E,n, de):
    s, xin, xout = s_effect(E)
    while np.abs(s - n) >= 0.005:
        E += de
        if E > 0:
            E = -1
            de /= 10.0
        s, xin, xout = s_effect(E)
    return E, xin, xout

nmax, x1min, x2min  = s_effect(-epsE)
En = []
En.append(0)
main_calc(En, sxin, sxout)

kscale = 180/2/gamma

print(En)
