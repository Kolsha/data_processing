import numpy as np
from matplotlib import mlab

delta = 0.0001
epsX, epsE = 0.0005, 0.0005
xmin = 2**(1/6)
npts = 40 #узлы интегрирования
#gamma = 600#21.7
gamma = 1/np.pi
#точки поворота
sxin, sxout = [], []

def potential(x):
    #return 4*(x**(-12) - x**(-6))
    return x*x/4
'''
def simpson_method(func, mim_lim, max_lim, delta):
    def integrate(func, mim_lim, max_lim, n):
        integral = 0.0
        step = (max_lim - mim_lim) / n
        for x in mlab.frange(mim_lim + step / 2, max_lim - step / 2, step):
            integral += step / 6 * ( func(x - step / 2) + 4 * func(x) + func(x + step / 2) )
        return integral

    d, n = 1, 2
    a1, a2 = 0, 0
    while np.fabs(d) > delta:
        a1 = integrate(func, mim_lim, max_lim, (n + 2))
        a2 = integrate(func, mim_lim, max_lim, n)
        d = (a1 - a2) / 15
        n += 2

    return a1

def x1_x2(E):
    return (1 - np.sqrt(1 - E)) / 2, (1 - np.sqrt(1 + E)) / 2

x1min, x2min = x1_x2(-epsE)
nmax = gamma * simpson_method(potential, x1min, x2min, delta) - 0.5
print(nmax)
'''

def s_effect(E):
    #ищем пределы интеграла
    xin = xmin
    dx = 0.1
    while dx > epsX:
        xin -= dx
        if potential(xin) < E:
            break
        xin += dx
        dx /= 2

    xout = xmin
    dx = 0.1
    while dx > epsX:
        xout += dx
        if potential(xout) < E:
            break
        xout -= dx
        dx /= 2

    H = (xout - xin)/npts
    #simpson
    '''
    def my_func(x):
        return np.sqrt(E - potential((x)))

    integral = simpson_method(my_func, xin - H, xin + H, delta)
    '''
    integral = np.sqrt(np.abs(E - potential(xin - H)))
    fac = 2
    for i in range(2, npts - 2):
        x = xin + i * H
        if fac == 2:
            fac = 4
        else:
            fac = 2
        integral += fac*np.sqrt(np.abs(E-potential(x)))
    integral += np.sqrt(np.abs(E - potential(xout - H)))
    integral *= H/3
    #учет точек поворота
    integral += np.sqrt(np.abs(E - potential(xin + H))) * 2 * H/3
    integral += np.sqrt(np.abs(E - potential(xout - H))) * 2 * H / 3
    s = gamma*integral
    return s, int(s - 0.5), xin, xout

def main_calc(En, ar_xin, ar_xout):
    s, nmax, xin, xout = s_effect(-0.0005)
    print(nmax)
    #start
    #En.append(-0.0005)
    #ar_xin.append(xin)
    #ar_xout.append(xout)
    E1 = -1
    F1 = - np.pi / 2
    for n in range(nmax):
        if n == 4:
            break
        E2 = E1 + np.abs(E1)/4
        de = 2*epsE
        F2, E2, xin, xout = calc_energy_level(de, E1, E2, F1, n, s, ar_xin, ar_xout)
        En.append(E2)
        ar_xin.append(xin)
        ar_xout.append(xout)
        E1 = E2
        F1 = F2 - np.pi

def calc_energy_level(de, E1, E2, F1, n, s, ar_xin, ar_xout):
    while np.abs(de) >= epsE:
        s, _, xin, xout = s_effect(E2)
        F2 = s - (n + 0.5)*np.pi
        if F2 == F1:
            return F2, E2, xin, xout
        de = -F2*(E2-E1)/(F2-F1)
        E1 = E2
        F1 = F2
        E2 = E1 + de
        if E2 > 0:
            E2 = -epsE
    return F2, E2, xin, xout
En = []
main_calc(En, sxin, sxout)

kscale = 180/2/gamma

print(En)
