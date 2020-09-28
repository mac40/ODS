'''
Standard Frank-Wolfe algorithm
'''

import sys
import time

import numpy as np


def FW(Q,  # noqa
       c,
       x,
       maxit,
       maxtime,
       eps,
       fstop,
       stopcr):
    '''
    implementation of standard Frank-Wolfe Algorithm
    '''
    gamma = 1e-4
    flagls = 0

    n = len(Q[0])

    fh = 0
    timeVec = 0

    if (maxit < np.inf):
        fh = np.zeros((1, maxit))
        timeVec = np.zeros((1, maxit))
    else:
        fh = np.zeros((1, 100*n))
        timeVec = np.zeros((1, 100*n))

    it = 0

    tstart = time.process_time()

    while (it < maxit and flagls == 0):
        Qx = Q.dot(x)
        xQx = x.T.dot(Qx)
        cx = c.T.dot(x)

        print("------ITERATION {}------".format(it))
        print("x\n", x)
        print("Qx", Qx)
        print("xQx", xQx)
        print("cx", cx)

        if (it == 0):
            fx = 0.5 * xQx - cx
            timeVec[0, it] = 0
        else:
            tstop = time.process_time()
            timeVec[0, it] = tstop - tstart

        fh[0, it] = fx

        # Gradient evaluation
        g = Qx - c

        # Time break
        if (timeVec[0, it] > maxtime):
            break

        # solution of FW problem
        istar = np.where(g == np.amin(g))[0][0]
        xstar = np.zeros((n, 1))
        xstar[istar] = 1

        print("g", g)
        print("g[istar]", g[istar])
        print("xstar", xstar)

        # direction calculation
        d = xstar - x
        gnr = g.T.dot(d)

        print("d", d)
        print("gnr", gnr)

        if (stopcr == 1):
            if (fx <= fstop):
                break
        elif (stopcr == 2):
            if (gnr >= -eps):
                break
        else:
            print("Unknown stopping criterion")
            sys.exit()

        # Armijo search
        alpha = 1
        ref = gamma * gnr

        while True:
            fz = 0.5 * ((1-alpha)**2 * xQx + 2 * alpha * (1 - alpha) *
                        Qx[istar] + alpha**2 * Q[istar, istar]) - (
                            (1 - alpha) * cx + alpha * c[istar])
            if (fz <= (fx + alpha * ref)):
                z = x + alpha*d
                break
            else:
                alpha = alpha * 0.5

            if (alpha <= 1e-20):
                z = x
                fz = fx
                flagls = 1
                it = it - 1
                break

        x = z
        fx = fz

        it = it + 1

    ttot = time.process_time() - tstart

    if (it < len(fh[0])):
        fh = fh[:, 0:it]
        timeVec = timeVec[:, 0:it]

    return x, it, fx[0][0], ttot, fh, timeVec
