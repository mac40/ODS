'''
Main function for Frank-Wolfe algorithms
'''

import matplotlib.pyplot as plt
import numpy as np

import Algorithms.FrankWolfe as FrankWolfe  # noqa


def main():
    '''
    main function
    '''

    # number of samples
    m = 2**5

    # number of variables
    # n = 2**12
    n = 3

    maxit = 3000
    maxtime = 100

    fstop = np.zeros((1, 1))

    x0 = np.zeros((n, 1))
    x0[0] = 1

    Q = np.random.randn(m, n)
    print("Q\n", Q)
    c = np.sum(np.square(Q), axis=0).T
    print("c\n", c)
    c = np.reshape(c, (n, 1))
    Q = 2*(Q.T.dot(Q))
    print("Q\n", Q)

    eps = 1e-6
    stopcr = 2

    result = FrankWolfe.FW(Q,
                           c,
                           x0,
                           maxit,
                           maxtime,
                           eps,
                           fstop[0],
                           stopcr)

    print("*****************\n*  FW STANDARD  *\n*****************\n")
    print("0.5 * xQx - cx = {}\n".format(result[2]))
    print("Number of non-zero component of x = {}\n".format(
        np.count_nonzero(result[0])))
    print("Number of iterations = {}\n".format(result[1]))
    print("Cpu time = {}\n".format(result[3]))

    plt.plot(result[5][0], result[4][0])
    plt.show()


if __name__ == "__main__":
    main()
