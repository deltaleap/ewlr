from re import L
import numpy as np


import numpy as np

def update(mu, beta, P, R, x, y):
    # compute vector k (with old R and x)
    k = np.einsum('ij,j->j', R, x) / (1 + np.dot(x, np.einsum('ij,j->j', R, x)))
    # update P (with old R, k and x)
    new_P = R - np.einsum('kj,ji->ki', np.einsum('k,j->kj',k,x), R)
    # update R (with old P)
    new_R = P / mu
    e = y - np.dot(beta, x)
    new_beta = beta + e * k

    return new_beta, new_P, new_R


if __name__ == "__main__":
    factor = 0.95
    beta = np.zeros(4)
    P = np.eye(4)
    mu = (1-factor) / factor
    R = P / mu

    x = np.array([10., 12., 14., 3.])
    y = 12.0
    print()
    for i in range(480):
        beta, P, R = update(mu, beta, P, R, x, y)
    print(beta)