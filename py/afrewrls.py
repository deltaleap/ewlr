from re import L
import numpy as np


import numpy as np

def update(alpha, l, l_min, l_max, P, S, beta, phi, x, y):

    I = np.eye(4, 4)
    e = y - np.dot(beta, x)
    print(f"{e=}")
    # udpate lambda(t)
    new_lambda = l + alpha * np.dot(phi, x) * e
    # truncate lambda(t)
    if new_lambda > 1.0:
        new_lambda = 1.0
    if new_lambda < 0.7:
        new_lambda = 0.7

    # update k(t)
    k = (np.dot(P, x) / new_lambda) / (1 - np.dot(x, np.dot(P, x)) / new_lambda)
    print(f"{k=}")

    # update beta(t)
    new_beta = beta + e * k

    # udpate P(t)
    new_P = P / new_lambda - np.einsum(
        "ij,jk->ik", np.einsum("i,j -> ij", k, x), P
    ) / new_lambda

    # update S(t)
    I_minus_kx = I - np.einsum("i,k->ik", k, x)
    new_S = ((
        np.einsum("ij,jk->ik", I_minus_kx, np.einsum("ij,jk->ik", S, I_minus_kx))
    ) - new_P + np.einsum("i,k->ik", k, k)) / new_lambda

    new_phi = np.einsum('j,ij->i', phi, I_minus_kx) + np.einsum('j,ji->i', x, new_S) * e
    print(f"phi part 1: {np.einsum('j,ij->i', phi, I_minus_kx)}")
    print(f"phi part 2: {np.dot(new_S, x) * e}")
    print(new_S)
    print(x)
    print(f"phi part 2.1: {np.dot(new_S, x)}")
    return new_lambda, new_P, new_S, new_beta, new_phi

if __name__ == "__main__":
    alpha = 0.95
    l_min = 0.7
    l_max = 1.0
    P = np.eye(4, 4)
    S = np.eye(4, 4)
    beta = np.zeros(4)
    phi = np.zeros(4)
    l = 0.0

    x = np.array([10., 12., 14., 3.])
    y = 12.0
    print("step 1:")
    l, P, S, beta, phi = update(alpha, l, l_min, l_max, P, S, beta, phi, x, y)
    print(f"{alpha=}")
    print(f"{l=}")
    print(f"{beta=}")
    print(f"{P=}")
    print(f"{S=}")
    print(f"{phi=}")

    print("step 2:")
    x = np.array([10.1, 11.8, 9., 3.5])
    y = 12.0
    l, P, S, beta, phi = update(alpha, l, l_min, l_max, P, S, beta, phi, x, y)
    print(f"{alpha=}")
    print(f"{l=}")
    print(f"{P=}")
    print(f"{S=}")
    print(f"{beta=}")
    print(f"{phi=}")

    print("step 3:")
    x = np.array([10.1, 11.8, 9., 3.5])
    y = 11.0
    l, P, S, beta, phi = update(alpha, l, l_min, l_max, P, S, beta, phi, x, y)
    print(f"{alpha=}")
    print(f"{l=}")
    print(f"{P=}")
    print(f"{S=}")
    print(f"{beta=}")
    print(f"{phi=}")

    for i in range(100):
        x = np.array([1., 1., 1., 1.])
        y = 11.0
        l, P, S, beta, phi = update(alpha, l, l_min, l_max, P, S, beta, phi, x, y)
    
    print(beta)
    print(np.dot(beta, x))