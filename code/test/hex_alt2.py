# --- external imports
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from math import gcd
import matplotlib.ticker as ticker
from fractions import Fraction
from matplotlib.ticker import MaxNLocator
from scipy import linalg
# --- internal imports
import functions.arguments as fa
from models.hofstadter import Hofstadter


def polyeig(*A):
    """
    Solve the polynomial eigenvalue problem:
        (A0 + e A1 +...+  e**p Ap)x=0 

    Return the eigenvectors [x_i] and eigenvalues [e_i] that are solutions.

    Usage:
        X,e = polyeig(A0,A1,..,Ap)

    Most common usage, to solve a second order system: (K + C e + M e**2) x =0
        X,e = polyeig(K,C,M)

    """
    if len(A)<=0:
        raise Exception('Provide at least one matrix')
    for Ai in A:
        if Ai.shape[0] != Ai.shape[1]:
            raise Exception('Matrices must be square')
        if Ai.shape != A[0].shape:
            raise Exception('All matrices must have the same shapes');

    n = A[0].shape[0]
    l = len(A)-1
    # Assemble matrices for generalized problem
    C = np.block([
        [np.zeros((n*(l-1),n)), np.eye(n*(l-1))],
        [-np.column_stack( A[0:-1])]
        ])
    D = np.block([
        [np.eye(n*(l-1)), np.zeros((n*(l-1), n))],
        [np.zeros((n, n*(l-1))), A[-1]          ]
        ]);
    # Solve generalized eigenvalue problem
    e, X = linalg.eig(C, D);
    if np.all(np.isreal(e)):
        e=np.real(e)
    X=X[:n,:]

    # Sort eigenvalues/vectors
    #I = np.argsort(e)
    #X = X[:,I]
    #e = e[I]

    # Scaling each mode by max
    X /= np.tile(np.max(np.abs(X),axis=0), (n,1))

    return X, e


def hamiltonian(p, M, k_val):

    # initialize the Hamiltonian
    Hamiltonian = np.zeros((M, M), dtype=np.complex128)

    nphi = p / q

    def A(nphi_val, m_val):
        value = 2 * np.cos(2 * np.pi * nphi_val * (m_val + 1/6)) + 3
        return value

    def B(nphi_val, m_val):
        value = 2 * np.cos(np.pi * nphi_val * (m_val + 1/6)) * np.exp(- 1j * np.pi * nphi * 0.5 * (2 * m_val + 1))
        return value

    upper_diag_array = np.array([A(nphi, m) for m in range(q)])
    Hamiltonian += np.roll(np.diag(upper_diag_array), 0, axis=1)

    upper_diag_array2 = np.array([B(nphi, m) for m in range(q)])
    Hamiltonian += np.roll(np.diag(upper_diag_array2), 1, axis=1)
    Hamiltonian += np.roll(np.diag(np.conj(upper_diag_array2)), 1, axis=0)

    # lmbda = np.real(np.linalg.eigvals(Hamiltonian))
    # eenergies = np.zeros(2 * len(lmbda))
    # for i in range(len(lmbda)):
    #     eenergies[i] = +np.sqrt(lmbda[i])
    #     eenergies[len(lmbda) + i] = -np.sqrt(lmbda[i])

    _, eenergies = polyeig(-Hamiltonian, np.zeros((q, q)), np.eye(q))

    return eenergies


if __name__ == '__main__':

    # input arguments
    q = 199

    # construct butterfly
    nphi_list, E_list = [], []
    values = []

    for p in range(1, q):

        if gcd(p, q) != 1:  # nphi must be a coprime fraction
            continue

        M = q

        nphi = p / q
        nphi_list = [nphi for i in range(2*M)]
        eham = hamiltonian(p, M, np.array([0, 0]))
        values.append((eham, nphi_list))

    # construct figure
    fig = plt.figure()
    ax = plt.subplot(111)
    ax.set_title(f"$n_\phi = p/{q}$")

    # nphi_list = list(np.concatenate(nphi_list).ravel())
    # E_list = list(np.concatenate(E_list).ravel())
    #
    # ax.scatter(nphi_list, E_list, s=1, marker='.')

    for eigenvalues, alphas in values:
        ax.plot(alphas, eigenvalues, '.', color='r', markersize=0.5)

    ax.set_ylabel('$E$')
    ax.set_xlabel('$n_\phi$')
    plt.show()
