# --- external imports
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from math import gcd
import matplotlib.ticker as ticker
from fractions import Fraction
from matplotlib.ticker import MaxNLocator
# --- internal imports
import functions.arguments as fa
from models.hofstadter import Hofstadter

# plt.rc('text', usetex=True)
# plt.rc('text.latex', preamble=r'\usepackage{amsmath}')


def hamiltonian(t, p, M, k_val):

    # initialize the Hamiltonian
    Hamiltonian = np.zeros((M, M), dtype=np.complex128)

    nphi = p / q

    def B(nphi_val, k_val_val, m_val):
        return (- 2 * t * np.exp(1j * np.pi * nphi_val/3)
                * np.cos(np.pi * nphi_val * (m_val + 1/2) + 3*k_val_val[1]*np.sqrt(3)/2))

    def C(nphi_val):
        return -t * np.exp(1j * np.pi * nphi_val / 3)

    for i in range(M - 1):
        Hamiltonian[i][i + 1] = B(nphi, k_val, i+1)
        Hamiltonian[i + 1][i] = np.conj(B(nphi, k_val, i+1))
    for i in range(M - 2):
        Hamiltonian[i][i + 2] = np.conj(C(nphi))
        Hamiltonian[i + 2][i] = C(nphi)

    # boundary terms
    Hamiltonian[0][M - 1] = np.conj(B(nphi, k_val, M))
    Hamiltonian[0][M - 2] = C(nphi)
    Hamiltonian[1][M - 1] = C(nphi)

    Hamiltonian[M - 1][0] = B(nphi, k_val, M)
    Hamiltonian[M - 2][0] = np.conj(C(nphi))
    Hamiltonian[M - 1][1] = np.conj(C(nphi))

    lmbda = np.real(np.linalg.eigvals(Hamiltonian))
    eenergies = np.zeros(2 * len(lmbda))
    for i in range(len(lmbda)):
        eenergies[i] = +np.sqrt(3 + lmbda[i])
        eenergies[len(lmbda) + i] = -np.sqrt(3 + lmbda[i])

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

        if p % 2 == 0:
            M = q
        else:
            M = 2 * q

        nphi = p / q
        nphi_list = [nphi for i in range(2*M)]
        eham = hamiltonian(-1, p, M, np.array([0, 0]))
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
