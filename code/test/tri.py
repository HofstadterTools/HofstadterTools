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


def hamiltonian(t, p, q, k_val):

    M = q

    # initialize the Hamiltonian
    Hamiltonian = np.zeros((M, M), dtype=np.complex128)

    nphi = p / q

    def B(k_val_val, m_val):
        return - 2 * t * np.cos(2*np.pi * nphi * (m_val + 1/2) + k_val_val[1]*np.sqrt(3)/2)

    for i in range(M - 1):
        Hamiltonian[i][i + 1] = B(k_val, i+1)
        Hamiltonian[i + 1][i] = B(k_val, i+1)

    for i in range(M - 2):
        Hamiltonian[i][i + 2] = -t * np.exp(+1j * k_val[0])
        Hamiltonian[i + 2][i] = -t * np.exp(-1j * k_val[0])

    # boundary terms
    Hamiltonian[0][M - 1] = B(k_val, M)
    Hamiltonian[M - 1][0] = B(k_val, M)

    Hamiltonian[0][M - 2] = -t * np.exp(-1j * k_val[0])
    Hamiltonian[1][M - 1] = -t * np.exp(-1j * k_val[0])
    Hamiltonian[M - 2][0] = -t * np.exp(+1j * k_val[0])
    Hamiltonian[M - 1][1] = -t * np.exp(+1j * k_val[0])

    return Hamiltonian


if __name__ == '__main__':

    # input arguments
    q = 299

    # construct butterfly
    nphi_list, E_list = [], []

    for p in range(1, q):

        if gcd(p, q) != 1:  # nphi must be a coprime fraction
            continue
        nphi = p / q

        ham = hamiltonian(1, p, q, np.array([0, 0]))

        M = q

        nphi_list.append([nphi] * M)
        E_list.append(np.sort(np.linalg.eigvalsh(ham)))

    # construct figure
    fig = plt.figure()
    ax = plt.subplot(111)
    ax.set_title(f"$n_\phi = p/{q}$")

    nphi_list = list(np.concatenate(nphi_list).ravel())
    E_list = list(np.concatenate(E_list).ravel())

    ax.scatter(nphi_list, E_list, s=1, marker='.')
    ax.set_ylabel('$E$')
    ax.set_xlabel('$n_\phi$')
    plt.show()
