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

factor = 1

def AA_block_func(p, q):

    # initialize the Hamiltonian
    Hamiltonian = np.zeros((q, q), dtype=np.complex128)

    nphi = p / q

    def A(nphi_val, m_val):
        value = -np.exp(-1j * 2 * np.pi * factor * nphi_val * m_val)
        return value

    def B_plus(nphi_val, m_val):
        value = -np.exp(-1j * np.pi * factor * nphi_val * (m_val + 1/2))
        return value

    diag_array = np.array([A(nphi, m) for m in range(q)])
    Hamiltonian += np.roll(np.diag(diag_array), 0, axis=1)

    upper_diag_array = np.array([B_plus(nphi, m) for m in range(q)])
    Hamiltonian += np.roll(np.diag(upper_diag_array), 1, axis=1)

    return Hamiltonian + np.conj(Hamiltonian).T


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
        nphi_list = [nphi for i in range(M-1)]

        AA_block = AA_block_func(p, q)
        eham = np.linalg.eigvalsh(AA_block)
        # int_dos = [np.sum(eham[:i+1]) for i in range(len(eham)-1)]
        gaps = [eham[i+1] - eham[i] for i in range(len(eham)-1)]
        scaled_gaps = [i/q for i in range(len(eham)-1)]
        # values.append((nphi_list, eham))
        values.append((nphi_list, scaled_gaps, gaps))

    # construct figure
    fig = plt.figure()
    ax = plt.subplot(111)
    ax.set_title(f"$n_\phi = p/{q}$")
    for alphas, eigenvalues, gaps in values:
        ax.scatter(alphas, eigenvalues, s=[5*i for i in gaps], c='r', linewidths=0)
    ax.set_ylabel('$E$')
    ax.set_xlabel('$n_\phi$')
    plt.show()
