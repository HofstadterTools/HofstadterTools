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

factor = 2

def AA_block_func(p, q):

    # initialize the Hamiltonian
    Hamiltonian = np.zeros((q, q), dtype=np.complex128)

    nphi = p / q

    def A(nphi_val, m_val):
        value = -2*np.cos(2 * np.pi * factor * nphi_val * m_val)
        return value

    def B_plus(nphi_val, m_val):
        value = -np.exp(-1j * np.pi * factor * nphi_val * (m_val + 1/2)) - np.exp(1j * np.pi * factor * nphi_val * (m_val + 1/2))
        return value

    # def B_minus(nphi_val, m_val):
    #     value = -np.exp(-1j * np.pi * nphi_val * (m_val - 1/2)) -np.exp(1j * np.pi * nphi_val * (m_val - 1/2))
    #     return value

    diag_array = np.array([A(nphi, m) for m in range(q)])
    Hamiltonian += np.roll(np.diag(diag_array), 0, axis=1)

    upper_diag_array = np.array([B_plus(nphi, m) for m in range(q)])
    Hamiltonian += np.roll(np.diag(upper_diag_array), 1, axis=1)
    Hamiltonian += np.roll(np.diag(np.conj(upper_diag_array)), 1, axis=0)

    # lower_diag_array = np.array([B_minus(nphi, m) for m in range(q)])
    # Hamiltonian += np.roll(np.diag(lower_diag_array), 1, axis=0)
    # Hamiltonian += np.roll(np.diag(np.conj(lower_diag_array)), 1, axis=1)

    return Hamiltonian


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
        nphi_list = [nphi for i in range(M)]

        AA_block = AA_block_func(p, q)
        eham = np.linalg.eigvalsh(AA_block)
        values.append((eham, nphi_list))

    # construct figure
    fig = plt.figure()
    ax = plt.subplot(111)
    ax.set_title(f"$n_\phi = p/{q}$")
    for eigenvalues, alphas in values:
        ax.plot(alphas, eigenvalues, '.', color='r', markersize=0.5)
    ax.set_ylabel('$E$')
    ax.set_xlabel('$n_\phi$')
    plt.show()
