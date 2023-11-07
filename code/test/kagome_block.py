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


factor = 8


def AB_block_func(p, q):

    # initialize the Hamiltonian
    Hamiltonian = np.zeros((q, q), dtype=np.complex128)

    nphi = p / q

    def A(nphi_val, m_val):
        value = -2 * np.cos(np.pi * factor * nphi_val * m_val)
        return value

    upper_diag_array = np.array([A(nphi, m) for m in range(q)])
    Hamiltonian += np.roll(np.diag(upper_diag_array), 0, axis=1)

    return Hamiltonian


def AC_block_func(p, q):

    # initialize the Hamiltonian
    Hamiltonian = np.zeros((q, q), dtype=np.complex128)

    nphi = p / q

    def A(nphi_val, m_val):
        value = -np.exp(- 1j * np.pi * factor * nphi_val * 0.5 * (m_val + 1/4))
        return value

    def B_minus(nphi_val, m_val):
        value = -np.exp(+ 1j * np.pi * factor * nphi_val * 0.5 * (m_val - 1/4))
        return value

    upper_diag_array = np.array([A(nphi, m) for m in range(q)])
    Hamiltonian += np.roll(np.diag(upper_diag_array), 0, axis=1)

    upper_diag_array2 = np.array([B_minus(nphi, m) for m in range(q)])
    Hamiltonian += np.roll(np.diag(upper_diag_array2), 1, axis=0)

    return Hamiltonian


def BA_block_func(p, q):

    # initialize the Hamiltonian
    Hamiltonian = np.zeros((q, q), dtype=np.complex128)

    nphi = p / q

    def A(nphi_val, m_val):
        value = -2 * np.cos(np.pi * factor * nphi_val * m_val)
        return value

    upper_diag_array = np.array([A(nphi, m) for m in range(q)])
    Hamiltonian += np.roll(np.diag(upper_diag_array), 0, axis=1)

    return Hamiltonian


def BC_block_func(p, q):

    # initialize the Hamiltonian
    Hamiltonian = np.zeros((q, q), dtype=np.complex128)

    nphi = p / q

    def A(nphi_val, m_val):
        value = -np.exp(+ 1j * np.pi * factor * nphi_val * 0.5 * (m_val + 1/4))
        return value

    def B_minus(nphi_val, m_val):
        value = -np.exp(- 1j * np.pi * factor * nphi_val * 0.5 * (m_val - 1/4))
        return value

    upper_diag_array = np.array([A(nphi, m) for m in range(q)])
    Hamiltonian += np.roll(np.diag(upper_diag_array), 0, axis=1)

    upper_diag_array2 = np.array([B_minus(nphi, m) for m in range(q)])
    Hamiltonian += np.roll(np.diag(upper_diag_array2), 1, axis=0)

    return Hamiltonian


def CA_block_func(p, q):

    # initialize the Hamiltonian
    Hamiltonian = np.zeros((q, q), dtype=np.complex128)

    nphi = p / q

    def A(nphi_val, m_val):
        value = -np.exp(+1j * np.pi * factor * nphi_val * 0.5 * (m_val+1/4))
        return value

    def B_plus(nphi_val, m_val):
        value = -np.exp(-1j * np.pi * factor * nphi_val * 0.5 * (m_val+3/4))
        return value

    upper_diag_array = np.array([A(nphi, m) for m in range(q)])
    Hamiltonian += np.roll(np.diag(upper_diag_array), 0, axis=1)

    upper_diag_array2 = np.array([B_plus(nphi, m) for m in range(q)])
    Hamiltonian += np.roll(np.diag(upper_diag_array2), 1, axis=1)

    return Hamiltonian


def CB_block_func(p, q):

    # initialize the Hamiltonian
    Hamiltonian = np.zeros((q, q), dtype=np.complex128)

    nphi = p / q

    def A(nphi_val, m_val):
        value = -np.exp(-1j * np.pi * factor * nphi_val * 0.5 * (m_val + 1/4))
        return value

    def B_plus(nphi_val, m_val):
        value = -np.exp(+1j * np.pi * factor * nphi_val * 0.5 * (m_val + 3/4))
        return value

    upper_diag_array = np.array([A(nphi, m) for m in range(q)])
    Hamiltonian += np.roll(np.diag(upper_diag_array), 0, axis=1)

    upper_diag_array2 = np.array([B_plus(nphi, m) for m in range(q)])
    Hamiltonian += np.roll(np.diag(upper_diag_array2), 1, axis=1)

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
        M = 3*q

        nphi = p / q
        nphi_list = [nphi for i in range(M)]

        AA_block = np.zeros((q, q))
        AB_block = AB_block_func(p, q)
        AC_block = AC_block_func(p, q)
        #
        BA_block = BA_block_func(p, q)
        BB_block = np.zeros((q, q))
        BC_block = BC_block_func(p, q)
        #
        CA_block = CA_block_func(p, q)
        CB_block = CB_block_func(p, q)
        CC_block = np.zeros((q, q))

        upper = np.concatenate((AA_block, AB_block, AC_block), axis=1)
        middle = np.concatenate((BA_block, BB_block, BC_block), axis=1)
        lower = np.concatenate((CA_block, CB_block, CC_block), axis=1)
        matrix = np.concatenate((upper, middle, lower), axis=0)

        eham = np.linalg.eigvalsh(matrix)
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
