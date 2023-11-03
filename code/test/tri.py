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


def hamiltonian(t, p, q, k_val):

    M = q

    # initialize the Hamiltonian
    Hamiltonian = np.zeros((M, M), dtype=np.complex128)

    nphi = p / q

    def B(k_val_val, m_val):
        return - 2 * np.cos(np.pi * nphi * (m_val + 1/2)) * np.exp(-1j * np.pi * nphi * (2*m_val + 1) / 2)

    # def B_lower(k_val_val, m_val):
    #     return - 2 * np.cos(np.pi * nphi * (m_val - 1/2))

    def C(k_val_val, m_val):
        return - np.exp(+1j * k_val_val[0]) * np.exp(-1j * np.pi * nphi * (4*m_val + 4) / 2)

    upper_diag_array = np.array([B(k_val, m) for m in range(q)])
    Hamiltonian += np.roll(np.diag(upper_diag_array), 1, axis=1)
    # lower_diag_array = np.array([B_lower(k_val, m+1) for m in range(q)])
    Hamiltonian += np.roll(np.diag(np.conj(upper_diag_array)), 1, axis=0)

    upper_diag_array2 = np.array([C(k_val, m) for m in range(q)])
    Hamiltonian += np.roll(np.diag(upper_diag_array2), 2, axis=1)
    Hamiltonian += np.roll(np.diag(np.conj(upper_diag_array2)), 2, axis=0)

    return Hamiltonian


if __name__ == '__main__':

    # input arguments
    q = 199

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
