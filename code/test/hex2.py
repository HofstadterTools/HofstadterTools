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
        value = - 2 * t * np.exp(1j * 2 * np.pi * nphi_val/3) * np.cos(2 * np.pi * nphi_val * (m_val + 1/2) + 3*k_val_val[1]*np.sqrt(3)/2)
        return value

    def C(nphi_val):
        value = -t * np.exp(-1j * 2 * np.pi * nphi_val / 3)
        return value

    upper_diag_array = np.array([B(nphi, k_val, (m + 1) % q) for m in range(q)])
    Hamiltonian += np.roll(np.diag(upper_diag_array), 1, axis=1)
    upper_diag_array2 = np.array([C(nphi) for m in range(q)])
    Hamiltonian += np.roll(np.diag(upper_diag_array2), 2, axis=1)

    lower_diag_array = np.array([np.conj(B(nphi, k_val, (m + 1) % q)) for m in range(q)])
    Hamiltonian += np.roll(np.diag(lower_diag_array), 1, axis=0)
    lower_diag_array2 = np.array([np.conj(C(nphi)) for m in range(q)])
    Hamiltonian += np.roll(np.diag(lower_diag_array2), 2, axis=0)

    print(Hamiltonian)
    1/0

    lmbda = np.real(np.linalg.eigvals(Hamiltonian))
    eenergies = np.zeros(2 * len(lmbda))
    for i in range(len(lmbda)):
        eenergies[i] = +np.sqrt(3 + lmbda[i])
        eenergies[len(lmbda) + i] = -np.sqrt(3 + lmbda[i])

    return eenergies


if __name__ == '__main__':

    # input arguments
    q = 5

    # construct butterfly
    nphi_list, E_list = [], []
    values = []

    for p in range(1, q):

        if gcd(p, q) != 1:  # nphi must be a coprime fraction
            continue

        # if p % 2 == 0:
        #     M = q
        # else:
        #     M = 2 * q

        M = q

        nphi = p / q
        #print(nphi)
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
