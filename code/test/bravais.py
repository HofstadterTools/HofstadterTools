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


def hamiltonian(p, q, k_val):

    # t0, t1, t2, t3 = 0.0345, 0.0345, 0.0345, 0
    t0, t1, t2, t3 = 0.5, 0.5, 0.2, 0.2
    M = q
    theta = (67/180) * np.pi

    # initialize the Hamiltonian
    Hamiltonian = np.zeros((M, M), dtype=np.complex128)

    nphi = p / q

    def A(k_val_val, m_val):
        return - 2 * t0 * np.cos(2*np.pi*nphi*m_val + k_val_val[0])

    def B(k_val_val, m_val, theta_val):
        return (-t1 * np.exp(1j*2*np.pi*nphi*np.cos(theta_val)*(m_val+1/2))
                -t2 * np.exp(-1j*2*np.pi*nphi*(1-np.cos(theta_val))*(m_val+1/2)-1j*k_val_val[0])
                -t3 * np.exp(1j*2*np.pi*nphi*(1+np.cos(theta_val))*(m_val+1/2)+1j*k_val_val[0]))

    diag_array = np.array([A(k_val, m) for m in range(q)])
    Hamiltonian += np.roll(np.diag(diag_array), 0, axis=1)

    upper_diag_array = np.array([B(k_val, m, theta) for m in range(q)])
    Hamiltonian += np.roll(np.diag(upper_diag_array), 1, axis=1)
    Hamiltonian += np.roll(np.diag(np.conj(upper_diag_array)), 1, axis=0)

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

        ham = hamiltonian(p, q, np.array([0, 0]))

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
