import numpy as np


class Hofstadter:

    def __init__(self, p, q, a0=1, t=1):
        self.p = p
        self.q = q
        self.a0 = a0
        self.t = t

    def unit_cell(self):

        num_bands_val = self.q
        # lattice vectors
        a1 = self.a0 * np.array([num_bands_val, 0])
        a2 = self.a0 * np.array([0, 1])
        avec_val = np.vstack((a1, a2))
        # reciprocal lattice vectors
        b1 = (2. * np.pi) / self.a0 * np.array([1 / num_bands_val, 0])
        b2 = (2. * np.pi) / self.a0 * np.array([0, 1])
        bvec_val = np.vstack((b1, b2))
        # symmetry points
        GA = np.array([0, 0])
        Y = np.array([0, 0.5])
        S = np.array([0.5, 0.5])
        X = np.array([0.5, 0])
        sym_points_val = [GA, Y, S, X]

        return num_bands_val, avec_val, bvec_val, sym_points_val

    def hamiltonian(self, k_val):

        # initialize the Hamiltonian
        num_bands_val = self.q
        Hamiltonian = np.zeros((num_bands_val, num_bands_val), dtype=np.complex128)

        # nearest neighbors
        delta = np.zeros((2, 2))
        delta[0, :] = self.a0 * np.array([1, 0])
        delta[1, :] = self.a0 * np.array([0, 1])

        nphi = self.p / self.q

        def h(k_val_val, m_val):
            return 2 * np.cos(2 * np.pi * nphi * m_val + k_val_val[1] * self.a0)

        for n in range(self.q):
            Hamiltonian[n][n] = self.t * h(k_val, n)

        for n in range(self.q - 1):
            Hamiltonian[n][n + 1] = self.t * np.exp(+1j * k_val[0] * self.a0)
            Hamiltonian[n + 1][n] = self.t * np.exp(-1j * k_val[0] * self.a0)

        Hamiltonian[0][self.q - 1] = self.t * np.exp(-1j * k_val[0] * self.a0)
        Hamiltonian[self.q - 1][0] = self.t * np.exp(+1j * k_val[0] * self.a0)

        return Hamiltonian
