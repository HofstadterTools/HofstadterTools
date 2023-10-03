"""The Hofstadter model classes."""

import numpy as np


class Hofstadter:
    r"""The Hofstadter model class.

    The Hamiltonian for the Hofstadter model is given as

    .. math::
        H = - \sum_{\braket{ij}} t e^{i \theta_{ij}} c_i^\dagger c_j + \mathrm{H.c.},

    where :math:`\braket{ij}` denotes nearest neighbors on a lattice in the xy-plane, :math:`t` are the hopping amplitudes, :math:`\theta_{ij}` are the Peierls phases, and :math:`c^{(\dagger)}` are the (creation)annihilation operators for bosons / spinless fermions.
    """

    def __init__(self, p, q, a0=1, t=1):
        """Constructor for the Hofstadter class.

        Parameters
        ----------
        p: int
            The numerator of the coprime flux density fraction.
        q: int
            The denominator of the coprime flux density fraction.
        a0: float
            The lattice constant (default=1).
        t: float
            The units of the hopping amplitudes (default=1).
        """

        self.p = p  #: int : The numerator of the coprime flux density fraction.
        self.q = q  #: int : The denominator of the coprime flux density fraction.
        self.a0 = a0  #: float :The lattice constant (default=1).
        self.t = t  #: float : The units of the hopping amplitudes (default=-1).

    def unit_cell(self):
        """The unit cell of the Hofstadter model.

        Returns
        -------
        num_bands_val: int
            The number of sites.
        avec_val: ndarray
            The lattice vectors.
        bvec_val: ndarray
            The reciprocal lattice vectors.
        sym_points_val: ndarray
            The high symmetry points.
        """

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
        """The Hamiltonian of the Hofstadter model.

        Parameters
        ----------
        k_val: ndarray
            The momentum vector.

        Returns
        -------
        Hamiltonian: ndarray
            The Hofstadter Hamiltonian matrix of dimension (num_bands, num_bands).
        """

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
            Hamiltonian[n][n] = -self.t * h(k_val, n)

        for n in range(self.q - 1):
            Hamiltonian[n][n + 1] = -self.t * np.exp(+1j * k_val[0] * self.a0)
            Hamiltonian[n + 1][n] = -self.t * np.exp(-1j * k_val[0] * self.a0)

        Hamiltonian[0][self.q - 1] = -self.t * np.exp(-1j * k_val[0] * self.a0)
        Hamiltonian[self.q - 1][0] = -self.t * np.exp(+1j * k_val[0] * self.a0)

        return Hamiltonian
