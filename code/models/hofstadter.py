"""The Hofstadter model classes."""

import numpy as np
import functions.models as fm


class Hofstadter:
    r"""The Hofstadter model class.

    The Hamiltonian for the Hofstadter model is given as

    .. math::
        H = - \sum_{\braket{ij}_n} t_n e^{i \theta_{ij}} c_i^\dagger c_j + \mathrm{H.c.},

    where :math:`\braket{ij}_n` denotes nth nearest neighbors on a lattice in the xy-plane, :math:`t_n` are the corresponding hopping amplitudes, :math:`\theta_{ij}` are the Peierls phases, and :math:`c^{(\dagger)}` are the particle (creation)annihilation operators.
    """

    def __init__(self, p, q, a0=1, t=None, lat="square"):
        """Constructor for the Hofstadter class.

        Parameters
        ----------
        p: int
            The numerator of the coprime flux density fraction.
        q: int
            The denominator of the coprime flux density fraction.
        a0: float
            The lattice constant (default=1).
        t: list
            The list of hopping amplitudes in order of ascending NN (default=[1]).
        lat: str
            The name of the lattice (default="square")
        """

        if t is None:
            t = [1]
        self.p = p  #: int : The numerator of the coprime flux density fraction.
        self.q = q  #: int : The denominator of the coprime flux density fraction.
        self.a0 = a0  #: float :The lattice constant (default=1).
        self.t = t  #: float : The units of the hopping amplitudes (default=[1]).
        self.lat = lat  #: str : The name of the lattice (default="square").

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
        if self.lat == "square":
            # lattice vectors
            a1 = self.a0 * np.array([1, 0])
            a2 = self.a0 * np.array([0, 1])
            avec_val = np.vstack((a1, a2))
            # reciprocal lattice vectors
            b1 = (2. * np.pi) / self.a0 * np.array([1, 0])
            b2 = (2. * np.pi) / self.a0 * np.array([0, 1])
            bvec_val = np.vstack((b1, b2))
            # basis vector
            basis = 1
            abasis1 = self.a0 * np.array([0, 0])
            abasisvec_val = np.array([abasis1])
            # reciprocal basis vectors
            bbasis1 = (2. * np.pi) / self.a0 * np.array([0, 0])
            bbasisvec_val = np.array([bbasis1])
            # lattice vectors (MUC)
            aMUC1 = num_bands_val * a1
            aMUC2 = a2
            aMUCvec_val = np.vstack((aMUC1, aMUC2))
            # reciprocal lattice vectors (MUC)
            bMUC1 = b1 / num_bands_val
            bMUC2 = b2
            bMUCvec_val = np.vstack((bMUC1, bMUC2))
            # cartesian vectors (for Peierls substitution)
            acart1 = self.a0 * np.array([1, 0])
            acart2 = self.a0 * np.array([0, 1])
            acartvec_val = np.vstack((acart1, acart2))
            # symmetry points
            GA = np.array([0, 0])
            Y = np.array([0, 0.5])
            S = np.array([0.5, 0.5])
            X = np.array([0.5, 0])
            sym_points_val = [GA, Y, S, X]
        elif self.lat == "triangular":
            # lattice vectors
            a1 = self.a0 * np.array([1, 0])
            a2 = self.a0 * np.array([1 / 2, np.sqrt(3) / 2])
            avec_val = np.vstack((a1, a2))
            # reciprocal lattice vectors
            b1 = (2. * np.pi) / self.a0 * np.array([1, -1 / np.sqrt(3)])
            b2 = (2. * np.pi) / self.a0 * np.array([0, 2 / np.sqrt(3)])
            bvec_val = np.vstack((b1, b2))
            # basis vector
            basis = 1
            abasis1 = self.a0 * np.array([0, 0])
            abasisvec_val = np.array([abasis1])
            # reciprocal basis vectors
            bbasis1 = (2. * np.pi) / self.a0 * np.array([0, 0])
            bbasisvec_val = np.array([bbasis1])
            # lattice vectors (MUC)
            aMUC1 = num_bands_val * a1
            aMUC2 = a2
            aMUCvec_val = np.vstack((aMUC1, aMUC2))
            # reciprocal lattice vectors (MUC)
            bMUC1 = b1 / num_bands_val
            bMUC2 = b2
            bMUCvec_val = np.vstack((bMUC1, bMUC2))
            # cartesian vectors (for Peierls substitution)
            acart1 = self.a0 * np.array([1/2, 0])
            acart2 = self.a0 * np.array([0, np.sqrt(3)/2])
            acartvec_val = np.vstack((acart1, acart2))
            # symmetry points
            GA = np.array([0, 0])
            Y = np.array([0, 0.5])
            S = np.array([0.5, 0.5])
            X = np.array([0.5, 0])
            sym_points_val = [GA, Y, S, X]
        elif self.lat == "honeycomb":
            # lattice vectors
            a1 = self.a0 * np.array([1, 0])
            a2 = self.a0 * np.array([1 / 2, np.sqrt(3) / 2])
            avec_val = np.vstack((a1, a2))
            # reciprocal lattice vectors
            b1 = (2. * np.pi) / self.a0 * np.array([1, -1 / np.sqrt(3)])
            b2 = (2. * np.pi) / self.a0 * np.array([0, 2 / np.sqrt(3)])
            bvec_val = np.vstack((b1, b2))
            # basis vector
            basis = 2
            abasis1 = self.a0 * np.array([0, 0])
            abasis2 = self.a0 * np.array([1/2, np.sqrt(3)/6])
            abasisvec_val = np.vstack((abasis1, abasis2))
            # reciprocal basis vectors
            bbasis1 = (2. * np.pi) / self.a0 * np.array([0, 0])
            bbasis2 = (2. * np.pi) / self.a0 * np.array([1, np.sqrt(3)])
            bbasisvec_val = np.vstack((bbasis1, bbasis2))
            # lattice vectors (MUC)
            aMUC1 = num_bands_val * a1
            aMUC2 = a2
            aMUCvec_val = np.vstack((aMUC1, aMUC2))
            # reciprocal lattice vectors (MUC)
            bMUC1 = b1 / num_bands_val
            bMUC2 = b2
            bMUCvec_val = np.vstack((bMUC1, bMUC2))
            # cartesian vectors (for Peierls substitution)
            acart1 = self.a0 * np.array([1/2, 0])
            acart2 = self.a0 * np.array([0, np.sqrt(3)/6])
            acartvec_val = np.vstack((acart1, acart2))
            # symmetry points
            GA = np.array([0, 0])
            Y = np.array([0, 0.5])
            S = np.array([0.5, 0.5])
            X = np.array([0.5, 0])
            sym_points_val = [GA, Y, S, X]

        self.avec = avec_val

        return (basis, num_bands_val, avec_val, bvec_val, abasisvec_val, bbasisvec_val,
                aMUCvec_val, bMUCvec_val, acartvec_val, sym_points_val)

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

        basis, _, avec, _, abasisvec, _, _, _, acartvec, _ = self.unit_cell()

        # nearest neighbors
        data0 = fm.nearest_neighbor_finder(avec, acartvec, abasisvec, self.t, 0, 0)
        data = [data0]
        sublattice = False
        for i, val in enumerate(data0):
            if val[10] == 1:  # check if any jumps change sublattice
                sublattice = True
                data.append(fm.nearest_neighbor_finder(avec, acartvec, abasisvec, self.t, val[4], val[5]))
        if not sublattice:
            basis = 1

        data = np.vstack(data)

        # sorted groups
        vec_group = fm.nearest_neighbor_sorter(data)

        # compute A_UC in mn units
        A_UC = np.linalg.norm(np.cross(avec[0], avec[1])) / np.linalg.norm(np.cross(acartvec[0], acartvec[1]))
        # compute cos(angle) between basis vectors
        cos_angle = np.vdot(avec[0], avec[1]) / (np.linalg.norm(avec[0])*np.linalg.norm(avec[1]))

        Hamiltonian = fm.Hamiltonian(self.t, self.p, self.q, A_UC, vec_group, k_val, cos_angle)

        return Hamiltonian, basis
