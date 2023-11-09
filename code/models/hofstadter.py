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

    def __init__(self, p, q, a0=1, t=None, lat="square", alpha=1, theta=(1, 2), period=1):
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
            The name of the lattice (default="square").
        alpha: float
            The anisotropy of the Bravais lattice vectors (default=1).
        theta: float
            The angle between Bravais lattice vectors in units of pi (default=0.5).
        period: int
            The factor by which to divide A_UC in the flux density (default=1).
        """

        if t is None:
            t = [1]
        self.p = p  #: int : The numerator of the coprime flux density fraction.
        self.q = q  #: int : The denominator of the coprime flux density fraction.
        self.a0 = a0  #: float :The lattice constant (default=1).
        self.t = t  #: float : The units of the hopping amplitudes (default=[1]).
        self.lat = lat  #: str : The name of the lattice (default="square").
        self.alpha = alpha  #: float : The anisotropy of the Bravais lattice vectors (default=1).
        self.theta = (theta[0]/theta[1])*np.pi  #: float : The angle between Bravais lattice vectors (default=pi/2).
        self.period = period  #: int : The factor by which to divide A_UC in the flux density (default=1).


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
        elif self.lat == "bravais":
            # lattice vectors
            a1 = self.a0 * np.array([1, 0])
            a2 = self.a0 * self.alpha * np.array([np.cos(self.theta), np.sin(self.theta)])
            avec_val = np.vstack((a1, a2))
            # reciprocal lattice vectors
            rec_factor = (2. * np.pi) / (a1[0]*a2[1] - a1[1]*a2[0])
            b1 = rec_factor * np.array([a2[1], -a2[0]])
            b2 = rec_factor * np.array([-a1[1], a1[0]])
            bvec_val = np.vstack((b1, b2))
            # basis vector
            basis = 1
            abasis1 = np.array([0, 0])
            abasisvec_val = np.array([abasis1])
            # reciprocal basis vectors
            bbasis1 = np.array([0, 0])
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
            # acart1 = np.array([np.max([a1[0], a2[0]]) - np.min([a1[0], a2[0]]), 0])
            # acart2 = np.array([0, np.max([a1[1], a2[1]]) - np.min([a1[1], a2[1]])])
            # acartvec_val = np.vstack((acart1, acart2))
            acart1 = np.array([None, 0])  # impossible to define gcd in general
            acart2 = np.array([0, a2[0]])
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
            a2 = self.a0 * self.alpha * np.array([np.cos(self.theta), np.sin(self.theta)])
            avec_val = np.vstack((a1, a2))
            # reciprocal lattice vectors
            rec_factor = (2. * np.pi) / (a1[0] * a2[1] - a1[1] * a2[0])
            b1 = rec_factor * np.array([a2[1], -a2[0]])
            b2 = rec_factor * np.array([-a1[1], a1[0]])
            bvec_val = np.vstack((b1, b2))
            # basis vector
            basis = 2
            abasis1 = np.array([0, 0])
            abasis2 = np.array([(1/2) * a1[0], (1/3) * a2[1]])  # centroid is 1/3 of the way up
            abasisvec_val = np.vstack((abasis1, abasis2))
            # reciprocal basis vectors
            bbasis1 = np.array([0, 0])
            bbasis2 = (2. * np.pi) / self.a0 * np.array([1/abasis2[0], 3/(2*abasis2[1])])
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
            # acart1 = self.a0 * np.array([1/2, 0])
            # acart2 = self.a0 * np.array([0, np.sqrt(3)/6])
            # acartvec_val = np.vstack((acart1, acart2))
            acart1 = np.array([(1/2)*a1[0], 0])
            acart2 = np.array([0, (1/3)*a2[1]])
            acartvec_val = np.vstack((acart1, acart2))
            # symmetry points
            K1 = np.array([2/3, 1/3])
            GA = np.array([0., 0.])
            MM = np.array([0.5, 0.5])
            K2 = np.array([1/3, 2/3])
            sym_points_val = [K1, GA, MM, K2]
        elif self.lat == "kagome":
            # lattice vectors
            a1 = self.a0 * np.array([1, 0])
            a2 = self.a0 * self.alpha * np.array([np.cos(self.theta), np.sin(self.theta)])
            avec_val = np.vstack((a1, a2))
            # reciprocal lattice vectors
            rec_factor = (2. * np.pi) / (a1[0] * a2[1] - a1[1] * a2[0])
            b1 = rec_factor * np.array([a2[1], -a2[0]])
            b2 = rec_factor * np.array([-a1[1], a1[0]])
            bvec_val = np.vstack((b1, b2))
            # basis vector
            basis = 3
            abasis1 = np.array([0, 0])
            abasis2 = np.array([(1/2) * a1[0], 0])
            abasis3 = np.array([(1/2) * a2[0], (1/2) * a2[1]])
            abasisvec_val = np.vstack((abasis1, abasis2, abasis3))
            # reciprocal basis vectors
            bbasis1 = np.array([0, 0])
            rec_factor = (2. * np.pi) / (abasis2[0] * abasis3[1] - abasis2[1] * abasis3[0])
            bbasis2 = rec_factor * np.array([abasis3[1], -abasis3[0]])
            bbasis3 = rec_factor * np.array([-abasis2[1], abasis2[0]])
            bbasisvec_val = np.vstack((bbasis1, bbasis2, bbasis3))
            # lattice vectors (MUC)
            aMUC1 = num_bands_val * a1
            aMUC2 = a2
            aMUCvec_val = np.vstack((aMUC1, aMUC2))
            # reciprocal lattice vectors (MUC)
            bMUC1 = b1 / num_bands_val
            bMUC2 = b2
            bMUCvec_val = np.vstack((bMUC1, bMUC2))
            # cartesian vectors (for Peierls substitution)
            # acart1 = self.a0 * np.array([1/2, 0])
            # acart2 = self.a0 * np.array([0, np.sqrt(3)/6])
            # acartvec_val = np.vstack((acart1, acart2))
            acart1 = np.array([(1/2)*a2[0], 0])
            acart2 = np.array([0, (1/2)*a2[1]])
            acartvec_val = np.vstack((acart1, acart2))
            # symmetry points
            K1 = np.array([2/3, 1/3])
            GA = np.array([0., 0.])
            MM = np.array([0.5, 0.5])
            K2 = np.array([1/3, 2/3])
            sym_points_val = [K1, GA, MM, K2]

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

        data0, bases = fm.nearest_neighbor_finder(avec, acartvec, abasisvec, self.t, 0, 0, 0)
        data = [data0]

        len_bases = len(bases)
        if len_bases > 1:
            for i, val in enumerate(bases[1:]):
                data_set, _ = fm.nearest_neighbor_finder(avec, acartvec, abasisvec, self.t, abasisvec[i+1][0], abasisvec[i+1][1], val)
                data.append(data_set)
        data = np.vstack(data)

        vec_group_matrix = np.zeros((len_bases, len_bases), dtype=object)
        for i in range(len_bases):  # initial sublattice
            mask_i = (data[:, 10] == i)
            data_mask_i = data[mask_i, :]
            for j in range(len_bases):  # final sublattice
                mask_j = (data_mask_i[:, 11] == j)
                data_mask_ij = data_mask_i[mask_j, :]
                vec_group_list = fm.nearest_neighbor_sorter(data_mask_ij)
                vec_group_matrix[i, j] = vec_group_list

        # compute A_UC in units of a
        A_UC = np.linalg.norm(avec[1]) / self.period
        # compute cos(angle) between basis vectors
        cos_angle = np.vdot(avec[0], avec[1]) / (np.linalg.norm(avec[0])*np.linalg.norm(avec[1]))

        Hamiltonian = fm.Hamiltonian(self.t, self.p, self.q, A_UC, cos_angle, vec_group_matrix, k_val)

        return Hamiltonian
