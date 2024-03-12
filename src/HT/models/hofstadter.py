"""The Hofstadter model classes."""

# --- external imports
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.patches import Polygon
from math import gcd
# --- internal imports
from HT.functions import models as fm


class Hofstadter:
    r"""The Hofstadter model class.

    The Hamiltonian for the Hofstadter model is given as

    .. math::
        H = - \sum_\kappa \sum_{\braket{ij}_\kappa} t_\kappa e^{i \theta_{ij}} c_i^\dagger c_j + \mathrm{H.c.},

    where :math:`\braket{\dots}_\kappa` denotes :math:`\kappa`-th nearest neighbors on some regular Euclidean lattice in the xy-plane, :math:`t_\kappa` are the corresponding hopping amplitudes, :math:`\theta_{ij}` are the Peierls phases, and :math:`c^{(\dagger)}` are the particle (creation)annihilation operators.
    """

    def __init__(self, p, q, a0=1, t=None, lat="bravais", alpha=1, theta=(1, 3), period=1):
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
            The name of the lattice (default="bravais").
        alpha: float
            The anisotropy of the Bravais lattice vectors (default=1).
        theta: tuple
            The angle between Bravais lattice vectors in units of pi (default=(1, 3)).
        period: int
            The factor by which to divide A_UC in the flux density (default=1).
        """

        if t is None:
            t = [1]
        self.p = p  #: int : The numerator of the coprime flux density fraction.
        self.q = q  #: int : The denominator of the coprime flux density fraction.
        self.a0 = a0  #: float :The lattice constant (default=1).
        self.t = t  #: list : The list of hopping amplitudes in order of ascending NN (default=[1]).
        self.lat = lat  #: str : The name of the lattice (default="bravais").
        self.alpha = alpha  #: float : The anisotropy of the Bravais lattice vectors (default=1).
        self.theta0 = theta[0]  #: int : The numerator of the fractional angle between Bravais lattice vectors in units of pi (default=1).
        self.theta1 = theta[1]  #: int : The denominator of the fractional angle between Bravais lattice vectors in units of pi (default=3).
        self.theta = (theta[0]/theta[1])*np.pi  #: float : The angle between Bravais lattice vectors (default=pi/3).
        self.period = period  #: int : The factor by which to divide A_UC in the flux density (default=1).

    def unit_cell(self):
        """The unit cell of the Hofstadter model.

        Returns
        -------
        num_bands_val: int
            The number of bands in the spectrum.
        avec_val: ndarray
            The lattice vectors.
        abasisvec_val: ndarray
            The basis vectors.
        bMUCvec_val: ndarray
            The reciprocal MUC vectors.
        sym_points_val: ndarray
            The high-symmetry/reference points.
        """

        if self.lat == "square":
            # lattice vectors
            a1 = self.a0 * np.array([1, 0])
            a2 = self.a0 * np.array([0, 1])
            avec_val = np.vstack((a1, a2))
            # basis vectors
            abasis1 = np.array([0, 0])
            abasisvec_val = np.array([abasis1])
            # lattice vectors (MUC)
            aMUC1 = a1
            aMUC2 = self.q * a2
            aMUCvec_val = np.vstack((aMUC1, aMUC2))
            # reciprocal lattice vectors (MUC)
            bMUCvec_val = fm.reciprocal_vectors(aMUCvec_val)
            # symmetry points
            GA = np.array([0, 0])
            X = np.array([0.5, 0])
            M = np.array([0.5, 0.5])
            Y = np.array([0, 0.5])
            sym_points_val = [("$\\Gamma$", GA), ("$X$", X), ("$M$", M), ("$Y$", Y)]
        elif self.lat == "triangular":
            # lattice vectors
            a1 = self.a0 * np.array([1, 0])
            a2 = self.a0 * np.array([1 / 2, np.sqrt(3) / 2])
            avec_val = np.vstack((a1, a2))
            # basis vectors
            abasis1 = np.array([0, 0])
            abasisvec_val = np.array([abasis1])
            # lattice vectors (MUC)
            aMUC1 = a1
            aMUC2 = self.q * a2
            aMUCvec_val = np.vstack((aMUC1, aMUC2))
            # reciprocal lattice vectors (MUC)
            bMUCvec_val = fm.reciprocal_vectors(aMUCvec_val)
            # symmetry points
            GA = np.array([0., 0.])
            K = np.array([2/3, 1/3])
            M = np.array([0.5, 0.5])
            Kp = np.array([1/3, 2/3])
            sym_points_val = [("$\\Gamma$", GA), ("$K$", K), ("$M$", M), ("$K'$", Kp)]
        elif self.lat == "bravais":
            # lattice vectors
            a1 = self.a0 * np.array([1, 0])
            a2 = self.a0 * self.alpha * np.array([np.cos(self.theta), np.sin(self.theta)])
            avec_val = np.vstack((a1, a2))
            # basis vectors
            abasis1 = np.array([0, 0])
            abasisvec_val = np.array([abasis1])
            # lattice vectors (MUC)
            aMUC1 = a1
            aMUC2 = self.q * a2
            aMUCvec_val = np.vstack((aMUC1, aMUC2))
            # reciprocal lattice vectors (MUC)
            bMUCvec_val = fm.reciprocal_vectors(aMUCvec_val)
            # symmetry points
            GA = np.array([0., 0.])
            M = np.array([0.5, 0.5])
            theta_gcd = gcd(self.theta0, self.theta1)
            theta = (self.theta0/theta_gcd, self.theta1/theta_gcd)
            if theta[1] == 3:  # multiples of pi/3 but not pi/2
                K = np.array([2/3, 1/3])
                Kp = np.array([1/3, 2/3])
                sym_points_val = [("$\\Gamma$", GA), ("$K$", K), ("$M$", M), ("$K'$", Kp)]
            else:
                X = np.array([0.5, 0])
                Y = np.array([0, 0.5])
                sym_points_val = [("$\\Gamma$", GA), ("$X$", X), ("$M$", M), ("$Y$", Y)]
        elif self.lat == "honeycomb":
            # lattice vectors
            a1 = self.a0 * np.array([1, 0])
            a2 = self.a0 * self.alpha * np.array([np.cos(self.theta), np.sin(self.theta)])
            avec_val = np.vstack((a1, a2))
            # basis vectors
            abasis1 = np.array([0, 0])
            abasis2 = np.array([(1 / 2) * a1[0], (1 / 3) * a2[1]])  # centroid is 1/3 of the way up
            abasisvec_val = np.vstack((abasis1, abasis2))
            # lattice vectors (MUC)
            aMUC1 = a1
            aMUC2 = self.q * a2
            aMUCvec_val = np.vstack((aMUC1, aMUC2))
            # reciprocal lattice vectors (MUC)
            bMUCvec_val = fm.reciprocal_vectors(aMUCvec_val)
            # symmetry points
            GA = np.array([0., 0.])
            K = np.array([2/3, 1/3])
            M = np.array([0.5, 0.5])
            Kp = np.array([1/3, 2/3])
            sym_points_val = [("$\\Gamma$", GA), ("$K$", K), ("$M$", M), ("$K'$", Kp)]
        elif self.lat == "kagome":
            # lattice vectors
            a1 = self.a0 * np.array([1, 0])
            a2 = self.a0 * self.alpha * np.array([np.cos(self.theta), np.sin(self.theta)])
            avec_val = np.vstack((a1, a2))
            # basis vectors
            abasis1 = np.array([0, 0])
            abasis2 = np.array([(1 / 2) * a1[0], 0])
            abasis3 = np.array([(1 / 2) * a2[0], (1 / 2) * a2[1]])
            abasisvec_val = np.vstack((abasis1, abasis2, abasis3))
            # lattice vectors (MUC)
            aMUC1 = a1
            aMUC2 = self.q * a2
            aMUCvec_val = np.vstack((aMUC1, aMUC2))
            # reciprocal lattice vectors (MUC)
            bMUCvec_val = fm.reciprocal_vectors(aMUCvec_val)
            # symmetry points
            GA = np.array([0., 0.])
            K = np.array([2/3, 1/3])
            M = np.array([0.5, 0.5])
            Kp = np.array([1/3, 2/3])
            Mp = np.array([0, 0.5])
            sym_points_val = [("$\\Gamma$", GA), ("$K$", K), ("$M$", M), ("$K'$", Kp), ("$M'$", Mp)]
        elif self.lat == "custom":
            # lattice vectors
            a1 = self.a0 * np.array([1, 0])
            a2 = self.a0 * self.alpha * np.array([np.cos(self.theta), np.sin(self.theta)])
            avec_val = np.vstack((a1, a2))
            # basis vectors
            abasis1 = np.array([0, 0])
            abasis2 = np.array([(1 / 2) * a1[0], 0])
            abasis3 = np.array([(1 / 2) * a2[0], (1 / 2) * a2[1]])
            abasisvec_val = np.vstack((abasis1, abasis2, abasis3))
            # lattice vectors (MUC)
            aMUC1 = a1
            aMUC2 = self.q * a2
            aMUCvec_val = np.vstack((aMUC1, aMUC2))
            # reciprocal lattice vectors (MUC)
            bMUCvec_val = fm.reciprocal_vectors(aMUCvec_val)
            # symmetry points
            GA = np.array([0., 0.])
            K = np.array([2/3, 1/3])
            M = np.array([0.5, 0.5])
            Kp = np.array([1/3, 2/3])
            Mp = np.array([0, 0.5])
            sym_points_val = [("$\\Gamma$", GA), ("$K$", K), ("$M$", M), ("$K'$", Kp), ("$M'$", Mp)]

        num_bands_val = len(abasisvec_val) * self.q

        return num_bands_val, avec_val, abasisvec_val, bMUCvec_val, sym_points_val

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

        if self.lat == "square" and len(self.t) == 1:
            Hamiltonian = fm.BasicSquareHamiltonian(self.t, self.p, self.q, k_val, self.period)
        elif self.lat == "triangular" and len(self.t) == 1:
            Hamiltonian = fm.BasicTriangularHamiltonian(self.t, self.p, self.q, k_val, self.period)
        elif self.lat == "honeycomb" and len(self.t) == 1 and self.alpha == 1 and self.theta0 == 1 and self.theta1 == 3:
            Hamiltonian = fm.BasicHoneycombHamiltonian(self.t, self.p, self.q, k_val, self.period)
        elif self.lat == "kagome" and len(self.t) == 1 and self.alpha == 1 and self.theta0 == 1 and self.theta1 == 3:
            Hamiltonian = fm.BasicKagomeHamiltonian(self.t, self.p, self.q, k_val, self.period)
        else:  # general case
            _, avec, abasisvec, _, _ = self.unit_cell()

            data0, bases = fm.nearest_neighbor_finder(avec, abasisvec, self.t, 0, 0, 0)
            data = [data0]

            len_bases = len(bases)
            if len_bases > 1:
                for i, val in enumerate(bases[1:]):
                    data_set, _ = fm.nearest_neighbor_finder(avec, abasisvec, self.t, abasisvec[i+1][0], abasisvec[i+1][1], val)
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

            # compute A_UC in units of a (scaled by periodicity factor)
            A_UC = np.linalg.norm(avec[1]) / self.period

            Hamiltonian = fm.Hamiltonian(self.t, self.p, self.q, A_UC, vec_group_matrix, k_val)

        return Hamiltonian

    def plot_lattice(self):
        """Plot the lattice."""

        # extract parameters
        _, avec, abasisvec, _, _ = self.unit_cell()
        t_list = self.t

        # create list of NN to consider from t_list
        numb_list = []
        for i, t in enumerate(t_list):
            if t != 0:
                numb_list.append(i + 1)

        # create grid of basis vectors from [-numb_t_max, numb_t_max]
        vectors = []
        vectors_unit = []
        for i in range(-5*numb_list[-1], 5*numb_list[-1] + 1):
            for j in range(-5*numb_list[-1], 5*numb_list[-1] + 1):
                r_unit = np.array([i, j])
                vectors_unit.append(r_unit)
                r = np.matmul(r_unit, avec)
                for idx, k in enumerate(abasisvec):
                    vectors.append([r + k, i, j, 0, idx])  # add basis index

        # construct figure
        fig_lat = plt.figure()
        fig_lat.canvas.manager.set_window_title('Lattice')
        ax = fig_lat.add_subplot(111)
        ax.set_xlabel("$x$")
        ax.set_ylabel("$y$")
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
        ax.set_title(f"$t$ = {t_list}")

        # plot grid points
        r_list = []
        for i, val in enumerate(vectors):
            r = np.linalg.norm(val[0])
            if r != 0:
                r_list.append(r)
            ax.scatter(val[0][0], val[0][1], c='k', s=5)

        # plot nearest neighbors
        r_list = np.sort(list(set(r_list)))
        r_used = []
        for i, r in enumerate(r_list):
            if i >= len(t_list):
                break
            elif t_list[i] != 0:
                circle = plt.Circle((0, 0), r, color=f"C{i}", fill=False)
                ax.add_patch(circle)
                r_used.append(r)

        # adjust x and y range
        UC_radius = np.linalg.norm(np.add(avec[0], avec[1]))
        max_rad = max([max(r_used), UC_radius])
        ax.set_xlim([-max_rad, max_rad])
        ax.set_ylim([-max_rad, max_rad])

        # plot lattice vectors
        for i, val in enumerate(avec):
            ax.annotate("", xy=(val[0], val[1]), xytext=(0, 0), arrowprops=dict(arrowstyle="->"))

        # plot unit cell
        UC = Polygon(((0, 0), (avec[0][0], avec[0][1]),
                      (avec[0][0] + avec[1][0], avec[0][1] + avec[1][1]), (avec[1][0], avec[1][1])),
                     fc=(1, 0, 0, 0.5), ec=(0, 0, 0, 0), lw=0)
        ax.add_artist(UC)
        ax.set_aspect('equal')

        return None
