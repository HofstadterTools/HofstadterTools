"""Functions for the model classes."""

# --- external imports
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from copy import deepcopy


def reciprocal_vectors(avec):
    r"""Finds the reciprocal lattice vectors in 2D.

    .. math::
        \begin{pmatrix}
        b_{1x} & b_{2x} \\
        b_{1y} & b_{2y}
        \end{pmatrix} =
        \frac{2\\\pi}{a_{1x}a_{2y} - a_{1y}a_{2x}}
        \begin{pmatrix}
        a_{2y} & -a_{1y} \\
        -a_{2x} & a_{1x}
        \end{pmatrix}

    Parameters
    ----------
    avec: ndarray
        The array of real-space lattice vectors.

    Returns
    -------
    bvec: ndarray
        The array of reciprocal-space lattice vectors.
    """

    if len(avec) != 2:
        raise ValueError("Reciprocal lattice vector function only works in 2D.")

    a1 = avec[0]
    a2 = avec[1]

    rec_factor = (2*np.pi) / (a1[0] * a2[1] - a1[1] * a2[0])
    b1 = rec_factor * np.array([a2[1], -a2[0]])
    b2 = rec_factor * np.array([-a1[1], a1[0]])
    bvec = np.vstack((b1, b2))

    return bvec


def nearest_neighbor_finder(avec, abasisvec, t_list, x_init, y_init, basis_init):
    """Finds the relevant nearest neighbors for a given lattice.

    Parameters
    ----------
    avec: ndarray
        The lattice vectors.
    abasisvec: ndarray
        The basis vectors.
    t_list: list
        The list of hopping amplitudes in order of ascending NN.
    x_init: float
        The initial x-coordinate.
    y_init: float
        The initial y-coordinate.
    basis_init: int
        The initial sublattice index.

    Returns
    -------
    data: ndarray
        The array of relevant nearest neighbors.
    bases: list
        The list of unique sublattice indices.
    """

    # --- Create list of NN to consider from t_list
    numb_list = []
    for i, t in enumerate(t_list):
        if t != 0:
            numb_list.append(i+1)

    # --- Create grid of basis vectors from [-numb_t_max, numb_t_max]
    abasisvec = [abasisvec[i] - abasisvec[basis_init] for i in range(len(abasisvec))]  # shift basis
    vectors = []
    vectors_unit = []
    for i in range(-numb_list[-1], numb_list[-1]+1):
        for j in range(-numb_list[-1], numb_list[-1]+1):
            r_unit = np.array([i, j])
            vectors_unit.append(r_unit)
            r = np.matmul(r_unit, avec)
            for idx, k in enumerate(abasisvec):
                vectors.append([r+k, i, j, basis_init, idx])  # add basis index

    # --- Define data array with info on each vector
    data = np.zeros((len(vectors), 13), dtype=object)
    for i, val in enumerate(vectors):
        data[i, 0] = round(np.linalg.norm(val[0]), 10)  # r (round so that we can use it for comparison)
        data[i, 1] = np.angle(val[0][0]+1j*val[0][1])  # phi
        data[i, 2] = x_init  # x_init
        data[i, 3] = y_init  # y_init
        data[i, 4] = y_init / avec[1][1]  # y_init / a2_y
        data[i, 5] = val[0][0]  # dx
        data[i, 6] = val[0][1]  # dy
        data[i, 7] = val[0][1] / avec[1][1]  # dy / a2_y
        data[i, 8] = val[1]  # dI
        data[i, 9] = val[2]  # dJ
        data[i, 10] = val[3]  # sub_init
        data[i, 11] = val[4]  # sub_final

    # --- Extract the NN groups (filter data based on radius)
    data = data[data[:, 0].argsort()]  # sort by increasing r
    # delete the first r=0 row
    mask = (data[:, 0] != 0)
    data = data[mask, :]
    radii = np.sort(list(set(data[:, 0])))
    # label the NN group
    for i, row in enumerate(data):
        for j, r in enumerate(radii):
            if row[0] == r:
                data[i, 12] = j+1  # NN group
    select_radii = [radii[i - 1] for i in numb_list]
    # delete rows with other radii
    rows_to_delete = []
    for i, row in enumerate(data):
        if row[0] not in select_radii:
            rows_to_delete.append(i)
    data = np.delete(data, rows_to_delete, axis=0)

    # --- Extract bases set
    bases = []
    for i, val in enumerate(data):
        bases.append(val[10])
        bases.append(val[11])
    bases = np.sort(list(set(bases)))

    return data, bases


def nearest_neighbor_sorter(data_array):
    """Sorts the relevant nearest neighbors array.

    Parameters
    ----------
    data_array: ndarray
        The array of relevant nearest neighbors.

    Returns
    -------
    grouped_paths: ndarray
        The array of relevant nearest neighbors grouped by dJ.
    """

    # count the number of dJ
    dJ_list = []
    for i, val in enumerate(data_array):
        dJ_list.append(val[9])
    dJ_list = np.sort(list(set(dJ_list)))
    numb_dJ = len(dJ_list)

    # group paths by dJ
    grouped_paths = np.zeros(numb_dJ, dtype=object)
    for i in range(numb_dJ):
        grouped_paths[i] = []

    for i, dJval in enumerate(dJ_list):
        for j, val in enumerate(data_array):
            if val[9] == dJval:
                grouped_paths[i].append(val)
        grouped_paths[i].append(dJval)

    return grouped_paths


def peierls_factor(nphi, dx, y_cart, dy_cart, A_UC):
    r"""The Peierls factor.

    The Peierls factor in Landau gauge :math:`\\\mathbf{A}=-By\hat{\\\mathbf{e}}_x` is given by

    .. math::
        e^{\\\mathrm{i}\theta_{ij}} = \exp\left[ -\frac{2\pi\\\mathrm{i}n_\phi}{A} \Delta X \left( Y_i + \frac{\Delta Y}{2} \right) \right],

    where :math:`\theta_{ij}` is the Peierls phase from site :math:`i=(X_i, Y_i)`. to :math:`j=(X_j, Y_j)`, :math:`\Delta X = X_j - X_i`, :math:`\Delta Y = Y_j - Y_i`, :math:`n_\phi` is the flux density, and :math:`A` is the area factor to make the expression dimensionless. :cite:`Peierls33`

    Parameters
    ----------
    nphi: float
        The flux density.
    dx: float
        The change in x-coordinates.
    y_cart: int
        The initial y-coordinate in units of a2[1].
    dy_cart: int
        The change in y-coordinates in units of a2[1].
    A_UC: float
        The unit cell area in units of a2 (possibly scaled by a periodicity factor).

    Returns
    -------
    factor: complex
        The Peierls factor.
    """

    phase = - 2 * np.pi * nphi * dx * (y_cart + dy_cart/2) / A_UC
    factor = np.exp(1j * phase)

    return factor


def diag_func(t_val, p_val, q_val, A_UC_val, vec_group, k_val, dJ_val, J_idx_val):
    r"""The diagonal function.

    The function that populates the diagonals of the Harper matrix is given by

    .. math::
        \Lambda_{l, n} = - \sum_\kappa \sum_{\langle ij \rangle_{\kappa}^l} t_\kappa e^{\\\mathrm{i}\theta_{ij}} e^{\\\mathrm{i}\\\mathbf{k}\cdot\\\mathbf{r}},

    where :math:`\langle \dots \rangle^l_\kappa` denotes the subset of :math:`\kappa`-th nearest neighbors with a net :math:`y` unit cell displacement of :math:`l`, :math:`\theta_{ij}` is the Peierls phase, :math:`\\\mathbf{k}` is the momentum vector, and :math:`\\\mathbf{r}` is the displacement vector.

    Parameters
    ----------
    t_val: list
        The list of hopping amplitudes in order of ascending NN.
    p_val: int
        The numerator of the flux density.
    q_val: int
        The denominator of the flux density.
    A_UC_val: float
        The unit cell area in units of a2 (possibly scaled by a periodicity factor).
    vec_group: ndarray
        The array of relevant nearest neighbors grouped by dJ.
    k_val: ndarray
        The momentum vector.
    dJ_val: int
        The y-displacement in terms of unit cells.
    J_idx_val: int
        The y-position in terms of unit cells.

    Returns
    -------
    term: complex
        The diagonal term.
    """

    nphi = p_val/q_val
    term = 0
    for idx, val in enumerate(vec_group):
        if val[-1] == dJ_val:  # extract rows with appropriate dJ
            for k, val2 in enumerate(val[:-1]):  # for each vector in path group
                NN_group = int(val2[12])
                term += - t_val[NN_group - 1] * (peierls_factor(nphi, val2[5], J_idx_val + val2[4], val2[7], A_UC_val)
                                                 * np.exp(1j * np.vdot(np.array([val2[5], val2[6]]), k_val)))

    return term


def Hamiltonian(t, p, q, A_UC, vec_group_matrix, k):
    r"""The generalized Hofstadter Hamiltonian.

    The generalized Hofstadter Hamiltonian is given by the :math:`N_b\times N_b` block matrix

    .. math::
        H = \begin{pmatrix}
        H^{00} & H^{01} & \dots \\
        H^{10} & H^{11} & \dots \\
        \vdots & \vdots & \ddots
        \end{pmatrix},

    where :math:`N_b` is the number of sites in the basis. Each submatrix :math:`H^{\alpha\beta}` has dimensions :math:`q\times q` and represents the Harper equation for hoppings from the :math:`\alpha` to the :math:`\beta` sublatttice, such that

    .. math::
        H = \begin{pmatrix}
        \Lambda_{0,0} & \Lambda_{0,1} & \dots \\
        \Lambda_{1,0} & \Lambda_{1,1} & \dots \\
        \vdots & \vdots & \ddots
        \end{pmatrix} +
        \begin{pmatrix}
        \ddots & \Lambda_{0, q-1}^* & \Lambda_{0, q}^* \\
        \Lambda_{q-1, 0} & \ddots & \Lambda_{1, q}^* \\
        \Lambda_{q, 0} & \Lambda_{q, 1} & \ddots
        \end{pmatrix},

    where :math:`\Lambda_{l, n}` is the diagonal function, and we have dropped the sublattice superscripts for readability. Note that we only populate the first matrix with unique hoppings in the positive J direction. If there are inter-unit-cell hoppings that are related by Hermitian conjugation in the same Harper equation, then the lower triangular matrix will simply be the conjugate of the upper triangular matrix. The second matrix represents rolled over boundary terms.

    Parameters
    ----------
    t: list
        The list of hopping amplitudes in order of ascending NN.
    p: int
        The numerator of the flux density.
    q: int
        The denominator of the flux density.
    A_UC: float
        The unit cell area in units of a2 (possibly scaled by a periodicity factor).
    vec_group_matrix: ndarray
        The matrix of grouped nearest neighbor arrays.
    k: ndarray
        The momentum vector.

    Returns
    -------
    Hamiltonian: ndarray
        The Hofstadter Hamiltonian matrix of dimension :math:`N_b q \times N_b q`.
    """

    I = np.shape(vec_group_matrix)[0]
    J = np.shape(vec_group_matrix)[1]

    Ham_matrix = []
    for i in range(I):
        Ham_row = []
        for j in range(J):
            Hamiltonian = np.zeros((q, q), dtype=np.complex128)

            dJ_list = []
            for term in vec_group_matrix[i, j]:
                dJ_list.append(term[-1])

            HC_flag = False
            for k1, val in enumerate(dJ_list):  # remove negative unit cell hoppings (for H.c. cases)
                for k2, val2 in enumerate(dJ_list):
                    if val == -val2 and val != 0:
                        HC_flag = True
                        if val2 < 0:
                            del dJ_list[k2]
                        else:
                            del dJ_list[k1]

            for dJ in dJ_list:
                # upper_diag_array
                diag_array = np.array([diag_func(t, p, q, A_UC, vec_group_matrix[i, j], k, dJ, J_idx) for J_idx in range(q)])
                Hamiltonian += np.roll(np.diag(diag_array), abs(dJ), axis=int((np.sign(dJ)+1)/2))
                # lower_diag_array
                if HC_flag and dJ > 0:
                    Hamiltonian += np.roll(np.diag(np.conj(diag_array)), abs(dJ), axis=0)

            Ham_row.append(Hamiltonian)
        Ham_matrix.append(Ham_row)
    Hamiltonian = np.block(Ham_matrix)

    return Hamiltonian


def BasicSquareHamiltonian(t, p, q, k, period):
    r"""The basic square lattice Hofstadter Hamiltonian.

    The Hofstadter Hamiltonian for the square lattice with nearest-neighbor hopping.

    .. note::

        We have hardcoded this Hamiltonian because it is commonly used and this function is faster than calling the generic Hamiltonian method.

    Parameters
    ----------
    t: list
        The list of hopping amplitudes in order of ascending NN.
    p: int
        The numerator of the flux density.
    q: int
        The denominator of the flux density.
    k: ndarray
        The momentum vector.
    period: int
        The factor by which to divide A_UC in the flux density.

    Returns
    -------
    Hamiltonian: ndarray
        The Hofstadter Hamiltonian matrix of dimension :math:`q \times q`.
    """

    Hamiltonian = np.zeros((q, q), dtype=np.complex128)
    nphi = p / q

    def A(t_val, nphi_val, m_val, k_val):
        value = -t_val*np.exp(-1j*2*np.pi*period*nphi_val*m_val + 1j*k_val[0])
        return value

    def B_plus(t_val, k_val):
        value = -t_val*np.exp(1j*k_val[1])
        return value

    diag_array = np.array([A(t[0], nphi, m, k) for m in range(q)])
    Hamiltonian += np.roll(np.diag(diag_array), 0, axis=1)

    upper_diag_array = np.array([B_plus(t[0], k) for _ in range(q)])
    Hamiltonian += np.roll(np.diag(upper_diag_array), 1, axis=1)

    Hamiltonian = Hamiltonian + np.conj(Hamiltonian).T

    return Hamiltonian


def BasicTriangularHamiltonian(t, p, q, k, period):
    r"""The basic triangular lattice Hofstadter Hamiltonian.

    The Hofstadter Hamiltonian for the triangular lattice with nearest-neighbor hopping.

    .. note::

        We have hardcoded this Hamiltonian because it is commonly used and this function is faster than calling the generic Hamiltonian method.

    Parameters
    ----------
    t: list
        The list of hopping amplitudes in order of ascending NN.
    p: int
        The numerator of the flux density.
    q: int
        The denominator of the flux density.
    k: ndarray
        The momentum vector.
    period: int
        The factor by which to divide A_UC in the flux density.

    Returns
    -------
    Hamiltonian: ndarray
        The Hofstadter Hamiltonian matrix of dimension :math:`q \times q`.
    """

    Hamiltonian = np.zeros((q, q), dtype=np.complex128)
    nphi = p / q

    def A(t_val, nphi_val, m_val, k_val):
        value = -t_val*np.exp(-1j*2*np.pi*period*nphi_val*m_val + 1j*k_val[0])
        return value

    def B_plus(t_val, nphi_val, m_val, k_val):
        value = (-t_val*np.exp(-1j*np.pi*period*nphi_val*(m_val + 1/2) + 1j*(0.5*k_val[0] + np.sqrt(3)*k_val[1]/2))
                 -t_val*np.exp(+1j*np.pi*period*nphi_val*(m_val + 1/2) + 1j*(-0.5*k_val[0] + np.sqrt(3)*k_val[1]/2)))
        return value

    diag_array = np.array([A(t[0], nphi, m, k) for m in range(q)])
    Hamiltonian += np.roll(np.diag(diag_array), 0, axis=1)

    upper_diag_array = np.array([B_plus(t[0], nphi, m, k) for m in range(q)])
    Hamiltonian += np.roll(np.diag(upper_diag_array), 1, axis=1)

    Hamiltonian = Hamiltonian + np.conj(Hamiltonian).T

    return Hamiltonian


def BasicHoneycombHamiltonian(t, p, q, k, period):
    r"""The basic honeycomb lattice Hofstadter Hamiltonian.

    The Hofstadter Hamiltonian for the honeycomb lattice with nearest-neighbor hopping.

    .. note::

        We have hardcoded this Hamiltonian because it is commonly used and this function is faster than calling the generic Hamiltonian method.

    Parameters
    ----------
    t: list
        The list of hopping amplitudes in order of ascending NN.
    p: int
        The numerator of the flux density.
    q: int
        The denominator of the flux density.
    k: ndarray
        The momentum vector.
    period: int
        The factor by which to divide A_UC in the flux density.

    Returns
    -------
    Hamiltonian: ndarray
        The Hofstadter Hamiltonian matrix of dimension :math:`2q \times 2q`.
    """

    def AB_block_func(t_val, p_val, q_val, k_val):
        ham = np.zeros((q_val, q_val), dtype=np.complex128)
        nphi = p_val / q_val

        def A(t_val, nphi_val, m_val, k_val):
            value = (-t_val * np.exp(-1j*np.pi*period*nphi_val*(m_val + 1/6) + 1j*(+k_val[0]/2 + np.sqrt(3)*k_val[1]/6))
                     -t_val * np.exp(+1j*np.pi*period*nphi_val*(m_val + 1/6) + 1j*(-k_val[0]/2 + np.sqrt(3)*k_val[1]/6)))
            return value

        def B_minus(t_val, k_val):
            value = -t_val * np.exp(-1j*2*np.sqrt(3)*k_val[1]/6)
            return value

        upper_diag_array = np.array([A(t_val, nphi, m, k_val) for m in range(q_val)])
        ham += np.roll(np.diag(upper_diag_array), 0, axis=1)

        upper_diag_array2 = np.array([B_minus(t_val, k_val) for _ in range(q_val)])
        ham += np.roll(np.diag(upper_diag_array2), 1, axis=0)

        return ham

    def BA_block_func(t_val, p_val, q_val, k_val):
        ham = np.zeros((q_val, q_val), dtype=np.complex128)
        nphi = p_val / q_val

        def A(t_val, nphi_val, m_val, k_val):
            value = (-t_val * np.exp(-1j*np.pi*period*nphi_val*(m_val + 1/6) + 1j*(+k_val[0]/2 - np.sqrt(3)*k_val[1]/6))
                     -t_val * np.exp(+1j*np.pi*period*nphi_val*(m_val + 1/6) + 1j*(-k_val[0]/2 - np.sqrt(3)*k_val[1]/6)))
            return value

        def B_plus(t_val, k_val):
            value = -t_val * np.exp(+1j*2*np.sqrt(3)*k_val[1]/6)
            return value

        upper_diag_array = np.array([A(t_val, nphi, m, k_val) for m in range(q_val)])
        ham += np.roll(np.diag(upper_diag_array), 0, axis=1)

        upper_diag_array2 = np.array([B_plus(t_val, k_val) for _ in range(q_val)])
        ham += np.roll(np.diag(upper_diag_array2), 1, axis=1)

        return ham

    AA_block = np.zeros((q, q))
    AB_block = AB_block_func(t[0], p, q, k)
    BA_block = BA_block_func(t[0], p, q, k)
    BB_block = np.zeros((q, q))

    upper = np.concatenate((AA_block, AB_block), axis=1)
    lower = np.concatenate((BA_block, BB_block), axis=1)
    Hamiltonian = np.concatenate((upper, lower), axis=0)

    return Hamiltonian


def BasicKagomeHamiltonian(t, p, q, k, period):
    r"""The basic kagome lattice Hofstadter Hamiltonian.

    The Hofstadter Hamiltonian for the kagome lattice with nearest-neighbor hopping.

    .. note::

        We have hardcoded this Hamiltonian because it is commonly used and this function is faster than calling the generic Hamiltonian method.

    Parameters
    ----------
    t: list
        The list of hopping amplitudes in order of ascending NN.
    p: int
        The numerator of the flux density.
    q: int
        The denominator of the flux density.
    k: ndarray
        The momentum vector.
    period: int
        The factor by which to divide A_UC in the flux density.

    Returns
    -------
    Hamiltonian: ndarray
        The Hofstadter Hamiltonian matrix of dimension :math:`3q \times 3q`.
    """

    def AB_block_func(t_val, p_val, q_val, k_val):
        ham = np.zeros((q_val, q_val), dtype=np.complex128)
        nphi = p_val / q_val

        def A(t_val, nphi_val, m_val, k_val):
            value = (-t_val*np.exp(-1j*np.pi*period*nphi_val*m_val + 1j*k_val[0]/2)
                     -t_val*np.exp(+1j*np.pi*period*nphi_val*m_val - 1j*k_val[0]/2))
            return value

        upper_diag_array = np.array([A(t_val, nphi, m, k_val) for m in range(q_val)])
        ham += np.roll(np.diag(upper_diag_array), 0, axis=1)

        return ham

    def AC_block_func(t_val, p_val, q_val, k_val):
        ham = np.zeros((q_val, q_val), dtype=np.complex128)
        nphi = p_val / q_val

        def A(t_val, nphi_val, m_val, k_val):
            value = -t_val*np.exp(-1j*np.pi*period*nphi_val*0.5*(m_val + 1/4) + 1j*(k_val[0]/4 + np.sqrt(3)*k_val[1]/4))
            return value

        def B_minus(t_val, nphi_val, m_val, k_val):
            value = -t_val*np.exp(+1j*np.pi*period*nphi_val*0.5*(m_val - 1/4) - 1j*(k_val[0]/4 + np.sqrt(3)*k_val[1]/4))
            return value

        upper_diag_array = np.array([A(t_val, nphi, m, k_val) for m in range(q_val)])
        ham += np.roll(np.diag(upper_diag_array), 0, axis=1)

        upper_diag_array2 = np.array([B_minus(t_val, nphi, m, k_val) for m in range(q_val)])
        ham += np.roll(np.diag(upper_diag_array2), 1, axis=0)

        return ham

    def BA_block_func(t_val, p_val, q_val, k_val):
        ham = np.zeros((q_val, q_val), dtype=np.complex128)
        nphi = p_val / q_val

        def A(t_val, nphi_val, m_val, k_val):
            value = (-t_val*np.exp(-1j*np.pi*period*nphi_val*m_val + 1j*k_val[0]/2)
                     -t_val*np.exp(+1j*np.pi*period*nphi_val*m_val - 1j*k_val[0]/2))
            return value

        upper_diag_array = np.array([A(t_val, nphi, m, k_val) for m in range(q_val)])
        ham += np.roll(np.diag(upper_diag_array), 0, axis=1)

        return ham

    def BC_block_func(t_val, p_val, q_val, k_val):
        ham = np.zeros((q_val, q_val), dtype=np.complex128)
        nphi = p_val / q_val

        def A(t_val, nphi_val, m_val, k_val):
            value = -t_val*np.exp(+1j*np.pi*period*nphi_val*0.5*(m_val + 1/4) + 1j*(-k_val[0]/4+np.sqrt(3)*k_val[1]/4))
            return value

        def B_minus(t_val, nphi_val, m_val, k_val):
            value = -t_val*np.exp(-1j*np.pi*period*nphi_val*0.5*(m_val - 1/4) + 1j*(+k_val[0]/4-np.sqrt(3)*k_val[1]/4))
            return value

        upper_diag_array = np.array([A(t_val, nphi, m, k_val) for m in range(q_val)])
        ham += np.roll(np.diag(upper_diag_array), 0, axis=1)

        upper_diag_array2 = np.array([B_minus(t_val, nphi, m, k_val) for m in range(q_val)])
        ham += np.roll(np.diag(upper_diag_array2), 1, axis=0)

        return ham

    def CA_block_func(t_val, p_val, q_val, k_val):
        ham = np.zeros((q_val, q_val), dtype=np.complex128)
        nphi = p_val / q_val

        def A(t_val, nphi_val, m_val, k_val):
            value = -t_val*np.exp(+1j*np.pi*period*nphi_val*0.5*(m_val + 1/4) - 1j*(k_val[0]/4+np.sqrt(3)*k_val[1]/4))
            return value

        def B_plus(t_val, nphi_val, m_val, k_val):
            value = -t_val*np.exp(-1j*np.pi*period*nphi_val*0.5*(m_val + 3/4) + 1j*(k_val[0]/4+np.sqrt(3)*k_val[1]/4))
            return value

        upper_diag_array = np.array([A(t_val, nphi, m, k_val) for m in range(q_val)])
        ham += np.roll(np.diag(upper_diag_array), 0, axis=1)

        upper_diag_array2 = np.array([B_plus(t_val, nphi, m, k_val) for m in range(q_val)])
        ham += np.roll(np.diag(upper_diag_array2), 1, axis=1)

        return ham

    def CB_block_func(t_val, p_val, q_val, k_val):
        ham = np.zeros((q_val, q_val), dtype=np.complex128)
        nphi = p_val / q_val

        def A(t_val, nphi_val, m_val, k_val):
            value = -t_val*np.exp(-1j*np.pi*period*nphi_val*0.5*(m_val + 1/4) + 1j*(k_val[0]/4-np.sqrt(3)*k_val[1]/4))
            return value

        def B_plus(t_val, nphi_val, m_val, k_val):
            value = -t_val*np.exp(+1j*np.pi*period*nphi_val*0.5*(m_val + 3/4) + 1j*(-k_val[0]/4+np.sqrt(3)*k_val[1]/4))
            return value

        upper_diag_array = np.array([A(t_val, nphi, m, k_val) for m in range(q_val)])
        ham += np.roll(np.diag(upper_diag_array), 0, axis=1)

        upper_diag_array2 = np.array([B_plus(t_val, nphi, m, k_val) for m in range(q_val)])
        ham += np.roll(np.diag(upper_diag_array2), 1, axis=1)

        return ham

    AA_block = np.zeros((q, q))
    AB_block = AB_block_func(t[0], p, q, k)
    AC_block = AC_block_func(t[0], p, q, k)
    #
    BA_block = BA_block_func(t[0], p, q, k)
    BB_block = np.zeros((q, q))
    BC_block = BC_block_func(t[0], p, q, k)
    #
    CA_block = CA_block_func(t[0], p, q, k)
    CB_block = CB_block_func(t[0], p, q, k)
    CC_block = np.zeros((q, q))

    upper = np.concatenate((AA_block, AB_block, AC_block), axis=1)
    middle = np.concatenate((BA_block, BB_block, BC_block), axis=1)
    lower = np.concatenate((CA_block, CB_block, CC_block), axis=1)
    Hamiltonian = np.concatenate((upper, middle, lower), axis=0)

    return Hamiltonian


if __name__ == '__main__':

    t = [1, 0, -0.25]

    vec_group = nearest_neighbor_finder("square", t)

    print(vec_group)

    ###

    num_bands = 5
    num_samples = 101

    b1 = (2. * np.pi) * np.array([1 / num_bands, 0])
    b2 = (2. * np.pi) * np.array([0, 1])
    bvec = np.vstack((b1, b2))

    eigenvalues = np.zeros((num_bands, num_samples, num_samples))  # real
    eigenvectors = np.zeros((num_bands, num_bands, num_samples, num_samples), dtype=np.complex128)  # complex
    for band in range(num_bands):
        for idx_x in range(num_samples):
            frac_kx = idx_x / (num_samples - 1)
            for idx_y in range(num_samples):
                frac_ky = idx_y / (num_samples - 1)
                k = np.matmul(np.array([frac_kx, frac_ky]), bvec)
                # print("k = ", k)
                eigvals, eigvecs = np.linalg.eigh(Hamiltonian(t, 1, num_bands, vec_group, k))
                idx = np.argsort(eigvals)
                eigenvalues[band, idx_x, idx_y] = eigvals[idx[band]]
                eigenvectors[:, band, idx_x, idx_y] = eigvecs[:, idx[band]]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    idx_x = np.linspace(0, num_samples - 1, num_samples - 1, dtype=int)
    idx_y = np.linspace(0, num_samples - 1, num_samples - 1, dtype=int)
    kx, ky = np.meshgrid(idx_x, idx_y)
    for i in range(num_bands):
        ax.plot_surface(kx, ky, eigenvalues[i, kx, ky], alpha=0.5)
    ax.set_xlabel('$k_1/|\\mathbf{b}_1|$')
    ax.set_ylabel('$k_2/|\\mathbf{b}_2|$')
    ax.set_zlabel('$E$')

    def normalize(value, tick_number):
        if value == 0:
            return "$0$"
        elif value == num_samples - 1:
            return "$1$"
        else:
            return f"${value / (num_samples - 1):.1g}$"

    ax.xaxis.set_major_formatter(plt.FuncFormatter(normalize))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(normalize))

    plt.show()
