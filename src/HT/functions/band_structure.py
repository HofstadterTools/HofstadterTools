"""Functions for band structure calculations."""

# --- external imports
import numpy as np
import cmath as cm


def _principal(z):
    r"""
    Compute the principal branch of a complex number :math:`z`, such that :math:`-\\\pi<\text{Im}(z)\leq\\\pi`.

    Parameters
    ----------
    z: complex
        The original complex number.

    Returns
    -------
    z: complex
        The principal value of the complex number.
    """

    re = np.real(z)
    im = np.imag(z)
    if im <= -np.pi:
        z = re + 1j * (np.pi - np.abs((im + np.pi)) % (2 * np.pi))
    elif im > np.pi:
        z = re + 1j * (-np.pi + np.abs((im + np.pi)) % (2 * np.pi))

    return z


def U(var_num, _eigenvectors, _band, _idx_x, _idx_y, _group_size):
    r"""Compute the link variable.

    The normalized link variables are defined as

    .. math::
        \tilde{\mathcal{U}}_\gamma(\mathbf{k}_\alpha) = \frac{\det U_\gamma(\mathbf{k}_\alpha)}{|\det U_\gamma(\mathbf{k}_\alpha)|}, \;\;\; \gamma = \{1, 2\},

    with link matrices

    .. math::
        U_\gamma(\mathbf{k}_\alpha) =
        \begin{pmatrix}
        \braket{u_1(\mathbf{k}_\alpha)|u_1(\mathbf{k}_\alpha+\hat{\mathbf{e}_\gamma})} & \dots & \braket{u_1(\mathbf{k}_\alpha)|u_M(\mathbf{k}_\alpha+\hat{\mathbf{e}_\gamma})} \\
        \vdots & \ddots & \vdots \\
         \braket{u_M(\mathbf{k}_\alpha)|u_1(\mathbf{k}_\alpha+\hat{\mathbf{e}_\gamma})} & \dots & \braket{u_M(\mathbf{k}_\alpha)|u_M(\mathbf{k}_\alpha+\hat{\mathbf{e}_\gamma})}
        \end{pmatrix}.

    Here, :math:`\mathbf{k}_\alpha` is the discretized momentum vector, :math:`\{\hat{\mathbf{e}}_1, \hat{\mathbf{e}}_2\}` are linearly independent unit vectors in the momentum grid, and :math:`\ket{u(\mathbf{k}_\alpha)}` is the eigenvector at momentum :math:`\mathbf{k}_\alpha`. The link variables are constructed for :math:`M` touching bands. :cite:`Fukui05`

    .. note::
        Input eigenvectors are already normalized from :class:`numpy.linalg.eig`.

    Parameters
    ----------
    var_num: [1, 2]
        The link variable number.
    _eigenvectors: ndarray
        The array of eigenvectors with dimension (num_bands, num_bands, num_samples, num_samples).
    _band: int
        The band number. If part of a band group, this must refer to the lowest band of the group.
    _idx_x: int
        The x-momentum, with respect to the discretized grid.
    _idx_y: int
        The y-momentum, with respect to the discretized grid.
    _group_size: int
        The number of bands in the band group.

    Returns
    -------
    link_var: complex
        The U(1) link variable.
    """

    _num_samples = np.shape(_eigenvectors)[2]

    link_matrix = np.zeros((_group_size, _group_size), dtype=complex)
    for i in range(_group_size):
        for j in range(_group_size):
            vec1 = _eigenvectors[:, _band + i, _idx_x, _idx_y]
            if var_num == 1:
                vec2 = _eigenvectors[:, _band + j, (_idx_x + 1) % _num_samples, _idx_y]
            elif var_num == 2:
                vec2 = _eigenvectors[:, _band + j, _idx_x, (_idx_y + 1) % _num_samples]
            else:
                raise ValueError("Link variable number must be in [1, 2].")
            link_matrix[i, j] = np.conj(vec1).dot(vec2)
    link_var = np.linalg.det(link_matrix)
    return link_var


def berry_curv(_eigenvectors, _band, _idx_x, _idx_y, _group_size=1):
    r"""
    Compute the Berry curvature.

    The Berry curvature around a plaquette is computed using the formula from :cite:`Fukui05` (example applications in :cite:`SoluyanovPhD, AidelsburgerPhD`), such that

    .. math::
       \mathcal{B}_{12}(\mathbf{k}_\alpha) \equiv - \text{Im}\;\log\;(\tilde{\mathcal{U}}_1(\mathbf{k}_\alpha)\tilde{\mathcal{U}}_2(\mathbf{k}_\alpha+\hat{\mathbf{e}}_1)\tilde{\mathcal{U}}_1(\mathbf{k}_\alpha+\hat{\mathbf{e}}_2)^{-1}\tilde{\mathcal{U}}_2(\mathbf{k}_\alpha)^{-1}),

    where :math:`\tilde{\mathcal{U}}` are the normalized link variables. The Berry curvature at a point :math:`\mathbf{k}` can then be computed by taking the limit of small plaquette size.

    .. note::
        The Berry curvature is defined within the principal branch of the logarithm. For example, the corresponding log sum formula for a single band would be

        .. math::
            \mathcal{B}_{12}(\mathbf{k}_\alpha) = - \text{Im}\;\mathcal{P}\;(\log\tilde{\mathcal{U}}_1(\mathbf{k}_\alpha) +\log\tilde{\mathcal{U}}_2(\mathbf{k}_\alpha+\hat{\mathbf{e}}_1)-\log\tilde{\mathcal{U}}_1(\mathbf{k}_\alpha+\hat{\mathbf{e}}_2)-\log\tilde{\mathcal{U}}_2(\mathbf{k}_\alpha)),

        where :math:`\mathcal{P}` denotes the principal value of the complex number :math:`z`, such that :math:`-\\\pi<\text{Im}(z)\leq\\\pi`.

    Parameters
    ----------
    _eigenvectors: ndarray
        The array of eigenvectors with dimension (num_bands, num_bands, num_samples, num_samples).
    _band: int
        The band number. If part of a band group, this must refer to the lowest band of the group.
    _idx_x: int
        The x-momentum index, with respect to the discretized grid.
    _idx_y: int
        The y-momentum index, with respect to the discretized grid.
    _group_size: int
        The number of touching bands a.k.a. number of bands in the band group (default=1).

    Returns
    -------
    Berry_curv: float
        The Berry curvature around a square plaquette.
    """

    Berry_curv = - np.imag(np.log(U(1, _eigenvectors, _band, _idx_x, _idx_y, _group_size)
                                  * U(2, _eigenvectors, _band, _idx_x+1, _idx_y, _group_size)
                                  * U(1, _eigenvectors, _band, _idx_x, _idx_y+1, _group_size)**-1
                                  * U(2, _eigenvectors, _band, _idx_x, _idx_y, _group_size)**-1))

    return Berry_curv


def wilson_loop(_eigenvectors, _band, _idx_y, _group_size=1):
    r"""
    Compute the Wilson loop.

    The Wilson loop term is defined as the product of Berry phases around a cycle of the Brillouin zone, such that

    .. math::
        W = -\Im \log \prod_{\alpha} \tilde{\mathcal{U}}_2(\mathbf{k}_\alpha),

    where :math:`\tilde{\mathcal{U}}_2` is the normalized link variable, :math:`\mathbf{k}_\alpha` is the discretized momentum vector, and the product is taken on a Brillouin zone cycle in the :math:`\gamma=2` direction. :cite:`Gresch18`

    Parameters
    ----------
    _eigenvectors: ndarray
        The array of eigenvectors with dimension (num_bands, num_bands, num_samples, num_samples).
    _band: int
        The band number. If part of a band group, this must refer to the lowest band of the group.
    _idx_y: int
        The y-momentum index, with respect to the discretized grid.
    _group_size: int
        The number of touching bands a.k.a. number of bands in the band group (default=1).

    Returns
    -------
    Wilson_loop: complex
        The Wilson loop term.
    """

    # print("np.shape(_eigenvectors) = ", np.shape(_eigenvectors))
    numb_kx = np.shape(_eigenvectors)[2]
    product = 1
    for i in range(numb_kx):
        product *= U(1, _eigenvectors, _band, i, _idx_y, _group_size)
    Wilson_loop = -np.imag(np.log(product))

    return Wilson_loop


def geom_tensor(_eigenvectors, _eigenvectors_dkx, _eigenvectors_dky, _bvec, _band, _idx_x, _idx_y, _group_size=1):
    r"""
    Compute the quantum geometric tensor.

    The quantum geometric tensor is computed using

    .. math::
       \mathcal{R}_{\mu\nu}(\mathbf{k}) = \mathrm{tr} ( \mathcal{P}_\mathbf{k} \partial_{k_\mu} \mathcal{P}_\mathbf{k} \partial_{k_\nu} \mathcal{P}_\mathbf{k}),

    where :math:`\mathcal{P}_\mathbf{k} = \sum_n^{N_\mathrm{g}} \ket{u_n(\mathbf{k})}\bra{u_n(\mathbf{k})}` is the band projector and the sum is performed over all bands in the band group :math:`N_\mathrm{g}`. :cite:`Parameswaran13, Hirschmann24`

    Parameters
    ----------
    _eigenvectors: ndarray
        The array of eigenvectors with dimension (num_bands, num_bands, num_samples, num_samples).
    _eigenvectors_dkx: ndarray
        The array of eigenvectors at a dkx offset with dimension (num_bands, num_bands, num_samples, num_samples).
    _eigenvectors_dky: ndarray
        The array of eigenvectors at a dky offset with dimension (num_bands, num_bands, num_samples, num_samples).
    _bvec: ndarray
        The array of reciprocal lattice vectors.
    _band: int
        The band number. If part of a band group, this must refer to the lowest band of the group.
    _idx_x: int
        The x-momentum index, with respect to the discretized grid.
    _idx_y: int
        The y-momentum index, with respect to the discretized grid.
    _group_size: int
        The number of touching bands a.k.a. number of bands in the band group (default=1).

    Returns
    -------
    tensor: ndarray
        The quantum geometric tensor with dimension (2,2).
    """

    numb_samp_x = np.shape(_eigenvectors)[2]
    numb_samp_y = np.shape(_eigenvectors)[3]

    tot_proj, tot_proj_dkx, tot_proj_dky = 0, 0, 0
    for i in range(_group_size):
        tot_proj += np.outer(_eigenvectors[:, _band+i, _idx_x, _idx_y], np.conj(_eigenvectors[:, _band+i, _idx_x, _idx_y]))
        tot_proj_dkx += np.outer(_eigenvectors_dkx[:, _band+i, _idx_x, _idx_y], np.conj(_eigenvectors_dkx[:, _band+i, _idx_x, _idx_y]))
        tot_proj_dky += np.outer(_eigenvectors_dky[:, _band+i, _idx_x, _idx_y], np.conj(_eigenvectors_dky[:, _band+i, _idx_x, _idx_y]))

    dkx = np.linalg.norm(_bvec[0]) / (1000 * (numb_samp_x - 1))
    dky = np.linalg.norm(_bvec[1]) / (1000 * (numb_samp_y - 1))

    grad_kx = np.subtract(tot_proj_dkx, tot_proj) / dkx
    grad_ky = np.subtract(tot_proj_dky, tot_proj) / dky

    tensor = np.zeros((2, 2), dtype=np.complex128)
    tensor[0][0] = np.trace(np.matmul(tot_proj, np.matmul(grad_kx, grad_kx)))
    tensor[0][1] = np.trace(np.matmul(tot_proj, np.matmul(grad_kx, grad_ky)))
    tensor[1][0] = np.trace(np.matmul(tot_proj, np.matmul(grad_ky, grad_kx)))
    tensor[1][1] = np.trace(np.matmul(tot_proj, np.matmul(grad_ky, grad_ky)))

    return tensor
