"""Set of functions for band structure calculations."""

import numpy as np
import cmath as cm


def _principal(log_sum):
    """
    Compute the principal branch of a complex number :math:`z`, such that :math:`-\pi<\text{Im}(z)\leq\pi`.

    Parameters
    ----------
    log_sum: complex
        The original complex number.

    Returns
    -------
    log_sum: complex
        The principal value of the complex number.
    """

    re = np.real(log_sum)
    im = np.imag(log_sum)
    if im <= -np.pi:
        log_sum = re + 1j * (np.pi - np.abs((im + np.pi)) % (2 * np.pi))
    elif im > np.pi:
        log_sum = re + 1j * (-np.pi + np.abs((im + np.pi)) % (2 * np.pi))

    return log_sum


def berry_curv(ev, ev_alpha, ev_beta, ev_alpha_beta):
    r"""
    Compute the Berry curvature around a plaquette, for a single band.

    The Berry curvature around a square plaquette is computed using the formula from :cite:`Fukui05` (example applications in :cite:`SoluyanovPhD, AidelsburgerPhD`), such that

    .. math::
       \mathcal{B}_{\square}(\mathbf{k}) = - \text{Im}\;\log (\braket{u_{\mathbf{k}} | u_{\mathbf{k}+\Delta_\alpha}} \braket{u_{\mathbf{k}+\Delta_\alpha}|u_{\mathbf{k}+\Delta_\alpha+\Delta_\beta}} \braket{u_{\mathbf{k}+\Delta_\alpha+\Delta_\beta}|u_{\mathbf{k}+\Delta_\beta}} \braket{u_{\mathbf{k}+\Delta_\beta}|u_{\mathbf{k}}}),

    where :math:`\ket{u_{\mathbf{k}}}` is the eigenvector at momentum :math:`\mathbf{k}` and :math:`\Delta_\alpha` is shorthand for :math:`\Delta \mathbf{k}_\alpha`. The Berry curvature at the point :math:`\mathbf{k}` can then be computed by taking the limit of small plaquette size.

    .. note::
        Input eigenvectors are already normalized from :class:`numpy.linalg.eig`.

    .. note::
        We may equivalently compute the Berry curvature around a plaquette using the log sum formula

        .. math::
            \mathcal{B}_{\square}(\mathbf{k}) = - \text{Im}\;\mathcal{P}\;(\log \braket{u_{\mathbf{k}} | u_{\mathbf{k}+\Delta_\alpha}} + \log\braket{u_{\mathbf{k}+\Delta_\alpha} | u_{\mathbf{k}+\Delta_\alpha+\Delta_\beta}} + \log\braket{u_{\mathbf{k}+\Delta_\alpha+\Delta_\beta} | u_{\mathbf{k}+\Delta_\beta}} + \log\braket{u_{\mathbf{k}+\Delta_\beta} | u_{\mathbf{k}}}),

        where :math:`\mathcal{P}` denotes the principal value of the complex number :math:`z`, such that :math:`-\pi<\text{Im}(z)\leq\pi`.

    Parameters
    ----------
    ev: ndarray
        The normalized eigenvector :math:`\ket{u_{\mathbf{k}}}`.
    ev_alpha: ndarray
        The normalized eigenvector :math:`\ket{u_{\mathbf{k}+\Delta_\alpha}}`
    ev_beta: ndarray
        The normalized eigenvector :math:`\ket{u_{\mathbf{k}+\Delta_\beta}}`.
    ev_alpha_beta: ndarray
        The normalized eigenvector :math:`\ket{u_{\mathbf{k}+\Delta_\alpha+\Delta_\beta}}`

    Returns
    -------
    Berry_curv: float
        The Berry curvature around a square plaquette.
    """

    U12 = np.conj(ev).dot(ev_alpha)
    U23 = np.conj(ev_alpha).dot(ev_alpha_beta)
    U34 = np.conj(ev_alpha_beta).dot(ev_beta)
    U41 = np.conj(ev_beta).dot(ev)

    Berry_curv = - np.imag(np.log(U12*U23*U34*U41))  # origin of minus sign?
    # Berry_curv = - np.imag(_principal(np.log(U12) + np.log(U23) + np.log(U34) + np.log(U41)))

    return Berry_curv


def berry_curv_2(ev, ev_alpha, ev_beta):
    r"""
    Compute the Berry curvature around a plaquette, for a single band, using the quantum metric.

    The Berry curvature around a square plaquette is computed from the quantum metric :cite:`Parameswaran13` (example applications in :cite:`Claassen15`), such that

    .. math::
       \mathcal{B}_{\square}(\mathbf{k}) = - 2 \text{Im}[\mathcal{R}_{\alpha \beta}(\mathbf{k})],

    where the quantum metric tensor is defined as

    .. math::
       \mathcal{R}_{\alpha \beta} = \braket{u_{\mathbf{k}+\Delta_\alpha}|u_{\mathbf{k}+\Delta_\beta}} - \braket{u_{\mathbf{k}+\Delta_\alpha}|u_\mathbf{k}}\braket{u_\mathbf{k}|u_{\mathbf{k}+\Delta_\beta}},

    where :math:`\ket{u_{\mathbf{k}}}` is the eigenvector at momentum :math:`\mathbf{k}` and :math:`\Delta_\alpha` is shorthand for :math:`\Delta \mathbf{k}_\alpha`. The Berry curvature at the point :math:`\mathbf{k}` can then be computed by taking the limit of small plaquette size.

    .. note::
        This method converges significantly slower than :class:`berry_curv`.

    Parameters
    ----------
    ev: ndarray
        The normalized eigenvector :math:`\ket{u_{\mathbf{k}}}`.
    ev_alpha: ndarray
        The normalized eigenvector :math:`\ket{u_{\mathbf{k}+\Delta_\alpha}}`
    ev_beta: ndarray
        The normalized eigenvector :math:`\ket{u_{\mathbf{k}+\Delta_\beta}}`.

    Returns
    -------
    Berry_curv: float
        The Berry curvature around a square plaquette.
    """

    chi = np.conj(ev_alpha).dot(ev_beta) - np.conj(ev_alpha).dot(ev)*np.conj(ev).dot(ev_beta)
    Berry_curv = -2 * np.imag(chi)

    return Berry_curv


def multi_berry_curv(ev1, ev1_alpha, ev1_beta, ev1_alpha_beta, ev2, ev2_alpha, ev2_beta, ev2_alpha_beta):

    matrix1 = np.zeros((2, 2), dtype=np.complex128)
    matrix1[0][0] = np.conj(ev1).dot(ev1_alpha)
    matrix1[0][1] = np.conj(ev1).dot(ev2_alpha)
    matrix1[1][0] = np.conj(ev2).dot(ev1_alpha)
    matrix1[1][1] = np.conj(ev2).dot(ev2_alpha)

    matrix2 = np.zeros((2, 2), dtype=np.complex128)
    matrix2[0][0] = np.conj(ev1_alpha).dot(ev1_alpha_beta)
    matrix2[0][1] = np.conj(ev1_alpha).dot(ev2_alpha_beta)
    matrix2[1][0] = np.conj(ev2_alpha).dot(ev1_alpha_beta)
    matrix2[1][1] = np.conj(ev2_alpha).dot(ev2_alpha_beta)

    matrix3 = np.zeros((2, 2), dtype=np.complex128)
    matrix3[0][0] = np.conj(ev1_alpha_beta).dot(ev1_beta)
    matrix3[0][1] = np.conj(ev1_alpha_beta).dot(ev2_beta)
    matrix3[1][0] = np.conj(ev2_alpha_beta).dot(ev1_beta)
    matrix3[1][1] = np.conj(ev2_alpha_beta).dot(ev2_beta)

    matrix4 = np.zeros((2, 2), dtype=np.complex128)
    matrix4[0][0] = np.conj(ev1_beta).dot(ev1)
    matrix4[0][1] = np.conj(ev1_beta).dot(ev2)
    matrix4[1][0] = np.conj(ev2_beta).dot(ev1)
    matrix4[1][1] = np.conj(ev2_beta).dot(ev2)

    multi_bc = - np.imag(np.log(np.linalg.det(matrix1) * np.linalg.det(matrix2) * np.linalg.det(matrix3) * np.linalg.det(matrix4)))

    return multi_bc


def geom_tensor(ev, ev_alpha, ev_beta):
    r"""
    Compute the quantum geometric tensor, for a single band.

    The quantum geometric tensor is computed using the formula from :cite:`Parameswaran13` (example applications in :cite:`Claassen15`), such that

    .. math::
       \mathcal{R}_{\alpha \beta} = \braket{u_{\mathbf{k}+\Delta_\alpha}|u_{\mathbf{k}+\Delta_\beta}} - \braket{u_{\mathbf{k}+\Delta_\alpha}|u_\mathbf{k}}\braket{u_\mathbf{k}|u_{\mathbf{k}+\Delta_\beta}},

    where :math:`\ket{u_{\mathbf{k}}}` is the eigenvector at momentum :math:`\mathbf{k}` and :math:`\Delta_\alpha` is shorthand for :math:`\Delta \mathbf{k}_\alpha`. The quantum geometric tensor at the point :math:`\mathbf{k}` can then be computed by taking the limit of small discretization.

    Parameters
    ----------
    ev: ndarray
        The normalized eigenvector :math:`\ket{u_{\mathbf{k}}}`.
    ev_alpha: ndarray
        The normalized eigenvector :math:`\ket{u_{\mathbf{k}+\Delta_\alpha}}`
    ev_beta: ndarray
        The normalized eigenvector :math:`\ket{u_{\mathbf{k}+\Delta_\beta}}`.

    Returns
    -------
    tensor: ndarray
        The (2,2) quantum geometric tensor.
    """

    tensor = np.zeros((2, 2), dtype=np.complex128)
    tensor[0][0] = np.conj(ev_alpha).dot(ev_alpha) - np.conj(ev_alpha).dot(ev)*np.conj(ev).dot(ev_alpha)
    tensor[0][1] = np.conj(ev_alpha).dot(ev_beta) - np.conj(ev_alpha).dot(ev)*np.conj(ev).dot(ev_beta)
    tensor[1][0] = np.conj(ev_beta).dot(ev_alpha) - np.conj(ev_beta).dot(ev)*np.conj(ev).dot(ev_alpha)
    tensor[1][1] = np.conj(ev_beta).dot(ev_beta) - np.conj(ev_beta).dot(ev)*np.conj(ev).dot(ev_beta)

    return tensor
