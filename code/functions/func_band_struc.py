import numpy as np
import cmath as cm


def _principal(log_sum):
    # take the principal branch of a log sum, such that -pi<Im(z)<=pi
    re = np.real(log_sum)
    im = np.imag(log_sum)
    if im <= -np.pi:
        log_sum = re + 1j * (np.pi - np.abs((im + np.pi)) % (2 * np.pi))
    elif im > np.pi:
        log_sum = re + 1j * (-np.pi + np.abs((im + np.pi)) % (2 * np.pi))

    return log_sum


def berry_curv(ev, ev_alpha, ev_beta, ev_alpha_beta):

    # Method from Fukui et al., J. Phys. Soc. Jpn. 74 (2005) pp. 1674-1677
    # Link variables (eigenvectors are already normalized from numpy.linalg.eig)
    U12 = np.conj(ev).dot(ev_alpha)
    U23 = np.conj(ev_alpha).dot(ev_alpha_beta)
    U34 = np.conj(ev_alpha_beta).dot(ev_beta)
    U41 = np.conj(ev_beta).dot(ev)

    Berry_curv = - np.imag(np.log(U12*U23*U34*U41))  # origin of minus sign?
    # Berry_curv = - np.imag(_principal(np.log(U12) + np.log(U23) + np.log(U34) + np.log(U41)))

    return Berry_curv


def berry_curv_2(ev, ev_alpha, ev_beta):

    # Method using imaginary part of quantum geometric tensor
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

    tensor = np.zeros((2, 2), dtype=np.complex128)
    tensor[0][0] = np.conj(ev_alpha).dot(ev_alpha) - np.conj(ev_alpha).dot(ev)*np.conj(ev).dot(ev_alpha)
    tensor[0][1] = np.conj(ev_alpha).dot(ev_beta) - np.conj(ev_alpha).dot(ev)*np.conj(ev).dot(ev_beta)
    tensor[1][0] = np.conj(ev_beta).dot(ev_alpha) - np.conj(ev_beta).dot(ev)*np.conj(ev).dot(ev_alpha)
    tensor[1][1] = np.conj(ev_beta).dot(ev_beta) - np.conj(ev_beta).dot(ev)*np.conj(ev).dot(ev_beta)

    return tensor
