import numpy as np


def berry_curv(eigvec, eigvec_alpha, eigvec_beta, eigvec_alpha_beta):

    Berry_curv = - np.imag(np.log(np.conj(eigvec).dot(eigvec_alpha) * np.conj(eigvec_alpha).dot(eigvec_alpha_beta)
                           * np.conj(eigvec_alpha_beta).dot(eigvec_beta) * np.conj(eigvec_beta).dot(eigvec)))

    return Berry_curv


def explicit_berry_curv(eigvec, eigvec_alpha, eigvec_beta, eigvec_alpha_beta):

    term1 = -np.angle(np.conj(eigvec).dot(eigvec_alpha))
    term2 = -np.angle(np.conj(eigvec_alpha).dot(eigvec_alpha_beta))

    term3 = -np.angle(np.conj(eigvec_alpha_beta).dot(eigvec_beta))
    term4 = -np.angle(np.conj(eigvec_beta).dot(eigvec))

    Berry_curv = term1 + term2 + term3 + term4

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


def fubini_study(ev1, ev1_alpha, ev1_beta):

    FS_matrix = np.zeros((2, 2))
    FS_matrix[0][0] = np.real(np.conj(ev1_alpha).dot(ev1_alpha) - np.conj(ev1_alpha).dot(ev1)*np.conj(ev1).dot(ev1_alpha))
    FS_matrix[0][1] = np.real(0.5 * (np.conj(ev1_alpha).dot(ev1_beta) + np.conj(ev1_beta).dot(ev1_alpha)
                              - np.conj(ev1_alpha).dot(ev1)*np.conj(ev1).dot(ev1_beta) - np.conj(ev1_beta).dot(ev1)*np.conj(ev1).dot(ev1_alpha)))
    FS_matrix[1][0] = np.real(0.5 * (np.conj(ev1_beta).dot(ev1_alpha) + np.conj(ev1_alpha).dot(ev1_beta)
                              - np.conj(ev1_beta).dot(ev1)*np.conj(ev1).dot(ev1_alpha) - np.conj(ev1_alpha).dot(ev1)*np.conj(ev1).dot(ev1_beta)))
    FS_matrix[1][1] = np.real(np.conj(ev1_beta).dot(ev1_beta) - np.conj(ev1_beta).dot(ev1)*np.conj(ev1).dot(ev1_beta))

    return FS_matrix
