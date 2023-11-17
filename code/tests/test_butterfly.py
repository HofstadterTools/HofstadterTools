import numpy as np
from math import gcd
from models.hofstadter import Hofstadter


def butterfly(q, t, lat, alpha=1, theta=(1, 3), period=1):

    nphi_list, E_list = [], []

    for p in range(1, q):

        # construct model
        model = Hofstadter(p, q, t=t, lat=lat, alpha=alpha, theta=theta, period=period)

        # define flux density
        if gcd(p, q) != 1:  # nphi must be a coprime fraction
            continue
        nphi = p / q

        # diagonalize Hamiltonian
        ham = model.hamiltonian(np.array([0, 0]))
        M = len(ham)
        nphi_list.append([nphi] * M)
        lmbda = np.sort(np.linalg.eigvalsh(ham))
        E_list.append(lmbda)

    return nphi_list, E_list


def test_square():

    # current
    nphi_list, E_list = butterfly(199, [1], "square")
    current = np.array([nphi_list, E_list])
    # reference
    filename = "code/tests/butterfly/butterfly_square_q_199_t_1.npy"
    file_data = np.load(filename, allow_pickle=True)
    data = file_data[2]
    nphi_list_ref = data['nphi_list']
    E_list_ref = data['E_list']
    reference = np.array([nphi_list_ref, E_list_ref])

    assert np.allclose(current, reference)


def test_triangular():

    # current
    nphi_list, E_list = butterfly(199, [1], "triangular")
    current = np.array([nphi_list, E_list])
    # reference
    filename = "code/tests/butterfly/butterfly_triangular_q_199_t_1.npy"
    file_data = np.load(filename, allow_pickle=True)
    data = file_data[2]
    nphi_list_ref = data['nphi_list']
    E_list_ref = data['E_list']
    reference = np.array([nphi_list_ref, E_list_ref])

    assert np.allclose(current, reference)


def test_bravais():

    # current
    nphi_list, E_list = butterfly(199, [0.5, 0.2], "bravais", alpha=1, theta=(67, 180))
    current = np.array([nphi_list, E_list])
    # reference
    filename = "code/tests/butterfly/butterfly_bravais_q_199_t_0.5_0.2_alpha_1_theta_67_180.npy"
    file_data = np.load(filename, allow_pickle=True)
    data = file_data[2]
    nphi_list_ref = data['nphi_list']
    E_list_ref = data['E_list']
    reference = np.array([nphi_list_ref, E_list_ref])

    assert np.allclose(current, reference)


def test_honeycomb():

    # current
    nphi_list, E_list = butterfly(199, [1], "honeycomb")
    current = np.array([nphi_list, E_list])
    # reference
    filename = "code/tests/butterfly/butterfly_honeycomb_q_199_t_1_alpha_1_theta_1_3.npy"
    file_data = np.load(filename, allow_pickle=True)
    data = file_data[2]
    nphi_list_ref = data['nphi_list']
    E_list_ref = data['E_list']
    reference = np.array([nphi_list_ref, E_list_ref])

    assert np.allclose(current, reference)


def test_kagome():

    # current
    nphi_list, E_list = butterfly(199, [1], "kagome", alpha=1, theta=(1, 3), period=8)
    current = np.array([nphi_list, E_list])
    # reference
    filename = "code/tests/butterfly/butterfly_kagome_q_199_t_1_alpha_1_theta_1_3_period_8.npy"
    file_data = np.load(filename, allow_pickle=True)
    data = file_data[2]
    nphi_list_ref = data['nphi_list']
    E_list_ref = data['E_list']
    reference = np.array([nphi_list_ref, E_list_ref])

    assert np.allclose(current, reference)
