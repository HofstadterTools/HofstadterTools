"""Unit tests for the butterfly program."""

# --- external imports
import numpy as np
from math import gcd
# --- internal imports
from HT.models.hofstadter import Hofstadter


def butterfly(q, t, lat, alpha=1, theta=(1, 3), period=1):
    """Minimal function to generate the butterfly spectrum.

    Parameters
    ----------
    q: int
        The denominator of the coprime flux density fraction.
    t: list
        The list of hopping amplitudes in order of ascending NN.
    lat: str
        The name of the lattice.
    alpha: float
        The anisotropy of the Bravais lattice vectors (default=1).
    theta: tuple
        The angle between Bravais lattice vectors in units of pi (default=(1, 3)).
    period: int
        The factor by which to divide A_UC in the flux density (default=1).

    Returns
    -------
    nphi_list: list
        The list of flux densities.
    E_list: list
        The list of corresponding energies.
    """

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


def test_square(test_dir):
    """Benchmark the square Hofstadter butterfly. This butterfly is plotted in many papers, e.g. Fig.4(d) of :cite:`Yilmaz17`."""

    # current
    nphi_list, E_list = butterfly(97, [1], "square")
    current = np.array([nphi_list, E_list])
    # reference
    filename = test_dir.joinpath('butterfly', 'butterfly_square_q_97_t_1.npz')
    file_data = np.load(filename, allow_pickle=True)
    data = file_data['data'].item()
    nphi_list_ref = data['nphi_list']
    E_list_ref = data['E_list']
    reference = np.array([nphi_list_ref, E_list_ref])

    assert np.allclose(current, reference)


def test_square_2(test_dir):
    """Benchmark the zero-quadratic model butterfly against Fig.2 of :cite:`Bauer22`."""

    # current
    nphi_list, E_list = butterfly(97, [1, 0, -0.25], "square")
    current = np.array([nphi_list, E_list])
    # reference
    filename = test_dir.joinpath('butterfly', 'butterfly_square_q_97_t_1_0_-0.25.npz')
    file_data = np.load(filename, allow_pickle=True)
    data = file_data['data'].item()
    nphi_list_ref = data['nphi_list']
    E_list_ref = data['E_list']
    reference = np.array([nphi_list_ref, E_list_ref])

    assert np.allclose(current, reference)


def test_triangular(test_dir):
    """Benchmark the triangular Hofstadter butterfly. This butterfly is plotted in many papers, e.g. Fig.4(a) of :cite:`Yilmaz17`."""

    # current
    nphi_list, E_list = butterfly(97, [1], "triangular")
    current = np.array([nphi_list, E_list])
    # reference
    filename = test_dir.joinpath('butterfly', 'butterfly_triangular_q_97_t_1.npz')
    file_data = np.load(filename, allow_pickle=True)
    data = file_data['data'].item()
    nphi_list_ref = data['nphi_list']
    E_list_ref = data['E_list']
    reference = np.array([nphi_list_ref, E_list_ref])

    assert np.allclose(current, reference)


def test_triangular_2(test_dir):
    """Benchmark the 2nd-NN triangular Hofstadter butterfly. This butterfly is plotted in many papers, e.g. Fig.4 of :cite:`Oh00`."""

    # current
    nphi_list, E_list = butterfly(97, [0, 1], "triangular", period=2)
    current = np.array([nphi_list, E_list])
    # reference
    filename = test_dir.joinpath('butterfly', 'butterfly_triangular_q_97_t_0_1_period_2.npz')
    file_data = np.load(filename, allow_pickle=True)
    data = file_data['data'].item()
    nphi_list_ref = data['nphi_list']
    E_list_ref = data['E_list']
    reference = np.array([nphi_list_ref, E_list_ref])

    assert np.allclose(current, reference)


def test_bravais(test_dir):
    """Benchmark the Bravais Hofstadter butterfly against Fig.4(b-c) of :cite:`Yilmaz17`."""

    # current
    nphi_list, E_list = butterfly(97, [0.5, 0.2], "bravais", alpha=1, theta=(67, 180))
    current = np.array([nphi_list, E_list])
    # reference
    filename = test_dir.joinpath('butterfly', 'butterfly_bravais_q_97_t_0.5_0.2_alpha_1_theta_67_180.npz')
    file_data = np.load(filename, allow_pickle=True)
    data = file_data['data'].item()
    nphi_list_ref = data['nphi_list']
    E_list_ref = data['E_list']
    reference = np.array([nphi_list_ref, E_list_ref])

    assert np.allclose(current, reference)


def test_bravais_2(test_dir):
    """Benchmark the Bravais Hofstadter butterfly in the triangular lattice limit."""

    # current
    nphi_list, E_list = butterfly(97, [1], "bravais", alpha=1, theta=(1, 3))
    current = np.array([nphi_list, E_list])
    # reference
    filename = test_dir.joinpath('butterfly', 'butterfly_triangular_q_97_t_1.npz')
    file_data = np.load(filename, allow_pickle=True)
    data = file_data['data'].item()
    nphi_list_ref = data['nphi_list']
    E_list_ref = data['E_list']
    reference = np.array([nphi_list_ref, E_list_ref])

    assert np.allclose(current, reference)


def test_bravais_3(test_dir):
    """Benchmark the Bravais Hofstadter butterfly in the square lattice limit."""

    # current
    nphi_list, E_list = butterfly(97, [1], "bravais", alpha=1, theta=(1, 2))
    current = np.array([nphi_list, E_list])
    # reference
    filename = test_dir.joinpath('butterfly', 'butterfly_square_q_97_t_1.npz')
    file_data = np.load(filename, allow_pickle=True)
    data = file_data['data'].item()
    nphi_list_ref = data['nphi_list']
    E_list_ref = data['E_list']
    reference = np.array([nphi_list_ref, E_list_ref])

    assert np.allclose(current, reference)


def test_honeycomb(test_dir):
    """Benchmark the honeycomb Hofstadter butterfly. This butterfly is plotted in many papers, e.g. Fig.5 of :cite:`Agazzi14`."""

    # current
    nphi_list, E_list = butterfly(97, [1], "honeycomb")
    current = np.array([nphi_list, E_list])
    # reference
    filename = test_dir.joinpath('butterfly', 'butterfly_honeycomb_q_97_t_1_alpha_1_theta_1_3.npz')
    file_data = np.load(filename, allow_pickle=True)
    data = file_data['data'].item()
    nphi_list_ref = data['nphi_list']
    E_list_ref = data['E_list']
    reference = np.array([nphi_list_ref, E_list_ref])

    assert np.allclose(current, reference)


def test_honeycomb_2(test_dir):
    """Benchmark the 2nd-NN honeycomb Hofstadter butterfly. This is equivalent to the triangular Hofstadter butterfly."""

    # current
    nphi_list, E_list = butterfly(97, [0, 1], "honeycomb")
    current = np.array([nphi_list, E_list])
    # reference
    filename = test_dir.joinpath('butterfly', 'butterfly_triangular_q_97_t_1.npz')
    file_data = np.load(filename, allow_pickle=True)
    data = file_data['data'].item()
    nphi_list_ref = data['nphi_list']
    E_list_ref = data['E_list']
    reference = np.array([nphi_list_ref, E_list_ref])

    assert np.allclose(current, reference)


def test_kagome(test_dir):
    """Benchmark the kagome Hofstadter butterfly. This butterfly is plotted in many papers, e.g. Fig.3 of :cite:`Jing-Min09`."""

    # current
    nphi_list, E_list = butterfly(97, [1], "kagome", alpha=1, theta=(1, 3), period=8)
    current = np.array([nphi_list, E_list])
    # reference
    filename = test_dir.joinpath('butterfly', 'butterfly_kagome_q_97_t_1_alpha_1_theta_1_3_period_8.npz')
    file_data = np.load(filename, allow_pickle=True)
    data = file_data['data'].item()
    nphi_list_ref = data['nphi_list']
    E_list_ref = data['E_list']
    reference = np.array([nphi_list_ref, E_list_ref])

    assert np.allclose(current, reference)
