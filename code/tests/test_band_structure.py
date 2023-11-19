import numpy as np
import functions.models as fm
from models.hofstadter import Hofstadter


def band_structure(nphi, t, lat, alpha=1, theta=(1, 3), period=1, samp=11):
    """Minimal function to generate the band structure.

    Parameters
    ----------
    nphi: tuple
        The flux density.
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
    samp: int
        The number of samples along an axis of the momentum grid (default=11).

    Returns
    -------
    eigenvalues: ndarray
        The array of eigenvalues.
    eigenvectors: ndarray
        The array of eigenvectors.
    """

    # define model
    model = Hofstadter(nphi[0], nphi[1], t=t, lat=lat, alpha=alpha, theta=theta, period=period)

    # define unit cell
    _, avec, abasisvec, bMUCvec, sym_points = model.unit_cell()
    _, bases = fm.nearest_neighbor_finder(avec, abasisvec, t, 0, 0, 0)
    num_bands = nphi[1] * len(bases)

    # construct bands
    eigenvalues = np.zeros((num_bands, samp, samp))  # real
    eigenvectors = np.zeros((num_bands, num_bands, samp, samp), dtype=np.complex128)  # complex
    for band in range(num_bands):
        for idx_x in range(samp):
            frac_kx = idx_x / (samp - 1)
            for idx_y in range(samp):
                frac_ky = idx_y / (samp - 1)
                k = np.matmul(np.array([frac_kx, frac_ky]), bMUCvec)
                ham = model.hamiltonian(k)
                eigvals, eigvecs = np.linalg.eigh(ham)
                idx = np.argsort(eigvals)
                eigenvalues[band, idx_x, idx_y] = eigvals[idx[band]]
                eigenvectors[:, band, idx_x, idx_y] = eigvecs[:, idx[band]]

    return eigenvalues, eigenvectors


def test_square():

    # current
    eigenvalues, eigenvectors = band_structure((1, 5), [1], "square")
    # reference
    filename = "code/tests/band_structure/band_structure_3D_square_nphi_1_5_t_1_samp_11.npz"
    file_data = np.load(filename, allow_pickle=True)
    data = file_data['data'].item()
    eigenvalues_ref = data['eigenvalues']
    eigenvectors_ref = data['eigenvectors']

    assert np.allclose(eigenvalues, eigenvalues_ref)
    assert np.allclose(np.abs(eigenvectors), np.abs(eigenvectors_ref))


def test_triangular():

    # current
    eigenvalues, eigenvectors = band_structure((1, 5), [1], "triangular")
    # reference
    filename = "code/tests/band_structure/band_structure_3D_triangular_nphi_1_5_t_1_samp_11.npz"
    file_data = np.load(filename, allow_pickle=True)
    data = file_data['data'].item()
    eigenvalues_ref = data['eigenvalues']
    eigenvectors_ref = data['eigenvectors']

    assert np.allclose(eigenvalues, eigenvalues_ref)
    assert np.allclose(np.abs(eigenvectors), np.abs(eigenvectors_ref))


def test_bravais():

    # current
    eigenvalues, eigenvectors = band_structure((1, 5), [0.5, 0.2], "bravais", alpha=1, theta=(67, 180))
    # reference
    filename = "code/tests/band_structure/band_structure_3D_bravais_nphi_1_5_t_0.5_0.2_alpha_1_theta_67_180_samp_11.npz"
    file_data = np.load(filename, allow_pickle=True)
    data = file_data['data'].item()
    eigenvalues_ref = data['eigenvalues']
    eigenvectors_ref = data['eigenvectors']

    assert np.allclose(eigenvalues, eigenvalues_ref)
    assert np.allclose(np.abs(eigenvectors), np.abs(eigenvectors_ref))


def test_honeycomb():

    # current
    eigenvalues, eigenvectors = band_structure((1, 5), [1], "honeycomb")
    # reference
    filename = "code/tests/band_structure/band_structure_3D_honeycomb_nphi_1_5_t_1_alpha_1_theta_1_3_samp_11.npz"
    file_data = np.load(filename, allow_pickle=True)
    data = file_data['data'].item()
    eigenvalues_ref = data['eigenvalues']
    eigenvectors_ref = data['eigenvectors']

    assert np.allclose(eigenvalues, eigenvalues_ref)
    assert np.allclose(np.abs(eigenvectors), np.abs(eigenvectors_ref))


def test_kagome():

    # current
    eigenvalues, eigenvectors = band_structure((1, 5), [1], "kagome", period=8)
    # reference
    filename = "code/tests/band_structure/band_structure_3D_kagome_nphi_1_5_t_1_alpha_1_theta_1_3_samp_11.npz"
    file_data = np.load(filename, allow_pickle=True)
    data = file_data['data'].item()
    eigenvalues_ref = data['eigenvalues']
    eigenvectors_ref = data['eigenvectors']

    assert np.allclose(eigenvalues, eigenvalues_ref)
    assert np.allclose(np.abs(eigenvectors), np.abs(eigenvectors_ref))
