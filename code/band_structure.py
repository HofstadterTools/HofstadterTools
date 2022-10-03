# --- environment imports
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from prettytable import PrettyTable
from tqdm import tqdm
# --- package imports
import functions.func_band_struc as fbs
import functions.func_args as fa
from models.hofstadter import Hofstadter


if __name__ == '__main__':

    # input
    args = fa.parse_input_arguments("band_structure")
    num_samples = args['samp']
    if args['model'] == "Hofstadter":
        model = Hofstadter(args['nphi'][0], args['nphi'][1])
    else:
        raise ValueError("model is not defined")

    # define unit cell
    num_bands, avec, bvec, sym_points = model.unit_cell()

    # construct bands
    eigenvalues = np.zeros((num_bands, num_samples, num_samples))  # real
    eigenvectors = np.zeros((num_bands, num_bands, num_samples, num_samples), dtype=np.complex128)  # complex
    for band in tqdm(range(num_bands), desc="Band Construction", ascii=True):
        for idx_x in range(num_samples):
            frac_kx = idx_x / (num_samples-1)
            for idx_y in range(num_samples):
                frac_ky = idx_y / (num_samples-1)
                k = np.matmul(np.array([frac_kx, frac_ky]), bvec)
                eigvals, eigvecs = np.linalg.eig(model.hamiltonian(k))
                idx = np.argsort(eigvals)
                eigenvalues[band][idx_x][idx_y] = np.real(eigvals[idx[band]])
                eigenvectors[:, band, idx_x, idx_y] = eigvecs[:, idx[band]]

    # band gap and isolated
    band_gap = np.zeros(num_bands)
    isolated = np.full(num_bands, True)
    for band_idx, band in enumerate(np.arange(num_bands)[::-1]):
        if band_idx == 0:
            band_gap[band] = "NaN"
        else:
            band_gap[band] = np.min(eigenvalues[band + 1]) - np.max(eigenvalues[band])
        if band_gap[band] < 1e-10:
            isolated[band] = False
            isolated[band + 1] = False

    # compute Berry fluxes
    berry_fluxes = np.zeros((num_bands, num_samples - 1, num_samples - 1))  # real
    TISM = np.zeros((num_bands, num_samples - 1, num_samples - 1))  # real
    DISM = np.zeros((num_bands, num_samples - 1, num_samples - 1))  # real
    for band in tqdm(range(num_bands), desc="Band Properties", ascii=True):
        for idx_x in range(num_samples - 1):
            for idx_y in range(num_samples - 1):
                if isolated[band]:
                    berry_fluxes[band, idx_x, idx_y] = fbs.berry_curv(eigenvectors[:, band, idx_x, idx_y],
                                                                      eigenvectors[:, band, idx_x + 1, idx_y],
                                                                      eigenvectors[:, band, idx_x, idx_y + 1],
                                                                      eigenvectors[:, band, idx_x + 1, idx_y + 1])
                    TISM[band, idx_x, idx_y] = np.trace(fbs.fubini_study(eigenvectors[:, band, idx_x, idx_y],
                                                                         eigenvectors[:, band, idx_x + 1, idx_y],
                                                                         eigenvectors[:, band, idx_x, idx_y + 1])) \
                                               - np.abs(berry_fluxes[band, idx_x, idx_y])
                    DISM[band, idx_x, idx_y] = np.linalg.det(fbs.fubini_study(eigenvectors[:, band, idx_x, idx_y],
                                                                              eigenvectors[:, band, idx_x + 1, idx_y],
                                                                              eigenvectors[:, band, idx_x, idx_y + 1])) \
                                               - 0.25*np.abs(berry_fluxes[band, idx_x, idx_y])**2
                elif band == 0 or isolated[band-1]:
                    berry_fluxes[band, idx_x, idx_y] = fbs.multi_berry_curv(eigenvectors[:, band, idx_x, idx_y],
                                                                            eigenvectors[:, band, idx_x+1, idx_y],
                                                                            eigenvectors[:, band, idx_x, idx_y+1],
                                                                            eigenvectors[:, band, idx_x+1, idx_y+1],
                                                                            eigenvectors[:, band+1, idx_x, idx_y],
                                                                            eigenvectors[:, band+1, idx_x+1, idx_y],
                                                                            eigenvectors[:, band+1, idx_x, idx_y+1],
                                                                            eigenvectors[:, band+1, idx_x+1, idx_y+1])
                else:
                    berry_fluxes[band, idx_x, idx_y] = berry_fluxes[band-1, idx_x, idx_y]

    # band properties
    chern_numbers = np.zeros(num_bands)
    berry_fluc = np.zeros(num_bands)
    band_width = np.zeros(num_bands)
    TISM_average = np.zeros(num_bands)
    DISM_average = np.zeros(num_bands)
    for band_idx, band in enumerate(np.arange(num_bands)[::-1]):
        chern_numbers[band] = np.sum(berry_fluxes[band, :, :]) / (2 * np.pi)
        berry_fluc[band] = np.std(berry_fluxes[band, :, :])/np.abs(np.average(berry_fluxes[band, :, :]))
        band_width[band] = np.max(eigenvalues[band]) - np.min(eigenvalues[band])
        TISM_average[band] = np.mean(TISM[band])
        DISM_average[band] = np.mean(DISM[band])

    # table
    table = PrettyTable()
    table.field_names = ["band", "isolated", "C", "sigma_B/|mu_B|", "width", "gap", "gap/width", "<T>", "<D>"]
    for band in np.arange(num_bands)[::-1]:
        table.add_row([band, isolated[band], round(chern_numbers[band]), berry_fluc[band], band_width[band],
                       band_gap[band], band_gap[band]/band_width[band], TISM_average[band], DISM_average[band]])
    print(table)

    # construct figure
    if args['display'] == "3D":
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        idx_x = np.linspace(0, num_samples - 1, num_samples, dtype=int)
        idx_y = np.linspace(0, num_samples - 1, num_samples, dtype=int)
        kx, ky = np.meshgrid(idx_x, idx_y)
        for i in range(num_bands):
            ax.plot_surface(kx, ky, eigenvalues[i, kx, ky])
        ax.set_xlabel('$k_x$')
        ax.set_ylabel('$k_y$')
        ax.set_zlabel('$E$')
    elif args['display'] == "2D":
        # construct bands
        num_paths = len(sym_points)
        points_per_path = int(num_samples/num_paths)
        num_points = num_paths*points_per_path
        eigenvalues = np.zeros((num_bands, num_points))  # real
        count = 0
        for i in range(num_paths):
            for j in range(points_per_path):
                k = sym_points[i] + (sym_points[(i+1) % num_paths] - sym_points[i]) * float(j) / float(points_per_path - 1)
                k = np.matmul(k, bvec)
                eigvals = np.linalg.eigvals(model.hamiltonian(k))
                idx = np.argsort(eigvals)
                for band in range(num_bands):
                    eigenvalues[band, count] = np.real(eigvals[idx[band]])
                count += 1

        # construct figure
        fig = plt.figure()
        ax = plt.subplot(111)
        for i in range(num_bands):
            ax.plot(eigenvalues[i])
        for i in range(1, num_paths):
            ax.axvline(i*points_per_path, color='k', linewidth=0.5, ls='--')
        ax.set_xlim([0, num_points])
        ax.set_xlabel('symmetry path')
        ax.set_ylabel('$E$')

    plt.show()
