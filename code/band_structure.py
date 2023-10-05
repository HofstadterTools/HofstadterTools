# --- external imports
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from prettytable import PrettyTable
from tqdm import tqdm
# --- internal imports
import functions.band_structure as fbs
import functions.arguments as fa
from models.hofstadter import Hofstadter
from configuration.band_structure import *

# plt.rc('text', usetex=True)
# plt.rc('text.latex', preamble=r'\usepackage{amsmath}')

if __name__ == '__main__':

    # input arguments
    args = fa.parse_input_arguments("band_structure")
    t = args['t']
    lat = args['lattice']
    if args['model'] == "Hofstadter":
        model = Hofstadter(args['nphi'][0], args['nphi'][1], t=t, lat=lat)
    else:
        raise ValueError("model is not defined")
    num_samples = args['samp']
    band_gap_threshold = args['bgt']
    wilson = args['wilson']

    # define unit cell
    num_bands, _, _, aMUCvec, bMUCvec, sym_points = model.unit_cell()

    # construct bands
    eigenvalues = np.zeros((num_bands, num_samples, num_samples))  # real
    eigenvectors = np.zeros((num_bands, num_bands, num_samples, num_samples), dtype=np.complex128)  # complex
    if any(geometry_columns):
        eigenvectors_dkx = np.zeros((num_bands, num_bands, num_samples, num_samples), dtype=np.complex128)  # complex
        eigenvectors_dky = np.zeros((num_bands, num_bands, num_samples, num_samples), dtype=np.complex128)  # complex
        Dkx = np.dot(np.array([1/(num_samples-1), 0]), bMUCvec[0])
        Dky = np.dot(np.array([0, 1/(num_samples-1)]), bMUCvec[1])
    for band in tqdm(range(num_bands), desc="Band Construction", ascii=True):
        for idx_x in range(num_samples):
            frac_kx = idx_x / (num_samples-1)
            if any(geometry_columns):
                frac_kx_dkx = (frac_kx + 1/(1000*(num_samples-1))) % 1
            for idx_y in range(num_samples):
                frac_ky = idx_y / (num_samples-1)
                k = np.matmul(np.array([frac_kx, frac_ky]), bMUCvec)
                eigvals, eigvecs = np.linalg.eigh(model.hamiltonian(k))
                idx = np.argsort(eigvals)
                eigenvalues[band, idx_x, idx_y] = eigvals[idx[band]]
                eigenvectors[:, band, idx_x, idx_y] = eigvecs[:, idx[band]]
                if any(geometry_columns):
                    frac_ky_dky = (frac_ky + 1/(1000*(num_samples-1))) % 1
                    k_dkx = np.matmul(np.array([frac_kx_dkx, frac_ky]), bMUCvec)
                    k_dky = np.matmul(np.array([frac_kx, frac_ky_dky]), bMUCvec)
                    eigvals_dkx, eigvecs_dkx = np.linalg.eig(model.hamiltonian(k_dkx))
                    eigvals_dky, eigvecs_dky = np.linalg.eig(model.hamiltonian(k_dky))
                    idx_dkx = np.argsort(eigvals_dkx)
                    idx_dky = np.argsort(eigvals_dky)
                    eigenvectors_dkx[:, band, idx_x, idx_y] = eigvecs_dkx[:, idx_dkx[band]]
                    eigenvectors_dky[:, band, idx_x, idx_y] = eigvecs_dky[:, idx_dky[band]]

    # band gap and isolated
    band_gap = np.zeros(num_bands)
    isolated = np.full(num_bands, True)
    for band_idx, band in enumerate(np.arange(num_bands)[::-1]):
        if band_idx == 0:
            band_gap[band] = "NaN"
        else:
            band_gap[band] = np.min(eigenvalues[band + 1]) - np.max(eigenvalues[band])
        if band_gap[band] < band_gap_threshold:
            isolated[band] = False
            isolated[band + 1] = False

    # band group (requires band gaps)
    band_group = np.zeros(num_bands, dtype=int)
    band_group_val = 0
    for band in range(num_bands):
        if band == 0:
            band_group[band] = 0
        elif band_gap[band - 1] > band_gap_threshold:
            band_group_val = band_group_val + 1
            band_group[band] = band_group_val
        else:
            band_group[band] = band_group_val

    # compute Berry fluxes
    berry_fluxes = np.zeros((num_bands, num_samples - 1, num_samples - 1))  # real
    if wilson:
        wilson_loops = np.zeros((num_bands, num_samples))  # real
    if any(geometry_columns):
        fs_metric = np.zeros((num_bands, num_samples - 1, num_samples - 1, 2, 2))  # real
        berry_fluxes_2 = np.zeros((num_bands, num_samples - 1, num_samples - 1))  # real
        TISM, DISM = np.zeros((2, num_bands, num_samples - 1, num_samples - 1))  # real
    for band, group in tqdm(enumerate(band_group), desc="Band Properties", ascii=True):
        group_size = np.count_nonzero(band_group == group)
        for idx_x in range(num_samples - 1):
            if wilson:
                if group != band_group[band - 1]:
                    wilson_loops[band, idx_x] = fbs.wilson_loop(eigenvectors, band, idx_x, group_size)
                else:
                    wilson_loops[band, idx_x] = wilson_loops[band - 1, idx_x]
            for idx_y in range(num_samples - 1):
                if group != band_group[band - 1]:
                    berry_fluxes[band, idx_x, idx_y] = fbs.berry_curv(eigenvectors, band, idx_x, idx_y, group_size)
                    # quantum geometry
                    if any(geometry_columns):
                        if group_size == 1:
                            geom_tensor = fbs.geom_tensor(eigenvectors, eigenvectors_dkx, eigenvectors_dky, bMUCvec, band, idx_x, idx_y, group_size)
                            fs_metric[band, idx_x, idx_y] = np.real(geom_tensor)
                            berry_curv = -2*np.imag(geom_tensor)
                            ###
                            berry_fluxes_2[band, idx_x, idx_y] = berry_curv[0][1]
                            TISM[band, idx_x, idx_y] = np.trace(fs_metric[band, idx_x, idx_y]) \
                                - np.abs(berry_fluxes_2[band, idx_x, idx_y])
                            DISM[band, idx_x, idx_y] = np.linalg.det(fs_metric[band, idx_x, idx_y]) \
                                - 0.25*np.abs(berry_fluxes_2[band, idx_x, idx_y])**2
                else:
                    berry_fluxes[band, idx_x, idx_y] = berry_fluxes[band-1, idx_x, idx_y]
        if wilson:
            if group != band_group[band - 1]:
                wilson_loops[band, num_samples-1] = fbs.wilson_loop(eigenvectors, band, num_samples-1, group_size)
            else:
                wilson_loops[band, num_samples-1] = wilson_loops[band - 1, num_samples-1]

    # band properties
    band_width = np.zeros(num_bands)
    std_B_norm = np.zeros(num_bands)
    chern_numbers, chern_numbers_2 = np.zeros((2, num_bands))
    av_gxx, std_gxx, av_gxy, std_gxy = np.zeros((4, num_bands))
    av_TISM, av_DISM = np.zeros((2, num_bands))
    for band_idx, band in enumerate(np.arange(num_bands)[::-1]):
        band_width[band] = np.max(eigenvalues[band]) - np.min(eigenvalues[band])
        std_B_norm[band] = np.std(berry_fluxes[band, :, :])/np.abs(np.average(berry_fluxes[band, :, :]))
        chern_numbers[band] = np.sum(berry_fluxes[band, :, :]) / (2 * np.pi)
        if any(geometry_columns):
            chern_numbers_2[band] = np.sum(berry_fluxes_2[band, :, :]) * Dkx * Dky / (2 * np.pi)
            av_gxx[band] = np.mean(fs_metric[band, :, :, 0, 0])
            std_gxx[band] = np.std(fs_metric[band, :, :, 0, 0])
            av_gxy[band] = np.mean(fs_metric[band, :, :, 0, 1])
            std_gxy[band] = np.std(fs_metric[band, :, :, 0, 1])
            av_TISM[band] = np.sum(TISM[band]) * Dkx * Dky / (2 * np.pi)
            av_DISM[band] = np.sum(DISM[band]) * Dkx * Dky / (2 * np.pi)

    # table
    headers = ["band", "group", "isolated", "width", "gap", "gap/width", "std_B/|av_B|", "C",
               "C (geom_tensor)", "av_gxx", "std_gxx", "av_gxy", "std_gxy", "<T>", "<D>"]
    bools = [show_band, show_group, show_isolated, show_width, show_gap, show_gap_width, show_std_B_norm, show_C,
             show_C_geom_tensor, show_av_gxx, show_std_gxx, show_av_gxy, show_std_gxy, show_T, show_D]
    table = PrettyTable()
    name_list = []
    table.field_names = [j for i, j in enumerate(headers) if bools[i]]
    for band in np.arange(num_bands)[::-1]:
        data = [band, band_group[band], isolated[band], band_width[band], band_gap[band], band_gap[band] / band_width[band],
                std_B_norm[band], round(chern_numbers[band]), chern_numbers_2[band],
                av_gxx[band], std_gxx[band], av_gxy[band], std_gxy[band], av_TISM[band], av_DISM[band]]
        table.add_row([j for i, j in enumerate(data) if bools[i]])
    print(table)

    if args['display'] == "3D":
        # construct figure
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.set_title(f"$n_\phi = {args['nphi'][0]}/{args['nphi'][1]}$")
        idx_x = np.linspace(0, num_samples - 1, num_samples, dtype=int)
        idx_y = np.linspace(0, num_samples - 1, num_samples, dtype=int)
        kx, ky = np.meshgrid(idx_x, idx_y)
        for i in range(num_bands):
            ax.plot_surface(kx, ky, eigenvalues[i, kx, ky], alpha=0.5)
        ax.set_xlabel('$k_1/|\mathbf{b}_1|$')
        ax.set_ylabel('$k_2/|\mathbf{b}_2|$')
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

        if wilson:
            fig2 = plt.figure()
            ax2 = fig2.add_subplot(111)
            ax2.set_title(f"$n_\phi = {args['nphi'][0]}/{args['nphi'][1]}$")
            idx_x = np.linspace(0, num_samples - 1, num_samples, dtype=int)
            for i in range(num_bands):
                ax2.scatter(idx_x, wilson_loops[i])
            ax2.set_xlabel('$k_1/|\mathbf{b}_1|$')
            ax2.set_ylabel('$\prod \\theta_\\mathrm{B}$')
            ax2.xaxis.set_major_formatter(plt.FuncFormatter(normalize))
            ax2.yaxis.set_major_formatter(plt.FuncFormatter(normalize))

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
                k = np.matmul(k, bMUCvec)
                eigvals = np.linalg.eigvals(model.hamiltonian(k))
                idx = np.argsort(eigvals)
                for band in range(num_bands):
                    eigenvalues[band, count] = np.real(eigvals[idx[band]])
                count += 1

        # construct figure
        fig = plt.figure()
        ax = plt.subplot(111)
        ax.set_title(f"$n_\phi = {args['nphi'][0]}/{args['nphi'][1]}$")
        for i in range(num_bands):
            ax.plot(eigenvalues[i])
        for i in range(1, num_paths):
            ax.axvline(i*points_per_path, color='k', linewidth=0.5, ls='--')
        ax.set_xlim([0, num_points])
        ax.set_xlabel('symmetry path')
        ax.set_ylabel('$E$')

    plt.show()
