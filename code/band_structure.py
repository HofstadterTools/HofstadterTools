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

# plt.rc('text', usetex=True)
# plt.rc('text.latex', preamble=r'\usepackage{amsmath}')


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
        if band_gap[band] < 0.1:
            isolated[band] = False
            isolated[band + 1] = False

    # compute Berry fluxes
    berry_fluxes = np.zeros((num_bands, num_samples - 1, num_samples - 1))  # real
    fs_metric = np.zeros((num_bands, num_samples - 1, num_samples - 1, 2, 2))  # real
    berry_fluxes_2 = np.zeros((num_bands, num_samples - 1, num_samples - 1))  # real
    tr_g = np.zeros((num_bands, num_samples - 1, num_samples - 1))  # real
    abs_B = np.zeros((num_bands, num_samples - 1, num_samples - 1))  # real
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
                    # quantum geometry
                    geom_tensor = fbs.geom_tensor(eigenvectors[:, band, idx_x, idx_y],
                                                  eigenvectors[:, band, idx_x + 1, idx_y],
                                                  eigenvectors[:, band, idx_x, idx_y + 1])
                    fs_metric[band, idx_x, idx_y] = np.real(geom_tensor)
                    berry_curv = -2*np.imag(geom_tensor)
                    ###
                    berry_fluxes_2[band, idx_x, idx_y] = berry_curv[0][1]
                    tr_g[band, idx_x, idx_y] = np.trace(fs_metric[band, idx_x, idx_y])
                    abs_B[band, idx_x, idx_y] = np.abs(berry_fluxes_2[band, idx_x, idx_y])
                    TISM[band, idx_x, idx_y] = np.trace(fs_metric[band, idx_x, idx_y]) \
                                               - np.abs(berry_fluxes_2[band, idx_x, idx_y])
                    DISM[band, idx_x, idx_y] = np.linalg.det(fs_metric[band, idx_x, idx_y]) \
                                               - 0.25*np.abs(berry_fluxes_2[band, idx_x, idx_y])**2
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
    chern_numbers_2 = np.zeros(num_bands)
    std_B_norm = np.zeros(num_bands)
    av_gxx = np.zeros(num_bands)
    std_gxx = np.zeros(num_bands)
    av_gxy = np.zeros(num_bands)
    std_gxy = np.zeros(num_bands)
    av_gyy = np.zeros(num_bands)
    std_gyy = np.zeros(num_bands)
    av_tr_g = np.zeros(num_bands)
    av_abs_B = np.zeros(num_bands)
    band_width = np.zeros(num_bands)
    av_TISM = np.zeros(num_bands)
    av_DISM = np.zeros(num_bands)
    for band_idx, band in enumerate(np.arange(num_bands)[::-1]):
        chern_numbers[band] = np.sum(berry_fluxes[band, :, :]) / (2 * np.pi)
        chern_numbers_2[band] = np.sum(berry_fluxes_2[band, :, :]) / (2 * np.pi)
        std_B_norm[band] = np.std(berry_fluxes[band, :, :])/np.abs(np.average(berry_fluxes[band, :, :]))
        av_gxx[band] = np.mean(fs_metric[band, :, :, 0, 0])
        std_gxx[band] = np.std(fs_metric[band, :, :, 0, 0])
        av_gxy[band] = np.mean(fs_metric[band, :, :, 0, 1])
        std_gxy[band] = np.std(fs_metric[band, :, :, 0, 1])
        av_gyy[band] = np.mean(fs_metric[band, :, :, 1, 1])
        std_gyy[band] = np.std(fs_metric[band, :, :, 1, 1])
        av_tr_g[band] = np.mean(tr_g[band])
        av_abs_B[band] = np.mean(abs_B[band])
        band_width[band] = np.max(eigenvalues[band]) - np.min(eigenvalues[band])
        av_TISM[band] = np.sum(TISM[band])
        av_DISM[band] = np.sum(DISM[band])

    # table
    show_band = True
    show_isolated = True
    show_C = False
    show_C_geom_tensor = False
    show_std_B_norm = True
    show_av_gxx = True
    show_std_gxx = False
    show_av_gxy = False
    show_std_gxy = False
    show_av_gyy = True
    show_std_gyy = False
    show_av_tr_g = True
    show_av_abs_B = True
    show_width = False
    show_gap = False
    show_gap_width = False
    show_T = True
    show_D = True

    headers = ["band", "isolated", "C", "C (geom_tensor)", "std_B/|av_B|", "av_gxx", "std_gxx", "av_gxy", "std_gxy",
               "av_gyy", "std_gyy", "av_tr_g", "av_abs_B", "width", "gap", "gap/width", "<T>", "<D>"]
    bools = [show_band, show_isolated, show_C, show_C_geom_tensor, show_std_B_norm, show_av_gxx, show_std_gxx,
             show_av_gxy, show_std_gxy, show_av_gyy, show_std_gyy, show_av_tr_g, show_av_abs_B, show_width, show_gap,
             show_gap_width, show_T, show_D]

    table = PrettyTable()
    name_list = []
    table.field_names = [j for i, j in enumerate(headers) if bools[i]]
    for band in np.arange(num_bands)[::-1]:
        data = [band, isolated[band], round(chern_numbers[band]), chern_numbers_2[band], std_B_norm[band],
                av_gxx[band], std_gxx[band], av_gxy[band], std_gxy[band], av_gyy[band], std_gyy[band],
                av_tr_g[band], av_abs_B[band], band_width[band], band_gap[band], band_gap[band] / band_width[band],
                av_TISM[band], av_DISM[band]]
        table.add_row([j for i, j in enumerate(data) if bools[i]])
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
        ax.set_xlabel('$k_1/|\mathbf{b}_1|$')
        ax.set_ylabel('$k_2/|\mathbf{b}_2|$')


        def normalize(value, tick_number):

            if value == 0:
                return "$0$"
            elif value == num_samples - 1:
                return "$1$"
            else:
                return f"${value / (num_samples - 1):.1g}$"

        ax.xaxis.set_major_formatter(plt.FuncFormatter(normalize))
        ax.yaxis.set_major_formatter(plt.FuncFormatter(normalize))

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
