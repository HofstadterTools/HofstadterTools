# --- external imports
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from prettytable import PrettyTable
from tqdm import tqdm
from math import gcd
import sys
import warnings
# --- internal imports
from HT.functions import band_structure as fbs
from HT.functions import arguments as fa
from HT.functions import utility as fu
from HT.functions import plotting as fp
from HT.functions import models as fm
from HT.models.hofstadter import Hofstadter


def main():
    # parse input arguments
    args = fa.parse_input_arguments("band_structure", "Plot the Hofstadter Band Structure.")
    # general arguments
    mod = args['model']
    a = args['a']
    t = fu.read_t_from_file() if args['input'] else args['t']
    lat = args['lattice']
    alpha = args['alpha']
    theta = args['theta']
    save = args['save']
    log = args['log']
    plt_lat = args["plot_lattice"]
    period = args['periodicity']
    dpi = args['dpi']
    ps = args['point_size']
    # band_structure arguments
    samp = args['samp']
    wil = args['wilson']
    disp = args['display']
    nphi = args['nphi']
    if gcd(nphi[0], nphi[1]) != 1:  # nphi must be a coprime fraction
        raise ValueError("nphi must be a coprime fraction.")
    bgt = args['bgt']
    load = args['load']
    topo = args['topology']
    geom = args['geometry']
    cols = args['columns']

    # table column selection
    if cols == ["band", "group", "isolated", "width", "gap", "gap_width"]:  # default
        show_band = show_group = show_isolated = show_width = show_gap = show_gap_width = True
        show_std_B = show_C = topo
        show_std_g = show_T = show_D = geom
        show_av_gxx = show_std_gxx = show_av_gxy = show_std_gxy = False
    else:
        show_band = True if "band" in cols else False
        show_group = True if "group" in cols else False
        show_isolated = True if "isolated" in cols else False
        show_width = True if "width" in cols else False
        show_gap = True if "gap" in cols else False
        show_gap_width = True if "gap_width" in cols else False
        show_std_B = True if "std_B" in cols else False
        show_C = True if "C" in cols else False
        show_std_g = True if "std_g" in cols else False
        show_av_gxx = True if "av_gxx" in cols else False
        show_std_gxx = True if "std_gxx" in cols else False
        show_av_gxy = True if "av_gxy" in cols else False
        show_std_gxy = True if "std_gxy" in cols else False
        show_T = True if "T" in cols else False
        show_D = True if "D" in cols else False
    topo_cols = [show_C, show_std_B]
    geom_cols = [show_std_g, show_av_gxx, show_std_gxx, show_av_gxy, show_std_gxy, show_T, show_D]

    # initialize logger
    if log:
        sys.stdout = sys.stderr = fu.Logger("band_structure", args)

    if not load:
        # initialize data array
        data = {'eigenvalues': None, 'eigenvalues_2D': None, 'eigenvectors': None,
                'eigenvectors_dkx': None, 'eigenvectors_dky': None, 'Dkx': None, 'Dky': None,
                'wilson_loops': None}

        # construct model
        if mod == "Hofstadter":
            model = Hofstadter(nphi[0], nphi[1], a0=a, t=t, lat=lat, alpha=alpha, theta=theta, period=period)
        else:
            raise ValueError("model is not defined.")

        # define unit cell
        _, avec, abasisvec, bMUCvec, sym_points = model.unit_cell()
        _, bases = fm.nearest_neighbor_finder(avec, abasisvec, t, 0, 0, 0)
        num_bands = nphi[1] * len(bases)

        if disp == "3D" or disp == "both" or wil or any(topo_cols) or any(geom_cols):
            # construct bands
            data['eigenvalues'] = np.zeros((num_bands, samp, samp))  # real
            data['eigenvectors'] = np.zeros((num_bands, num_bands, samp, samp), dtype=np.complex128)  # complex
            if any(geom_cols):
                data['eigenvectors_dkx'] = np.zeros((num_bands, num_bands, samp, samp), dtype=np.complex128)  # complex
                data['eigenvectors_dky'] = np.zeros((num_bands, num_bands, samp, samp), dtype=np.complex128)  # complex
                data['Dkx'] = np.dot(np.array([1 / (samp - 1), 0]), bMUCvec[0])
                data['Dky'] = np.dot(np.array([0, 1 / (samp - 1)]), bMUCvec[1])
            for band in tqdm(range(num_bands), desc="Band Construction", ascii=True):
                for idx_x in range(samp):
                    frac_kx = idx_x / (samp - 1)
                    if any(geom_cols):
                        frac_kx_dkx = (frac_kx + 1 / (1000 * (samp - 1))) % 1
                    for idx_y in range(samp):
                        frac_ky = idx_y / (samp - 1)
                        k = np.matmul(np.array([frac_kx, frac_ky]), bMUCvec)
                        ham = model.hamiltonian(k)
                        eigvals, eigvecs = np.linalg.eigh(ham)
                        idx = np.argsort(eigvals)
                        data['eigenvalues'][band, idx_x, idx_y] = eigvals[idx[band]]
                        data['eigenvectors'][:, band, idx_x, idx_y] = eigvecs[:, idx[band]]
                        if any(geom_cols):
                            frac_ky_dky = (frac_ky + 1 / (1000 * (samp - 1))) % 1
                            k_dkx = np.matmul(np.array([frac_kx_dkx, frac_ky]), bMUCvec)
                            k_dky = np.matmul(np.array([frac_kx, frac_ky_dky]), bMUCvec)
                            ham_dkx = model.hamiltonian(k_dkx)
                            eigvals_dkx, eigvecs_dkx = np.linalg.eigh(ham_dkx)
                            ham_dky = model.hamiltonian(k_dky)
                            eigvals_dky, eigvecs_dky = np.linalg.eigh(ham_dky)
                            idx_dkx = np.argsort(eigvals_dkx)
                            idx_dky = np.argsort(eigvals_dky)
                            data['eigenvectors_dkx'][:, band, idx_x, idx_y] = eigvecs_dkx[:, idx_dkx[band]]
                            data['eigenvectors_dky'][:, band, idx_x, idx_y] = eigvecs_dky[:, idx_dky[band]]

        if disp == "2D" or disp == "both":
            # construct bands
            num_paths = len(sym_points)
            points_per_path = int(samp / num_paths)
            num_points = num_paths * points_per_path
            data['eigenvalues_2D'] = np.zeros((num_bands, num_points))  # real
            count = 0
            for i in range(num_paths):
                for j in range(points_per_path):
                    k = sym_points[i][1] + (sym_points[(i + 1) % num_paths][1] - sym_points[i][1]) * float(j) / float(
                        points_per_path - 1)
                    k = np.matmul(k, bMUCvec)
                    ham = model.hamiltonian(k)
                    eigvals = np.linalg.eigvalsh(ham)
                    idx = np.argsort(eigvals)
                    for band in range(num_bands):
                        data['eigenvalues_2D'][band, count] = np.real(eigvals[idx[band]])
                    count += 1
    else:  # load from file
        model, args_load, data = fu.load_data("band_structure", load)

        # fix certain arguments
        # general arguments
        mod = args_load['model']
        t = fu.read_t_from_file() if args_load['input'] else args_load['t']
        lat = args_load['lattice']
        alpha = args_load['alpha']
        theta = args_load['theta']
        period = args_load['periodicity']
        # band_structure arguments
        samp = args_load['samp']
        nphi = args_load['nphi']

        # define unit cell
        num_bands, _, _, bMUCvec, sym_points = model.unit_cell()

    # band gap and isolated
    if disp == "3D" or disp == "both" or wil or any(topo_cols) or any(geom_cols):
        eigenvals = data['eigenvalues']
    else:
        eigenvals = data['eigenvalues_2D']
    band_gap = np.zeros(num_bands)
    isolated = np.full(num_bands, True)
    for band_idx, band in enumerate(np.arange(num_bands)[::-1]):
        if band_idx == 0:
            band_gap[band] = "NaN"
        else:
            band_gap[band] = np.min(eigenvals[band + 1]) - np.max(eigenvals[band])
        if band_gap[band] < bgt:
            isolated[band] = False
            isolated[band + 1] = False

    # band group (requires band gaps)
    band_group = np.zeros(num_bands, dtype=int)
    band_group_val = 0
    for band in range(num_bands):
        if band == 0:
            band_group[band] = 0
        elif band_gap[band - 1] > bgt:
            band_group_val = band_group_val + 1
            band_group[band] = band_group_val
        else:
            band_group[band] = band_group_val

    # compute Berry fluxes
    if wil or any(topo_cols) or any(geom_cols):
        berry_fluxes = np.zeros((num_bands, samp - 1, samp - 1))  # real
        if wil:
            data['wilson_loops'] = np.zeros((num_bands, samp))  # real
        if any(geom_cols):
            fs_metric = np.zeros((num_bands, samp - 1, samp - 1, 2, 2))  # real
            berry_fluxes_2 = np.zeros((num_bands, samp - 1, samp - 1))  # real
            TISM, DISM = np.zeros((2, num_bands, samp - 1, samp - 1))  # real
        for band, group in tqdm(enumerate(band_group), desc="Band Properties", ascii=True):
            group_size = np.count_nonzero(band_group == group)
            for idx_y in range(samp - 1):
                if wil:
                    if group != band_group[band - 1]:
                        data['wilson_loops'][band, idx_y] = fbs.wilson_loop(data['eigenvectors'], band, idx_y,
                                                                            group_size)
                    else:
                        data['wilson_loops'][band, idx_y] = data['wilson_loops'][band - 1, idx_y]
                for idx_x in range(samp - 1):
                    if group != band_group[band - 1]:
                        berry_fluxes[band, idx_x, idx_y] = fbs.berry_curv(data['eigenvectors'], band, idx_x, idx_y,
                                                                          group_size)
                        # quantum geometry
                        if any(geom_cols):
                            geom_tensor = fbs.geom_tensor(data['eigenvectors'], data['eigenvectors_dkx'],
                                                          data['eigenvectors_dky'], bMUCvec, band, idx_x, idx_y,
                                                          group_size)
                            fs_metric[band, idx_x, idx_y] = np.real(geom_tensor)
                            berry_curv = -2 * np.imag(geom_tensor)
                            ###
                            berry_fluxes_2[band, idx_x, idx_y] = berry_curv[0][1]
                            TISM[band, idx_x, idx_y] = np.trace(fs_metric[band, idx_x, idx_y]) \
                                                       - np.abs(berry_fluxes_2[band, idx_x, idx_y])
                            DISM[band, idx_x, idx_y] = np.linalg.det(fs_metric[band, idx_x, idx_y]) \
                                                       - 0.25 * np.abs(berry_fluxes_2[band, idx_x, idx_y]) ** 2
                    else:
                        berry_fluxes[band, idx_x, idx_y] = berry_fluxes[band - 1, idx_x, idx_y]
                        if any(geom_cols):
                            fs_metric[band, idx_x, idx_y] = fs_metric[band - 1, idx_x, idx_y]
                            berry_fluxes_2[band, idx_x, idx_y] = berry_fluxes_2[band - 1, idx_x, idx_y]
                            TISM[band, idx_x, idx_y] = TISM[band - 1, idx_x, idx_y]
                            DISM[band, idx_x, idx_y] = DISM[band - 1, idx_x, idx_y]
            if wil:
                if group != band_group[band - 1]:
                    data['wilson_loops'][band, samp - 1] = fbs.wilson_loop(data['eigenvectors'], band, samp - 1,
                                                                           group_size)
                else:
                    data['wilson_loops'][band, samp - 1] = data['wilson_loops'][band - 1, samp - 1]

    # band properties
    band_width = np.zeros(num_bands)
    std_B_norm = np.zeros(num_bands)
    chern_numbers = np.zeros(num_bands)
    std_g_norm, av_gxx, std_gxx, av_gxy, std_gxy = np.zeros((5, num_bands))
    av_TISM, av_DISM, av_TISM_proj = np.zeros((3, num_bands))
    for band_idx, band in enumerate(np.arange(num_bands)[::-1]):
        band_width[band] = np.max(eigenvals[band]) - np.min(eigenvals[band])
        if any(topo_cols):
            std_B_norm[band] = np.std(berry_fluxes[band, :, :]) / np.abs(np.average(berry_fluxes[band, :, :]))
            chern_numbers[band] = np.sum(berry_fluxes[band, :, :]) / (2 * np.pi)
        if any(geom_cols):
            g_var_sum = np.var(fs_metric[band, :, :][0, 0]) + np.var(fs_metric[band, :, :][0, 1])
            std_g_norm[band] = np.sqrt(g_var_sum)
            av_gxx[band] = np.mean(fs_metric[band, :, :, 0, 0])
            std_gxx[band] = np.std(fs_metric[band, :, :, 0, 0])
            av_gxy[band] = np.mean(fs_metric[band, :, :, 0, 1])
            std_gxy[band] = np.std(fs_metric[band, :, :, 0, 1])
            av_TISM[band] = np.sum(TISM[band]) * data['Dkx'] * data['Dky'] / (2 * np.pi)
            av_DISM[band] = np.sum(DISM[band]) * data['Dkx'] * data['Dky'] / (2 * np.pi)

    # table
    headers = ["band", "group", "isolated", "width", "gap", "gap/width", "std_B", "C",
               "std_g", "av_gxx", "std_gxx", "av_gxy", "std_gxy", "<T>", "<D>"]
    bools = [show_band, show_group, show_isolated, show_width, show_gap, show_gap_width, show_std_B, show_C,
             show_std_g, show_av_gxx, show_std_gxx, show_av_gxy, show_std_gxy, show_T, show_D]
    table = PrettyTable()
    table.field_names = [j for i, j in enumerate(headers) if bools[i]]
    for band in np.arange(num_bands)[::-1]:
        table_data = [band, band_group[band], isolated[band], band_width[band], band_gap[band],
                      band_gap[band] / band_width[band],
                      std_B_norm[band], round(chern_numbers[band]), std_g_norm[band],
                      av_gxx[band], std_gxx[band], av_gxy[band], std_gxy[band], av_TISM[band], av_DISM[band]]
        table.add_row([j for i, j in enumerate(table_data) if bools[i]])
    print(table)

    # save data
    if save:
        fu.save_data("band_structure", model, args, data)

    # construct figures
    fp.band_structure(model, args, data)
    plt.show()


if __name__ == '__main__':

    main()
