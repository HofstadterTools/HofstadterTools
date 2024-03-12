"""Functions for plotting data."""

# --- external imports
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os
from matplotlib.ticker import MaxNLocator
from copy import deepcopy
import matplotlib.colors as mcolors
from mpl_toolkits.mplot3d import axes3d
from matplotlib import rcParams
# --- internal imports
from HT.functions import utility as fu
from HT.functions import models as fm


def band_structure(model, args, data, plotting=False):
    """Plot the Hofstadter band structure.

    Parameters
    ----------
    model: Hofstadter.hamiltonian
        The Hofstadter Hamiltonian class attribute.
    args: dict
        The arguments passed to the program.
    data: ndarray
        The data array.
    plotting: bool
        The plotting script flag.
    """

    # read input arguments
    # general arguments
    mod = args['model']
    t = args['t']
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
    bgt = args['bgt']

    # read data entries
    eigenvalues = data['eigenvalues']
    eigenvalues_2D = data['eigenvalues_2D']
    wilson_loops = data['wilson_loops']

    # define unit cell
    _, avec, abasisvec, bMUCvec, sym_points = model.unit_cell()
    _, bases = fm.nearest_neighbor_finder(avec, abasisvec, t, 0, 0, 0)
    num_bands = nphi[1] * len(bases)

    if disp in ["3D", "both"]:
        # construct figure
        fig = plt.figure()
        fig.canvas.manager.set_window_title('Band Structure (3D)')
        ax = fig.add_subplot(111, projection='3d')
        ax.set_title(f"$n_\\phi = {nphi[0]}/{nphi[1]}$")
        idx_x = np.linspace(0, samp - 1, samp, dtype=int)
        idx_y = np.linspace(0, samp - 1, samp, dtype=int)
        kx, ky = np.meshgrid(idx_x, idx_y)
        for i in range(num_bands):
            ax.plot_surface(kx, ky, eigenvalues[i, kx, ky], alpha=0.5)
        ax.set_xlabel('$k_1/|\\mathbf{b}_1|$')
        ax.set_ylabel('$k_2/|\\mathbf{b}_2|$')
        ax.set_zlabel('$E$')

        def normalize(value, tick_number):
            if value == 0:
                return "$0$"
            elif value == samp - 1:
                return "$1$"
            else:
                return f"${value / (samp - 1):.1g}$"

        ax.xaxis.set_major_formatter(plt.FuncFormatter(normalize))
        ax.yaxis.set_major_formatter(plt.FuncFormatter(normalize))

        if wil:
            fig2 = plt.figure()
            fig2.canvas.manager.set_window_title('Wilson Loops')
            ax2 = fig2.add_subplot(111)
            ax2.set_title(f"$n_\\phi = {args['nphi'][0]}/{args['nphi'][1]}$")
            idx_y = np.linspace(0, samp - 1, samp, dtype=int)
            for i in range(num_bands):
                ax2.scatter(idx_y, wilson_loops[i], s=ps*rcParams['lines.markersize']**2)
            ax2.set_xlabel('$k_2/|\\mathbf{b}_2|$')
            ax2.set_ylabel('$\\prod \\theta_\\mathrm{B}$')
            ax2.xaxis.set_major_formatter(plt.FuncFormatter(normalize))
            ax2.yaxis.set_major_formatter(plt.FuncFormatter(normalize))

    if disp in ["2D", "both"]:
        # define path mesh size
        num_paths = len(sym_points)
        points_per_path = int(samp / num_paths)
        num_points = num_paths * points_per_path
        # construct figure
        fig3 = plt.figure()
        fig3.canvas.manager.set_window_title('Band Structure (2D)')
        ax3 = fig3.add_subplot(111)
        ax3.set_title(f"$n_\\phi = {nphi[0]}/{nphi[1]}$")
        for i in range(num_bands):
            ax3.plot(eigenvalues_2D[i])
        locations, labels = [0], [sym_points[0][0]]
        for i in range(1, num_paths):
            ax3.axvline(i*points_per_path, color='k', linewidth=0.5, ls='--')
            locations.append(i*points_per_path)
            labels.append(sym_points[i][0])
        locations.append(num_paths*points_per_path)
        labels.append(sym_points[0][0])
        ax3.set_xticks(locations, labels=labels)
        ax3.set_xlim([0, num_points])
        ax3.set_xlabel('reference path')
        ax3.set_ylabel('$E$')

    if save:
        rel_path = "../../.." if plotting else "../.."
        dir = f"{rel_path}/figs/band_structure/" if os.path.isdir(f'{rel_path}/figs/band_structure/') else ""
        if disp in ["2D", "both"]:
            filename_2D = fu.create_filename("band_structure", args, aux_text="2D")
            fig3.savefig(dir+filename_2D+".png", bbox_inches='tight', dpi=dpi)
        if disp in ["3D", "both"]:
            filename_3D = fu.create_filename("band_structure", args, aux_text="3D")
            fig.savefig(dir+filename_3D+".png", bbox_inches='tight', dpi=dpi)
        if wil:
            fig2.savefig(dir+filename_3D.replace("band_structure_3D", "wilson")+".png",
                         bbox_inches='tight', dpi=dpi)

    if plt_lat:
        model.plot_lattice()

    return None


def butterfly(model, args, data, plotting=False):
    """Plot the Hofstadter butterfly.

    Parameters
    ----------
    model: Hofstadter.hamiltonian
        The Hofstadter Hamiltonian class attribute.
    args: dict
        The arguments passed to the program.
    data: ndarray
        The data array.
    plotting: bool
        The plotting script flag.
    """

    # read input arguments
    # general arguments
    mod = args['model']
    t = args['t']
    lat = args['lattice']
    alpha = args['alpha']
    theta = args['theta']
    save = args['save']
    log = args['log']
    plt_lat = args["plot_lattice"]
    dpi = args['dpi']
    ps = args['point_size']
    # butterfly arguments
    q = args['q']
    col = args['color']
    pal = args['palette']
    wan = args['wannier']
    period = args['periodicity']
    art = args['art']

    # read data entries
    nphi_list = data['nphi_list']
    E_list = data['E_list']
    E_list_orig = data['E_list_orig']
    chern_list = data['chern_list']
    matrix = data['matrix']
    nphi_DOS_list = data['nphi_DOS_list']
    DOS_list = data['DOS_list']
    gaps_list = data['gaps_list']
    tr_DOS_list = data['tr_DOS_list']

    # construct figure
    fig = plt.figure()
    fig.canvas.manager.set_window_title('Butterfly Spectrum')
    ax = fig.add_subplot(111)

    if not art:
        ax.set_title(f"$n_\\phi = p/{q}$")
        transparent = False
    else:
        transparent = True
    if col:  # define color palette
        if pal == "jet":
            cmap = plt.get_cmap('jet', 21)
        elif pal == "red-blue":
            colors1 = plt.cm.Blues(np.linspace(0, 1, 10))
            colors2 = plt.cm.seismic([0.5])
            if art:
                colors2[:, -1] = 0  # set transparent region background
            colors3 = plt.cm.Reds_r(np.linspace(0, 1, 10))
            colors = np.vstack((colors1, colors2, colors3))
            cmap = mcolors.LinearSegmentedColormap.from_list('red-blue', colors, 21)
        else:  # avron
            colors1 = plt.cm.gist_rainbow(np.linspace(0.75, 1, 10)[::-1])
            colors2 = plt.cm.seismic([0.5])
            if art:
                colors2[:, -1] = 0  # set transparent region background
            colors3 = plt.cm.gist_rainbow(np.linspace(0., 0.5, 10))
            colors = np.vstack((colors1, colors2, colors3))
            cmap = mcolors.LinearSegmentedColormap.from_list('avron', colors, 21)
    if col == "point":
        sc = ax.scatter(nphi_list, E_list, c=chern_list, cmap=cmap, s=ps*7*(199/q), marker='.', vmin=-10, vmax=10, linewidths=0)
        if not art:
            cbar = plt.colorbar(sc, extend='both')
            cbar.set_label("$C$")
            tick_locs = np.linspace(-10, 10, 2 * 21 + 1)[1::2]
            cbar_tick_label = np.arange(-10, 10 + 1)
            cbar.set_ticks(tick_locs)
            cbar.set_ticklabels(cbar_tick_label)
    elif col == "plane":
        sc = ax.imshow(matrix.T, origin='lower', cmap=cmap,
                       extent=[0, 1, np.min(E_list_orig[0]), np.max(E_list_orig[0])],
                       aspect="auto", vmin=-10, vmax=10)
        if not art:
            cbar = plt.colorbar(sc)
            cbar.set_label("$t$")
            tick_locs = np.linspace(-10, 10, 2 * 21 + 1)[1::2]
            cbar_tick_label = np.arange(-10, 10 + 1)
            cbar.set_ticks(tick_locs)
            cbar.set_ticklabels(cbar_tick_label)
    else:
        nphi_list = list(np.concatenate(nphi_list).ravel())
        E_list = list(np.concatenate(E_list).ravel())
        ax.scatter(nphi_list, E_list, s=ps*7*(199/q), marker='.', linewidths=0)

    if not art:
        ax.set_ylabel('$E$')
        ax.set_xlabel('$n_\\phi$')
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))

    if wan:
        fig2 = plt.figure()
        fig2.canvas.manager.set_window_title('Wannier Diagram')
        ax2 = fig2.add_subplot(111)
        if not art:
            ax2.set_title(f"$n_\\phi = p/{q}$")
            ax2.set_ylabel('$N(E)$')
            ax2.set_xlabel('$n_\\phi$')
            ax2.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
            ax2.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
        else:
            ax2.set_xlim([0, 1])

        nphi_DOS_list = list(np.concatenate(nphi_DOS_list).ravel())
        DOS_list = list(np.concatenate(DOS_list).ravel())
        gaps_list = list(np.concatenate(gaps_list).ravel())

        if not col:
            ax2.scatter(nphi_DOS_list, DOS_list, s=[5*i for i in gaps_list], c='r', linewidths=0)
        else:
            tr_DOS_list = list(np.concatenate(tr_DOS_list).ravel())
            sc2 = ax2.scatter(nphi_DOS_list, DOS_list, s=[10 * i for i in gaps_list], c=tr_DOS_list, cmap=cmap,
                              linewidths=0, vmin=-10, vmax=10)
            if not art:
                cbar2 = plt.colorbar(sc2, extend='both')
                cbar2.set_label("$t$")
                tick_locs = np.linspace(-10, 10, 2 * 21 + 1)[1::2]
                cbar_tick_label = np.arange(-10, 10 + 1)
                cbar2.set_ticks(tick_locs)
                cbar2.set_ticklabels(cbar_tick_label)

    if art:
        ax.axis('off')
        # ax.set_box_aspect(1)
        if wan:
            ax2.axis('off')

    if save:
        filename = fu.create_filename("butterfly", args)
        rel_path = "../../../" if plotting else "../.."
        dir = f"{rel_path}/figs/butterfly/" if os.path.isdir(f'{rel_path}/figs/butterfly/') else ""
        fig.savefig(dir+filename+".png", bbox_inches='tight', dpi=dpi, transparent=transparent)
        if wan:
            fig2.savefig(dir+filename.replace("butterfly", "wannier")+".png",
                         bbox_inches='tight', dpi=dpi, transparent=transparent)

    if plt_lat:
        model.plot_lattice()

    return None
