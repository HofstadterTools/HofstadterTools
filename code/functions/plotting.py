"""Functions for plotting data."""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os
from matplotlib.ticker import MaxNLocator
from copy import deepcopy
import matplotlib.colors as mcolors


def butterfly(model, args, nphi_list, E_list, E_list_orig, chern_list, matrix, nphi_DOS_list, DOS_list, gaps_list, tr_DOS_list):

    # read input arguments
    if args['input']:
        t = fa.read_t_from_file()
    else:
        t = args['t']
    lat = args['lattice']
    alpha = args['alpha']
    theta = args['theta']
    save = args['save']
    plt_lat = args["plot_lattice"]
    q = args['q']
    color = args['color']
    pal = args['palette']
    wan = args['wannier']
    period = args['periodicity']
    art = args['art']
    dpi = args['dpi']

    # construct figure
    fig = plt.figure()
    ax = fig.add_subplot(111)

    if not art:
        ax.set_title(f"$n_\phi = p/{q}$")
        transparent = False
    else:
        transparent = True
    if color:  # define color palette
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
    if color == "point":
        sc = ax.scatter(nphi_list, E_list, c=chern_list, cmap=cmap, s=1, marker='.', vmin=-10, vmax=10)
        if not art:
            cbar = plt.colorbar(sc, extend='both')
            cbar.set_label("$C$")
            tick_locs = np.linspace(-10, 10, 2 * 21 + 1)[1::2]
            cbar_tick_label = np.arange(-10, 10 + 1)
            cbar.set_ticks(tick_locs)
            cbar.set_ticklabels(cbar_tick_label)
    elif color == "plane":
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
        ax.scatter(nphi_list, E_list, s=1, marker='.')

    if not art:
        ax.set_ylabel('$E$')
        ax.set_xlabel('$n_\phi$')
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))

    if wan:
        fig2 = plt.figure()
        ax2 = fig2.add_subplot(111)
        if not art:
            ax2.set_title(f"$n_\phi = p/{q}$")
            ax2.set_ylabel('$D(E)$')
            ax2.set_xlabel('$n_\phi$')
            ax2.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
            ax2.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
        else:
            ax2.set_xlim([0, 1])

        nphi_DOS_list = list(np.concatenate(nphi_DOS_list).ravel())
        DOS_list = list(np.concatenate(DOS_list).ravel())
        gaps_list = list(np.concatenate(gaps_list).ravel())

        if not color:
            ax2.scatter(nphi_DOS_list, DOS_list, s=[5 * i for i in gaps_list], c='r', linewidths=0)
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
        if wan:
            ax2.axis('off')

    if save:
        t_str = "t_"+'_'.join([f"{i:g}" for i in t])+"_"
        brav_str = f"alpha_{alpha:g}_theta_{theta[0]:g}_{theta[1]:g}_" if lat not in ["square", "triangular"] else ""
        col_str = f"col_{color}_{pal}_" if color else ""
        per_str = f"period_{period:g}_" if period != 1 else ""
        art_str = "art_" if art else ""
        dpi_str = f"dpi_{dpi:g}" if dpi != 300 else ""
        # create file name
        file_name = f"{lat}_q_{q:g}_{t_str}{brav_str}{col_str}{per_str}{art_str}{dpi_str}.png".replace("_.", ".")
        dir = "../figs/" if os.path.isdir('../figs') else ""
        fig.savefig(dir+f"butterfly_{file_name}", bbox_inches='tight', dpi=dpi, transparent=transparent)
        if wan:
            fig2.savefig(dir + f"wannier_{file_name}", bbox_inches='tight', dpi=dpi, transparent=transparent)

    if plt_lat:
        model.plot_lattice()

    return None
