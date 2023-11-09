# --- external imports
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from math import gcd
import matplotlib.ticker as ticker
from fractions import Fraction
from matplotlib.ticker import MaxNLocator
from copy import deepcopy
import matplotlib.colors as mcolors
import os
# --- internal imports
import functions.arguments as fa
import functions.band_structure as fb
from models.hofstadter import Hofstadter


def chern(pval, qval):

    nphi = Fraction(pval, qval)
    p = nphi.numerator
    q = nphi.denominator

    # determine r and s
    sr_list, tr_list = [], []

    for r in range(q+1):
        if q % 2 == 0 and r == q/2:
            continue
        for tr in range(-int(q/2), int(q/2)+1):
            for sr in range(-q, q+1):
                if r == q*sr + p*tr:
                    sr_list.append(sr)
                    tr_list.append(tr)
                    break
            else:
                continue  # only executed if the inner loop did NOT break
            break  # only executed if the inner loop DID break

    Chern_list = []
    if q % 2 != 0:
        numb_band_groups = q
    else:
        numb_band_groups = q-1

    for i in range(numb_band_groups):
        Chern_list.append(tr_list[i+1] - tr_list[i])

    if q % 2 == 0:
        Chern_list.insert(q//2-1, Chern_list[q//2-1])

    return Chern_list, tr_list


if __name__ == '__main__':

    # input arguments
    args = fa.parse_input_arguments("butterfly", "Plots the Hofstadter Butterfly.")
    t = args['t']
    lat = args['lattice']
    alpha = args['alpha']
    theta = args['theta']
    save = args['save']
    q = args['q']
    color = args['color']
    pal = args['palette']
    wan = args['wannier']
    period = args['periodicity']
    art = args['art']
    dpi = args['dpi']

    # construct butterfly
    nphi_list, E_list = [], []
    if color:
        chern_list, tr_list = [], []
    if wan:
        nphi_DOS_list, DOS_list, gaps_list, tr_DOS_list = [], [], [], []
    for p in tqdm(range(1, q), desc="Butterfly Construction", ascii=True):
        if args['model'] == "Hofstadter":
            model = Hofstadter(p, q, t=t, lat=lat, alpha=alpha, theta=theta, period=period)
        else:
            raise ValueError("model is not defined")
        if gcd(p, q) != 1:  # nphi must be a coprime fraction
            continue
        nphi = p / q

        ham = model.hamiltonian(np.array([0, 0]))

        M = len(ham)
        nphi_list.append([nphi] * M)

        lmbda = np.sort(np.linalg.eigvalsh(ham))
        E_list.append(lmbda)

        if wan:
            nphi_DOS_list.append([nphi] * (M-1))
            DOS_list.append([i/q for i in range(len(lmbda)-1)])
            gaps_list.append([lmbda[i+1] - lmbda[i] for i in range(len(lmbda)-1)])

        if color:
            cherns, trs = chern(p, q)

            chern_list.append(cherns)
            tr_list.append(trs)
            if wan:
                tr_DOS_list.append(trs[1:-1])

            # if model.lat == "honeycomb":
            #     cherns_double = cherns + [i for i in cherns[::-1]]
            #     chern_list.append(cherns_double)
            #     trs_double = trs + [-i for i in trs[::-1]]
            #     tr_list.append(trs_double)
            #     if wan:
            #         trs_double_wan = trs[1:] + [-i for i in trs[::-1]][1:]
            #         tr_DOS_list.append(trs_double_wan[:-1])
            # else:
            #     chern_list.append(cherns)
            #     tr_list.append(trs)
            #     if wan:
            #         tr_DOS_list.append(trs[1:-1])

    if color == "plane":

        E_list_orig = deepcopy(E_list)

        if round(M/q) == 2:
            half_len = int(np.shape(E_list)[1]/2)
            for i, val in enumerate(E_list):
                E_list[i] = val[:half_len]  # consider only lower half

        resx = np.shape(E_list)[0]
        resy = np.shape(E_list)[1]
        res = [resx, resy]

        E_vals = np.linspace(np.min(E_list), np.max(E_list), res[1])  # energy bins
        matrix = np.zeros((res[0], res[1]))

        for i, p in enumerate(range(1, res[0]+1)):  # p goes from 1 to 199
            for j, E in enumerate(E_vals):  # E goes through energy bins
                for k, El in enumerate(E_list[i]):  # for each energy bin, compare it to the E_list
                    if E <= El:  # if energy is lower than sorted E_list value
                        matrix[i][j] = tr_list[i][k]  # assign the corresponding tr of that E_list value
                        break

        if round(M/q) == 2:
            matrix = np.concatenate((matrix, -matrix[:, ::-1]), axis=1)  # double the spectrum

    # construct figure
    fig = plt.figure()
    ax = plt.subplot(111)

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
            tick_locs = np.linspace(-10, 10, 2*21 + 1)[1::2]
            cbar_tick_label = np.arange(-10, 10 + 1)
            cbar.set_ticks(tick_locs)
            cbar.set_ticklabels(cbar_tick_label)
    elif color == "plane":
        sc = ax.imshow(matrix.T, origin='lower', cmap=cmap, extent=[0, 1, np.min(E_list_orig[0]), np.max(E_list_orig[0])],
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

        nphi_DOS_list = list(np.concatenate(nphi_DOS_list).ravel())
        DOS_list = list(np.concatenate(DOS_list).ravel())
        gaps_list = list(np.concatenate(gaps_list).ravel())

        print(np.shape(nphi_DOS_list), np.shape(DOS_list), np.shape(gaps_list), np.shape(tr_DOS_list))

        if not color:
            ax2.scatter(nphi_DOS_list, DOS_list, s=[5*i for i in gaps_list], c='r', linewidths=0)
        else:
            tr_DOS_list = list(np.concatenate(tr_DOS_list).ravel())
            sc2 = ax2.scatter(nphi_DOS_list, DOS_list, s=[10*i for i in gaps_list], c=tr_DOS_list, cmap=cmap, linewidths=0, vmin=-10, vmax=10)
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
    plt.show()
