# --- external imports
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from math import gcd
import matplotlib.ticker as ticker
from fractions import Fraction
from matplotlib.ticker import MaxNLocator
from copy import deepcopy
# --- internal imports
import functions.arguments as fa
import functions.band_structure as fb
from models.hofstadter import Hofstadter

# plt.rc('text', usetex=True)
# plt.rc('text.latex', preamble=r'\usepackage{amsmath}')


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
    args = fa.parse_input_arguments("butterfly")
    t = args['t']
    lat = args['lattice']
    alpha = args['alpha']
    theta = args['theta']
    q = args['q']
    color = args['color']

    # construct butterfly
    nphi_list, E_list = [], []
    if color:
        chern_list, tr_list = [], []
    for p in tqdm(range(1, q), desc="Butterfly Construction", ascii=True):
        if args['model'] == "Hofstadter":
            model = Hofstadter(p, q, t=t, lat=lat, alpha=alpha, theta=theta)
        else:
            raise ValueError("model is not defined")
        if gcd(p, q) != 1:  # nphi must be a coprime fraction
            continue
        nphi = p / q

        ham_list, basis = model.hamiltonian(np.array([0, 0]))
        if basis == 2:
            M = 2 * q
            # M = q
        else:
            M = q
        nphi_list.append([nphi] * M)
        #lmbda = np.sort(np.linalg.eigvalsh(ham))

        if basis == 2:
            # eenergies = np.zeros(2 * len(lmbda))
            # for i in range(len(lmbda)):
            #     if lmbda[i] < -3:  # avoid taking sqrt of negative number
            #         lmbda[i] = -3
            #     eenergies[i] = +np.sqrt(3 + lmbda[i])
            #     eenergies[len(lmbda) + i] = -np.sqrt(3 + lmbda[i])
            _, eenergies = fb.polyeig(-ham_list[1], -ham_list[0], np.eye(q))
            E_list.append(np.sort(eenergies))
        else:
            lmbda = np.sort(np.linalg.eigvalsh(ham_list[0]))
            E_list.append(lmbda)
        if color:
            cherns, trs = chern(p, q)
            if model.lat == "honeycomb":
                chern_list.append(cherns + [i for i in cherns[::-1]])
                tr_list.append(trs + [-i for i in trs[::-1]])
            else:
                chern_list.append(cherns)
                tr_list.append(trs)

    if color == "plane":

        E_list_orig = deepcopy(E_list)

        if basis == 2:
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

        if basis == 2:
            matrix = np.concatenate((matrix, -matrix[:, ::-1]), axis=1)  # double the spectrum

    # construct figure
    fig = plt.figure()
    ax = plt.subplot(111)

    ax.set_title(f"$n_\phi = p/{q}$")
    if color == "point":
        cmap = plt.get_cmap('jet', 21)
        sc = ax.scatter(nphi_list, E_list, c=chern_list, cmap=cmap, s=1, marker='.', vmin=-10, vmax=10)
        cbar = plt.colorbar(sc, extend='both')
        cbar.set_label("$C$")
        tick_locs = np.linspace(-10, 10, 2*21 + 1)[1::2]
        cbar_tick_label = np.arange(-10, 10 + 1)
        cbar.set_ticks(tick_locs)
        cbar.set_ticklabels(cbar_tick_label)
    elif color == "plane":
        cmap = plt.get_cmap('jet', 21)
        sc = ax.imshow(matrix.T, origin='lower', cmap=cmap, extent=[0, 1, np.min(E_list_orig[0]), np.max(E_list_orig[0])],
                       aspect="auto", vmin=-10, vmax=10)
        cbar = plt.colorbar(sc)
        cbar.set_label("$C$")
        tick_locs = np.linspace(-10, 10, 2 * 21 + 1)[1::2]
        cbar_tick_label = np.arange(-10, 10 + 1)
        cbar.set_ticks(tick_locs)
        cbar.set_ticklabels(cbar_tick_label)
    else:
        nphi_list = list(np.concatenate(nphi_list).ravel())
        E_list = list(np.concatenate(E_list).ravel())

        ax.scatter(nphi_list, E_list, s=1, marker='.')

    ax.set_ylabel('$E$')
    ax.set_xlabel('$n_\phi$')
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('$%g$'))

    plt.show()
