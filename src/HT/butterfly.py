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
import sys
import warnings
# --- internal imports
from HT.functions import arguments as fa
from HT.functions import butterfly as fb
from HT.functions import utility as fu
from HT.functions import plotting as fp
from HT.models.hofstadter import Hofstadter


def main():
    # parse input arguments
    args = fa.parse_input_arguments("butterfly", "Plot the Hofstadter Butterfly.")
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
    # butterfly arguments
    q = args['q']
    col = args['color']
    pal = args['palette']
    wan = args['wannier']
    art = args['art']

    # initialize logger
    if log:
        sys.stdout = sys.stderr = fu.Logger("butterfly", args)

    # initialize data array
    data = {'nphi_list': [], 'E_list': [],
            'chern_list': [], 'tr_list': [],
            'nphi_DOS_list': [], 'DOS_list': [], 'gaps_list': [], 'tr_DOS_list': [],
            'E_list_orig': [], 'matrix': None}

    # construct butterfly
    for p in tqdm(range(1, q), desc="Butterfly Construction", ascii=True):

        # construct model
        if mod == "Hofstadter":
            model = Hofstadter(p, q, a0=a, t=t, lat=lat, alpha=alpha, theta=theta, period=period)
        else:
            raise ValueError("model is not defined.")

        # define flux density
        if gcd(p, q) != 1:  # nphi must be a coprime fraction
            continue
        nphi = p / q

        # diagonalize Hamiltonian
        ham = model.hamiltonian(np.array([0, 0]))
        M = len(ham)
        data['nphi_list'].append([nphi] * M)
        lmbda = np.sort(np.linalg.eigvalsh(ham))
        data['E_list'].append(lmbda)

        # Wannier diagram data lists
        if wan:
            data['nphi_DOS_list'].append([nphi] * (M - 1))
            data['DOS_list'].append([i / M for i in range(M - 1)])
            data['gaps_list'].append([lmbda[i + 1] - lmbda[i] for i in range(M - 1)])

        # color data lists
        if col:
            cherns, trs = fb.chern(p, q)
            if round(M / q) == 2 and len(t) == 1:  # NN honeycomb
                cherns_double = cherns + [i for i in cherns[::-1]]
                data['chern_list'].append(cherns_double)
                trs_double = trs + [-i for i in trs[::-1]]
                data['tr_list'].append(trs_double)
                if wan:
                    trs_double_wan = trs[1:] + [-i for i in trs[::-1]][1:]
                    data['tr_DOS_list'].append(trs_double_wan[:-1])
            elif round(M / q) == 1:
                data['chern_list'].append(cherns)
                data['tr_list'].append(trs)
                if wan:
                    data['tr_DOS_list'].append(trs[1:-1])
            else:
                warnings.warn(
                    "Color and wannier are only implemented for square/triangular/bravais/[honeycomb+1NN] models. Continuing without color and wannier...")
                args['color'] = False
                col = args['color']
                args['wannier'] = False
                wan = args['wannier']

    # color plane data matrix
    if col == "plane":
        data['E_list_orig'] = deepcopy(data['E_list'])

        if round(M / q) == 2 and len(t) == 1:  # NN honeycomb
            half_len = int(np.shape(data['E_list'])[1] / 2)
            for i, val in enumerate(data['E_list']):
                data['E_list'][i] = val[:half_len]  # consider only lower half

        resx = np.shape(data['E_list'])[0]
        resy = np.shape(data['E_list'])[1]
        res = [resx, resy]

        E_vals = np.linspace(np.min(data['E_list']), np.max(data['E_list']), res[1])  # energy bins
        data['matrix'] = np.zeros((res[0], res[1]))

        for i, p in enumerate(range(1, res[0] + 1)):  # p goes from 1 to 199
            for j, E in enumerate(E_vals):  # E goes through energy bins
                for k, El in enumerate(data['E_list'][i]):  # for each energy bin, compare it to the E_list
                    if E <= El:  # if energy is lower than sorted E_list value
                        data['matrix'][i][j] = data['tr_list'][i][k]  # assign the corresponding tr of that E_list value
                        break

        if round(M / q) == 2 and len(t) == 1:  # NN honeycomb
            data['matrix'] = np.concatenate((data['matrix'], -data['matrix'][:, ::-1]), axis=1)  # double the spectrum

    # save data
    if save:
        fu.save_data("butterfly", model, args, data)

    # construct figure(s)
    fp.butterfly(model, args, data)
    plt.show()


if __name__ == '__main__':

    main()
