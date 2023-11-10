"""Functions for input/output and miscellaneous use."""

import os
import csv
import numpy as np


def read_t_from_file():

    dir = "configuration/" if os.path.isdir('configuration') else ""
    with open(dir+"hopping_input.txt", 'r') as csvfile:
        data = csv.reader(csvfile, delimiter='\t')
        NN, t = [], []
        for i, row in enumerate(data):
            NN.append(int(row[0]))
            t.append(float(row[1]))

    t_list = np.zeros(max(NN))
    for i, val in enumerate(NN):
        t_list[val-1] = t[i]

    return t_list


def butterfly_filename(args):

    # read input arguments
    if args['input']:
        t = fa.read_t_from_file()
    else:
        t = args['t']
    lat = args['lattice']
    alpha = args['alpha']
    theta = args['theta']
    q = args['q']
    color = args['color']
    pal = args['palette']
    period = args['periodicity']
    art = args['art']
    dpi = args['dpi']

    t_str = "t_" + '_'.join([f"{i:g}" for i in t]) + "_"
    brav_str = f"alpha_{alpha:g}_theta_{theta[0]:g}_{theta[1]:g}_" if lat not in ["square", "triangular"] else ""
    col_str = f"col_{color}_{pal}_" if color else ""
    per_str = f"period_{period:g}_" if period != 1 else ""
    art_str = "art_" if art else ""
    dpi_str = f"dpi_{dpi:g}" if dpi != 300 else ""

    filename = f"butterfly_{lat}_q_{q:g}_{t_str}{brav_str}{col_str}{per_str}{art_str}{dpi_str}"[:-1]  # remove last underscore

    return filename


def save_data(model, args, nphi_list, E_list, E_list_orig, chern_list, matrix, nphi_DOS_list, DOS_list, gaps_list, tr_DOS_list):

    filename = butterfly_filename(args)
    save_list = np.array([model, args, nphi_list, E_list,
                          E_list_orig, chern_list, matrix, nphi_DOS_list, DOS_list, gaps_list, tr_DOS_list], dtype=object)
    dir = "../data/" if os.path.isdir('../data') else ""
    np.save(dir+filename, save_list)

    return None
