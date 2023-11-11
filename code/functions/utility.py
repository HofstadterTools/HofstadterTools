"""Functions for input/output and miscellaneous use."""

import os
import csv
import numpy as np
import sys


def read_t_from_file():
    """Reads the hopping amplitudes from file.

    Returns
    -------
    t_list: list
        The list of hopping amplitudes in order of ascending NN.
    """

    directory = "configuration/" if os.path.isdir('configuration') else ""
    with open(directory+"hopping_input.txt", 'r') as csvfile:
        data = csv.reader(csvfile, delimiter='\t')
        NN, t = [], []
        for i, row in enumerate(data):
            NN.append(int(row[0]))
            t.append(float(row[1]))

    t_list = np.zeros(max(NN))
    for i, val in enumerate(NN):
        t_list[val-1] = t[i]

    return t_list


def create_filename(program, args):
    """Create the filename string.

    Parameters
    ----------
    program: str
        The name of the program.
    args: dict
        The arguments parsed to the program.

    Returns
    -------
    filename: str
        The filename string.
    """

    # read input arguments
    mod = args['model']
    t = fa.read_t_from_file() if args['input'] else args['t']
    lat = args['lattice']
    alpha = args['alpha']
    theta = args['theta']
    save = args['save']
    log = args['log']

    t_str = "t_" + '_'.join([f"{i:g}" for i in t]) + "_"
    brav_str = f"alpha_{alpha:g}_theta_{theta[0]:g}_{theta[1]:g}_" if lat not in ["square", "triangular"] else ""

    if program == "band_structure":
        samp = args['samp']
        wil = args['wilson']
        disp = args['display']
        nphi = args['nphi']
        bgt = args['bgt']

        disp_str = f"{disp}_"
        mod_str = f"{mod}_" if mod != "Hofstadter" else ""
        nphi_str = f"nphi_{nphi[0]}_{nphi[1]}_"
        bgt_str = f"bgt_{bgt:g}_"
        samp_str = f"samp_{samp:g}_" if samp != 101 else ""

        filename = f"band_structure_{disp_str}{mod_str}{lat}_{nphi_str}{t_str}{brav_str}{bgt_str}{samp_str}"[:-1]  # remove last underscore

    elif program == "butterfly":
        plt_lat = args["plot_lattice"]
        q = args['q']
        color = args['color']
        pal = args['palette']
        wan = args['wannier']
        period = args['periodicity']
        art = args['art']
        dpi = args['dpi']

        col_str = f"col_{color}_{pal}_" if color else ""
        per_str = f"period_{period:g}_" if period != 1 else ""
        art_str = "art_" if art else ""
        dpi_str = f"dpi_{dpi:g}_" if dpi != 300 else ""

        filename = f"butterfly_{lat}_q_{q:g}_{t_str}{brav_str}{col_str}{per_str}{art_str}{dpi_str}"[:-1]  # remove last underscore
    else:
        raise ValueError("program is not defined")

    return filename


def save_data(program, model, args, data):
    """Save data to file.

    Parameters
    ----------
    program: str
        The name of the program.
    model: Hofstadter.hamiltonian
        The Hamiltonian class attribute.
    args: dict
        The arguments parsed to the program.
    data: ndarray
        The data array.
    """

    filename = create_filename(program, args)
    save_list = [model, args, data]
    directory = "../data/" if os.path.isdir('../data') else ""
    np.save(directory+filename, save_list)

    return None


class Logger(object):
    """Stream stdout and stderr to file."""

    def __init__(self, program, args):
        self.terminal = sys.stdout or sys.stderr
        directory = "../logs/" if os.path.isdir('../logs') else ""
        filename = create_filename(program, args)
        self.log = open(directory+filename+".log", 'w', buffering=1)

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        # this flush method is needed for python 3 compatibility.
        # this handles the flush command by doing nothing.
        # you might want to specify some extra behavior here.
        pass
