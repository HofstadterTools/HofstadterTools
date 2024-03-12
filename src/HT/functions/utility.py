"""Functions for input/output and miscellaneous use."""

# --- external imports
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


def create_filename(program, args, aux_text=""):
    """Create the filename string.

    Parameters
    ----------
    program: str
        The name of the program.
    args: dict
        The arguments passed to the program.
    aux_text: str
        The auxiliary text passed to the filename.

    Returns
    -------
    filename: str
        The filename string.
    """

    # read input arguments
    mod = args['model']
    a = args['a']
    t = read_t_from_file() if args['input'] else args['t']
    lat = args['lattice']
    alpha = args['alpha']
    theta = args['theta']
    save = args['save']
    log = args['log']
    period = args['periodicity']
    dpi = args['dpi']

    aux_str = aux_text if aux_text == "" else aux_text+"_"
    a_str = f"a_{a:g}_" if a != 1 else ""
    t_str = "t_" + '_'.join([f"{i:g}" for i in t]) + "_"
    brav_str = f"alpha_{alpha:g}_theta_{theta[0]:g}_{theta[1]:g}_" if lat not in ["square", "triangular"] else ""
    per_str = f"period_{period:g}_" if period != 1 else ""
    dpi_str = f"dpi_{dpi:g}_" if dpi != 300 else ""

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

        filename = f"band_structure_{aux_str}{disp_str}{mod_str}{lat}_{nphi_str}{a_str}{t_str}{brav_str}{per_str}{samp_str}{dpi_str}"[:-1]

    elif program == "butterfly":
        plt_lat = args["plot_lattice"]
        q = args['q']
        color = args['color']
        pal = args['palette']
        wan = args['wannier']
        art = args['art']

        q_str = f"q_{q:g}_"
        col_str = f"col_{color}_{pal}_" if color else ""
        art_str = "art_" if art else ""

        filename = f"butterfly_{aux_str}{lat}_{q_str}{a_str}{t_str}{brav_str}{col_str}{per_str}{art_str}{dpi_str}"[:-1]

    else:
        raise ValueError("program is not defined.")

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
        The arguments passed to the program.
    data: ndarray
        The data array.
    """

    directory = f"../../data/{program}/" if os.path.isdir(f"../../data/{program}/") else ""
    filename = create_filename(program, args)
    np.savez_compressed(directory+filename, model=model, args=args, data=data)

    return None


def load_data(program, filename, plotting=False):
    """Load data from file.

    Parameters
    ----------
    program: str
        The name of the program.
    filename: str
        The name of the file to be loaded.
    plotting: bool
        The plotting script flag.
    """

    rel_path = "../../.." if plotting else "../.."
    directory = f"{rel_path}/data/{program}/" if os.path.isdir(f"{rel_path}/data/{program}/") else ""
    file_data = np.load(directory+filename, allow_pickle=True)
    model = file_data['model'].item()  # .item() unpacks 0-dim array
    args = file_data['args'].item()
    data = file_data['data'].item()

    return model, args, data


class Logger(object):
    """Stream stdout and stderr to file."""

    def __init__(self, program, args):
        self.terminal = sys.stdout or sys.stderr
        directory = f"../../logs/{program}/" if os.path.isdir(f"../../logs/{program}/") else ""
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
