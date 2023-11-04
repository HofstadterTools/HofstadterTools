"""Functions for argument parsing."""

import argparse


def parse_input_arguments(program):
    """
    Parse the input arguments to a given program.

    Each program may be run with a given set of flags. This function parses those command-line arguments and returns them as a dictionary.

    Parameters
    ----------
    program: string
        The name of the program.

    Returns
    -------
    args: dict
        The dictionary of input arguments.
    """

    parser = argparse.ArgumentParser(prog=program)
    models = ["Hofstadter"]
    parser.add_argument("-mod", "--model", type=str, default="Hofstadter", choices=models, help="name of model")
    parser.add_argument("-t", nargs='+', type=float, default=[1], help="list of hopping amplitudes in order of ascending NN")
    lattices = ["square", "triangular", "bravais", "honeycomb"]
    parser.add_argument("-lat", "--lattice", type=str, default="bravais", choices=lattices, help="name of lattice")
    parser.add_argument("-alpha", type=float, default=1, help="length of a2 basis vector relative to a1 (Bravais lattice anisotropy)")
    parser.add_argument("-theta", nargs=2, type=float, default=[1, 3], help="angle between Bravais basis vectors (in units of pi)")

    if program == "band_structure":
        parser.add_argument("-samp", type=int, default=101, help="number of samples in linear direction")
        parser.add_argument("-wil", "--wilson", default=False, action='store_true', help="plot the wilson loops")
        displays = ["3D", "2D"]
        parser.add_argument("-disp", "--display", type=str, default="3D", choices=displays, help="how to display band structure")
        parser.add_argument("-nphi", nargs=2, type=int, default=[1, 4], help="flux density")
        parser.add_argument("-bgt", type=float, default=0.01, help="band gap threshold")
    if program == "butterfly":
        parser.add_argument("-q", type=int, default=199, help="denominator of flux density (prime)")
        colors = [False, "point", "plane"]
        parser.add_argument("-col", "--color", type=str, default=False, choices=colors, help="how to color the Hofstadter butterfly")

    args = vars(parser.parse_args())

    return args
