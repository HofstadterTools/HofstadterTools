"""Functions for argument parsing."""

import argparse


def parse_input_arguments(program, description):
    """
    Parse the input arguments to a given program.

    Each program may be run with a given set of flags. This function parses those command-line arguments and returns them as a dictionary.

    Parameters
    ----------
    program: string
        The name of the program.
    description: string
        The description of the program.

    Returns
    -------
    args: dict
        The dictionary of input arguments.
    """

    parser = argparse.ArgumentParser(prog=program, description=description, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    models = ["Hofstadter"]
    parser.add_argument("-mod", "--model", type=str, default="Hofstadter", choices=models, help="name of model")
    parser.add_argument("-t", nargs='+', type=float, default=[1], help="list of hopping amplitudes in ascending order [1NN, 2NN, 3NN, ...]")
    lattices = ["square", "triangular", "bravais", "honeycomb", "kagome", "custom"]
    parser.add_argument("-lat", "--lattice", type=str, default="bravais", choices=lattices, help="name of lattice")
    parser.add_argument("-alpha", type=float, default=1, help="length of a2 Bravais vector relative to a1 (Bravais lattice anisotropy)")
    parser.add_argument("-theta", nargs=2, type=float, default=[1, 3], help="angle between Bravais basis vectors as a fraction of pi (Bravais lattice obliqueness)")
    parser.add_argument("-input", default=False, action='store_true', help="read hopping parameters from hopping_input.txt file")
    parser.add_argument("-save", default=False, action='store_true', help="save the output data and plots")
    parser.add_argument("-log", default=False, action='store_true', help="save the output logs")
    parser.add_argument("-plt_lat", "--plot_lattice", default=False, action='store_true', help="plot the lattice")
    parser.add_argument("-period", "--periodicity", type=int, default=1, help="factor by which to divide A_UC in the flux density")
    parser.add_argument("-dpi", type=int, default=300, help="dots-per-inch resolution of the saved output image")

    if program == "band_structure":
        parser.add_argument("-samp", type=int, default=101, help="number of samples in linear direction")
        parser.add_argument("-wil", "--wilson", default=False, action='store_true', help="plot the wilson loops")
        displays = ["3D", "2D", "both"]
        parser.add_argument("-disp", "--display", type=str, default="3D", choices=displays, help="how to display band structure")
        parser.add_argument("-nphi", nargs=2, type=int, default=[1, 4], help="flux density")
        parser.add_argument("-bgt", type=float, default=0.01, help="band gap threshold")
        parser.add_argument("-load", type=str, default=False, help="load data from file")
    if program == "butterfly":
        parser.add_argument("-q", type=int, default=199, help="denominator of flux density (prime integer)")
        colors = [False, "point", "plane"]
        parser.add_argument("-col", "--color", type=str, default=False, choices=colors, help="how to color the Hofstadter butterfly")
        palettes = ["avron", "jet", "red-blue"]
        parser.add_argument("-pal", "--palette", type=str, default="avron", choices=palettes, help="color palette")
        parser.add_argument("-wan", "--wannier", default=False, action='store_true', help="plot the Wannier diagram")
        parser.add_argument("-art", default=False, action='store_true', help="remove all plot axes and labels")

    args = vars(parser.parse_args())

    return args
