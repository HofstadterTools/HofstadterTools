"""Functions for argument parsing."""

# --- external imports
import argparse


def parse_input_arguments(program, description):
    """
    Parse the input arguments for a given program.

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

    general = parser.add_argument_group("general options")
    band_structure = parser.add_argument_group("band_structure options")
    butterfly = parser.add_argument_group("butterfly options")

    if program in ["plot_band_structure", "plot_butterfly"]:
        general.add_argument("-load", type=str, required=True, help="filename of data to load")
    if program in ["band_structure", "butterfly"]:
        models = ["Hofstadter"]
        general.add_argument("-mod", "--model", type=str, default="Hofstadter", choices=models, help="name of model")
        general.add_argument("-a", type=float, default=1, help="lattice constant")
        general.add_argument("-t", nargs='+', type=float, default=[1], help="list of hopping amplitudes in ascending order [1NN, 2NN, 3NN, ...]")
        lattices = ["square", "triangular", "bravais", "honeycomb", "kagome", "custom"]
        general.add_argument("-lat", "--lattice", type=str, default="bravais", choices=lattices, help="name of lattice")
        general.add_argument("-alpha", type=float, default=1, help="length of a2 Bravais vector relative to a1 (Bravais lattice anisotropy)")
        general.add_argument("-theta", nargs=2, type=int, default=[1, 3], help="angle between Bravais basis vectors as a fraction of pi (Bravais lattice obliqueness)")
        general.add_argument("-input", default=False, action='store_true', help="read hopping parameters from hopping_input.txt file")
        general.add_argument("-period", "--periodicity", type=int, default=1, help="factor by which to divide A_UC in the flux density")
        general.add_argument("-log", default=False, action='store_true', help="save the output logs")
    if program in ["band_structure", "butterfly", "plot_band_structure", "plot_butterfly"]:
        general.add_argument("-save", default=False, action='store_true', help="save the output data and plots")
        general.add_argument("-plt_lat", "--plot_lattice", default=False, action='store_true', help="plot the lattice")
        general.add_argument("-dpi", type=int, default=300, help="dots-per-inch resolution of the saved output image")
        general.add_argument("-ps", "--point_size", type=float, default=1, help="scale factor by which to scale the default point size")

    if program == "band_structure":
        band_structure.add_argument("-samp", type=int, default=101, help="number of samples in linear direction")
        band_structure.add_argument("-nphi", nargs=2, type=int, default=[1, 4], help="flux density")
        band_structure.add_argument("-bgt", type=float, default=0.01, help="band gap threshold")
        band_structure.add_argument("-load", type=str, default=False, help="filename of data to load")
        band_structure.add_argument("-topo", "--topology", default=False, action='store_true', help="print the topology columns")
        band_structure.add_argument("-geom", "--geometry", default=False, action='store_true', help="print the basic quantum geometry columns")
        columns = ["band", "group", "isolated", "width", "gap", "gap_width", "std_B", "C", "std_g", "av_gxx", "std_gxx", "av_gxy", "std_gxy", "T", "D"]
        band_structure.add_argument("-cols", "--columns", nargs='+', type=str, default=["band", "group", "isolated", "width", "gap", "gap_width"], choices=columns, help="select the table columns individually")
    if program in ["band_structure", "plot_band_structure"]:
        band_structure.add_argument("-wil", "--wilson", default=False, action='store_true', help="plot the wilson loops")
        displays = ["3D", "2D", "both"]
        band_structure.add_argument("-disp", "--display", type=str, default="3D", choices=displays, help="how to display band structure")

    if program == "butterfly":
        butterfly.add_argument("-q", type=int, default=199, help="denominator of flux density (prime integer)")
    if program in ["butterfly", "plot_butterfly"]:
        colors = [False, "point", "plane"]
        butterfly.add_argument("-col", "--color", type=str, default=False, choices=colors, help="how to color the Hofstadter butterfly")
        palettes = ["avron", "jet", "red-blue"]
        butterfly.add_argument("-pal", "--palette", type=str, default="avron", choices=palettes, help="color palette")
        butterfly.add_argument("-wan", "--wannier", default=False, action='store_true', help="plot the Wannier diagram")
        butterfly.add_argument("-art", default=False, action='store_true', help="remove all plot axes and labels")

    args = vars(parser.parse_args())

    return args
