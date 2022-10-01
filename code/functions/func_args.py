import argparse


def parse_input_arguments(program):

    parser = argparse.ArgumentParser(prog=program)
    models = ["Hofstadter"]
    parser.add_argument("-mod", "--model", type=str, default="Hofstadter", choices=models, help="name of model")
    displays = ["3D", "2D"]
    parser.add_argument("-disp", "--display", type=str, default="3D", choices=displays, help="how to display band structure")
    parser.add_argument("-nphi", nargs=2, type=int, default=[1, 4], help="flux density")
    parser.add_argument("-samp", type=int, default=101, help="number of samples in linear direction")
    args = vars(parser.parse_args())

    return args
