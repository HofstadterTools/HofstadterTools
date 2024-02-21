# --- external imports
import matplotlib.pyplot as plt
# --- internal imports
from HT.functions import arguments as fa
from HT.functions import utility as fu
from HT.functions import plotting as fp


def main():
    # parse input arguments
    new_args = fa.parse_input_arguments("plot_band_structure", "Replot the Hofstadter Band Structure.")

    # load the file
    filename = new_args['load']
    model, old_args, data = fu.load_data("band_structure", filename, True)

    # overwrite old_args with new_args parameters, if specified
    old_args['save'] = new_args['save']
    old_args['plot_lattice'] = new_args['plot_lattice']
    old_args['dpi'] = new_args['dpi']
    old_args['point_size'] = new_args['point_size']

    old_args['wilson'] = new_args['wilson']
    old_args['display'] = new_args['display']

    # construct the figure(s)
    fp.band_structure(model, old_args, data, True)
    plt.show()


if __name__ == '__main__':

    main()
