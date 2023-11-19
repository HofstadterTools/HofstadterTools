import matplotlib.pyplot as plt
import functions.utility as fu
import functions.plotting as fp


if __name__ == '__main__':

    # load the file
    filename = "band_structure_3D_square_nphi_1_4_t_1.npz"
    model, args, data = fu.load_data("band_structure", filename)

    # overwrite args parameters
    # args['save'] = False

    # construct the figure(s)
    fp.band_structure(model, args, data)
    plt.show()
