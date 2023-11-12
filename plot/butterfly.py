import matplotlib.pyplot as plt
import functions.utility as fu
import functions.plotting as fp


if __name__ == '__main__':

    # load the file
    filename = "butterfly_square_q_399_t_1_col_plane_avron.npy"
    model, args, data = fu.load_data("butterfly", filename)

    # overwrite args parameters
    # args['save'] = False
    # args['dpi'] = 1000

    # construct the figure(s)
    fp.butterfly(model, args, data)
    plt.show()
