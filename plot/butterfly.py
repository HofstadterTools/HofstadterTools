import matplotlib.pyplot as plt
import functions.utility as fu
import functions.plotting as fp


if __name__ == '__main__':

    # load the file
    filename = "butterfly_honeycomb_q_99_t_1_alpha_1_theta_1_3_col_plane_avron.npy"
    model, args, data = fu.load_data("butterfly", filename)

    # construct the figure(s)
    fp.butterfly(model, args, data)
    plt.show()
