import matplotlib.pyplot as plt
import functions.utility as fu
import functions.plotting as fp


if __name__ == '__main__':

    # load the file
    filename = "butterfly_kagome_q_499_t_1_alpha_1_theta_1_3_period_8_dpi_1500.npz"
    model, args, data = fu.load_data("butterfly", filename)

    # overwrite args parameters
    # args['save'] = False
    args['dpi'] = 600

    # construct the figure(s)
    fp.butterfly(model, args, data)
    plt.show()
