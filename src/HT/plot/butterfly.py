# --- external imports
import matplotlib.pyplot as plt
# --- internal imports
from HT.functions import utility as fu
from HT.functions import plotting as fp


if __name__ == '__main__':

    # load the file
    filename = "butterfly_bravais_q_199_t_1_alpha_1_theta_1_3.npz"
    model, args, data = fu.load_data("butterfly", filename, True)

    # overwrite args parameters
    # args['save'] = False
    # args['dpi'] = 600

    # construct the figure(s)
    fp.butterfly(model, args, data, True)
    plt.show()
