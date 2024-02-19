# --- external imports
import matplotlib.pyplot as plt
# --- internal imports
from HT.functions import utility as fu
from HT.functions import plotting as fp


if __name__ == '__main__':

    # load the file
    filename = "band_structure_3D_bravais_nphi_1_4_t_1_alpha_1_theta_1_3.npz"
    model, args, data = fu.load_data("band_structure", filename, True)

    # overwrite args parameters
    # args['save'] = False

    # construct the figure(s)
    fp.band_structure(model, args, data, True)
    plt.show()
