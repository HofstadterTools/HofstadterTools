import numpy as np
import os
import matplotlib.pyplot as plt
import functions.utility as fu
import functions.plotting as fp


if __name__ == '__main__':

    # load the file
    filename = "butterfly_honeycomb_q_99_t_1_alpha_1_theta_1_3_col_point_avron.npy"
    dir = "../data/" if os.path.isdir('../data') else ""
    data = np.load(dir+filename, allow_pickle=True)

    # read/customize the data
    model = data[0]
    args = data[1]
    nphi_list = data[2]
    E_list = data[3]
    E_list_orig = data[4]
    chern_list = data[5]
    matrix = data[6]
    nphi_DOS_list = data[7]
    DOS_list = data[8]
    gaps_list = data[9]
    tr_DOS_list = data[10]

    # plot the butterfly
    fp.butterfly(model, args, nphi_list, E_list, E_list_orig, chern_list, matrix, nphi_DOS_list, DOS_list, gaps_list, tr_DOS_list)

    plt.show()
