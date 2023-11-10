"""Functions for the model classes."""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from copy import deepcopy


def nearest_neighbor_finder(avec, abasisvec, t_list, x_init, y_init, basis_init):

    # --- Create list of NN to consider from t_list
    numb_list = []
    for i, t in enumerate(t_list):
        if t != 0:
            numb_list.append(i+1)

    # --- Create grid of basis vectors from [-numb_t_max, numb_t_max]
    abasisvec = [abasisvec[i] - abasisvec[basis_init] for i in range(len(abasisvec))]  # shift basis
    vectors = []
    vectors_unit = []
    for i in range(-numb_list[-1], numb_list[-1]+1):
        for j in range(-numb_list[-1], numb_list[-1]+1):
            r_unit = np.array([i, j])
            vectors_unit.append(r_unit)
            r = np.matmul(r_unit, avec)
            for idx, k in enumerate(abasisvec):
                vectors.append([r+k, i, j, basis_init, idx])  # add basis index

    # --- Define data array with info on each vector
    data = np.zeros((len(vectors), 13), dtype=object)
    for i, val in enumerate(vectors):
        data[i, 0] = round(np.linalg.norm(val[0]), 10)  # r (round so that we can use it for comparison)
        data[i, 1] = np.angle(val[0][0]+1j*val[0][1])  # phi
        data[i, 2] = x_init  # x_init
        data[i, 3] = y_init  # y_init
        data[i, 4] = y_init / avec[1][1]  # y_init / a2_y
        data[i, 5] = val[0][0]  # dx
        data[i, 6] = val[0][1]  # dy
        data[i, 7] = val[0][1] / avec[1][1]  # dy / a2_y
        data[i, 8] = val[1]  # dI
        data[i, 9] = val[2]  # dJ
        data[i, 10] = val[3]  # sub_init
        data[i, 11] = val[4]  # sub_final

    # --- Extract the NN groups (filter data based on radius)
    data = data[data[:, 0].argsort()]  # sort by increasing r
    # delete the first r=0 row
    mask = (data[:, 0] != 0)
    data = data[mask, :]
    radii = np.sort(list(set(data[:, 0])))
    # label the NN group
    for i, row in enumerate(data):
        for j, r in enumerate(radii):
            if row[0] == r:
                data[i, 12] = j+1  # NN group
    select_radii = [radii[i - 1] for i in numb_list]
    # delete rows with other radii
    rows_to_delete = []
    for i, row in enumerate(data):
        if row[0] not in select_radii:
            rows_to_delete.append(i)
    data = np.delete(data, rows_to_delete, axis=0)

    # --- Extract bases set
    bases = []
    for i, val in enumerate(data):
        bases.append(val[10])
        bases.append(val[11])
    bases = np.sort(list(set(bases)))

    return data, bases


def nearest_neighbor_sorter(data_array):

    # count the number of dJ
    dJ_list = []
    for i, val in enumerate(data_array):
        dJ_list.append(val[9])
    dJ_list = np.sort(list(set(dJ_list)))
    numb_dJ = len(dJ_list)

    # group paths by dJ
    grouped_paths = np.zeros(numb_dJ, dtype=object)
    for i in range(numb_dJ):
        grouped_paths[i] = []

    for i, dJval in enumerate(dJ_list):
        for j, val in enumerate(data_array):
            if val[9] == dJval:
                grouped_paths[i].append(val)
        grouped_paths[i].append(dJval)

    return grouped_paths


def peierls_factor(nphi, dx, y_cart, dy_cart, A_UC):

    phase = - 2 * np.pi * nphi * dx * (y_cart + dy_cart/2) / A_UC
    factor = np.exp(1j * phase)

    return factor


def diag_func(t_val, p_val, q_val, A_UC_val, vec_group, k_val, dJ_val, J_idx_val):

    nphi = p_val/q_val
    term = 0
    for idx, val in enumerate(vec_group):
        if val[-1] == dJ_val:  # extract rows with appropriate dJ
            for k, val2 in enumerate(val[:-1]):  # for each vector in path group
                NN_group = int(val2[12])
                term += - t_val[NN_group - 1] * (peierls_factor(nphi, val2[5], J_idx_val + val2[4], val2[7], A_UC_val)
                                                 * np.exp(1j * np.vdot(np.array([val2[5], val2[6]]), k_val)))

    return term


def Hamiltonian(t, p, q, A_UC, vec_group_matrix, k):

    I = np.shape(vec_group_matrix)[0]
    J = np.shape(vec_group_matrix)[1]

    Ham_matrix = []
    for i in range(I):
        Ham_row = []
        for j in range(J):
            Hamiltonian = np.zeros((q, q), dtype=np.complex128)

            dJ_list = []
            for term in vec_group_matrix[i, j]:
                dJ_list.append(term[-1])

            HC_flag = False
            for k1, val in enumerate(dJ_list):  # remove negative unit cell hoppings (for H.c. cases)
                for k2, val2 in enumerate(dJ_list):
                    if val == -val2 and val != 0:
                        HC_flag = True
                        if val2 < 0:
                            del dJ_list[k2]
                        else:
                            del dJ_list[k1]

            for dJ in dJ_list:
                # upper_diag_array
                diag_array = np.array([diag_func(t, p, q, A_UC, vec_group_matrix[i, j], k, dJ, J_idx) for J_idx in range(q)])
                Hamiltonian += np.roll(np.diag(diag_array), abs(dJ), axis=int((np.sign(dJ)+1)/2))
                # lower_diag_array
                if HC_flag and dJ > 0:
                    Hamiltonian += np.roll(np.diag(np.conj(diag_array)), abs(dJ), axis=0)

            Ham_row.append(Hamiltonian)
        Ham_matrix.append(Ham_row)
    Hamiltonian = np.block(Ham_matrix)

    return Hamiltonian


if __name__ == '__main__':

    t = [1, 0, -0.25]

    vec_group = nearest_neighbor_finder("square", t)

    print(vec_group)

    ###

    num_bands = 5
    num_samples = 101

    b1 = (2. * np.pi) * np.array([1 / num_bands, 0])
    b2 = (2. * np.pi) * np.array([0, 1])
    bvec = np.vstack((b1, b2))

    eigenvalues = np.zeros((num_bands, num_samples, num_samples))  # real
    eigenvectors = np.zeros((num_bands, num_bands, num_samples, num_samples), dtype=np.complex128)  # complex
    for band in range(num_bands):
        for idx_x in range(num_samples):
            frac_kx = idx_x / (num_samples - 1)
            for idx_y in range(num_samples):
                frac_ky = idx_y / (num_samples - 1)
                k = np.matmul(np.array([frac_kx, frac_ky]), bvec)
                # print("k = ", k)
                eigvals, eigvecs = np.linalg.eigh(Hamiltonian(t, 1, num_bands, vec_group, k))
                idx = np.argsort(eigvals)
                eigenvalues[band, idx_x, idx_y] = eigvals[idx[band]]
                eigenvectors[:, band, idx_x, idx_y] = eigvecs[:, idx[band]]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    idx_x = np.linspace(0, num_samples - 1, num_samples - 1, dtype=int)
    idx_y = np.linspace(0, num_samples - 1, num_samples - 1, dtype=int)
    kx, ky = np.meshgrid(idx_x, idx_y)
    for i in range(num_bands):
        ax.plot_surface(kx, ky, eigenvalues[i, kx, ky], alpha=0.5)
    ax.set_xlabel('$k_1/|\mathbf{b}_1|$')
    ax.set_ylabel('$k_2/|\mathbf{b}_2|$')
    ax.set_zlabel('$E$')

    def normalize(value, tick_number):
        if value == 0:
            return "$0$"
        elif value == num_samples - 1:
            return "$1$"
        else:
            return f"${value / (num_samples - 1):.1g}$"

    ax.xaxis.set_major_formatter(plt.FuncFormatter(normalize))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(normalize))

    plt.show()
