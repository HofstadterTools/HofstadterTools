import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from copy import deepcopy


def nearest_neighbor_finder(avec, acartvec, abasisvec, t_list, x_init, y_init):

    # --- Create list of NN to consider from t_list
    numb_list = []
    for i, t in enumerate(t_list):
        if t != 0:
            numb_list.append(i+1)

    # --- Create grid of basis vectors from [-t_max, t_max]
    vectors = []
    vectors_unit = []
    for i in range(-numb_list[-1], numb_list[-1]+1):
        for j in range(-numb_list[-1], numb_list[-1]+1):
            r_unit = np.array([i, j])
            vectors_unit.append(r_unit)
            r = np.matmul(r_unit, avec)
            for idx, k in enumerate(abasisvec):
                vectors.append([r+k, idx])  # add basis index

    # --- Shift the grid of vectors relative to an initial point (x_init, y_init)
    # vectors = [np.subtract(i, np.array([x_init, y_init])) for i in vectors]
    shift_vectors = []
    for i, val in enumerate(vectors):
        shift_vectors.append([np.subtract(val[0], np.array([x_init, y_init])), val[1]])

    # --- Define data array with info on each vector
    data = np.zeros((len(shift_vectors), 13), dtype=object)
    for i, val in enumerate(shift_vectors):
        data[i, 0] = round(np.linalg.norm(val[0]), 10)  # round so that we can use it for comparison
        data[i, 1] = np.angle(val[0][0]+1j*val[0][1])
        data[i, 2] = val[0][0]
        data[i, 3] = val[0][1]
        data[i, 4] = None  # m
        data[i, 5] = val[0][1] / avec[1][1]  # n
        data[i, 10] = val[1]

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
                data[i, 6] = j+1  # NN group
                data[i, 7] = None  # m_init
                data[i, 8] = y_init / avec[1][1]  # n_init
                data[i, 9] = data[i, 8] + data[i, 5]  # n_total
                data[i, 11] = x_init
                data[i, 12] = y_init
    select_radii = [radii[i - 1] for i in numb_list]
    # delete rows with other radii
    rows_to_delete = []
    for i, row in enumerate(data):
        if row[0] not in select_radii:
            rows_to_delete.append(i)
    data = np.delete(data, rows_to_delete, axis=0)

    # --- Compute change sublattice flags
    # if len(abasisvec) > 1:  # if there is a multi-particle basis
    #     for i, val in enumerate(data):
    #         for j, val2 in enumerate(abasisvec):
    #             if j > 0:  # skip origin vector
    #                 dx = round(abs(val[2]), 10)
    #                 dy = round(abs(val[3]), 10)
    #                 if round(dx % acartvec[0][0], 10) == 0 and round(dy % avec[1][1], 10) == 0:
    #                     data[i, 10] = 0
    #                 elif round(dx % val2[0], 10) == 0 and round(dy % val2[1], 10) == 0:
    #                     data[i, 10] = j  # label sublattice by number
    #                 else:
    #                     data[i, 10] = None  # point not on any sublattice
    # print("data = ", data)
    # 1/0

    return data


def nearest_neighbor_sorter(data_array):

    # --- Count backtrack vectors
    delete_list = []
    backtrack_list = []
    for i, val in enumerate(data_array):
        if round(val[11] + val[2], 10) == 0 and round(val[12] + val[3], 10) == 0:  # backtrack
            delete_list.append(i)
            backtrack_list.append(val[6])
    data_array = np.delete(data_array, delete_list, axis=0)
    # print("filter data array = ", data_array)

    # --- Count the independent paths
    numb_single_paths, numb_double_paths = 0, 0
    # double paths
    double_idx = []
    for i, val in enumerate(data_array):
        if round(val[11]) == 0 and round(val[12]) == 0:  # origin
            for j, val2 in enumerate(data_array):
                if round(val[11] + val[2], 10) == round(val2[11], 10) and round(val[12] + val[3], 10) == round(val2[12], 10):
                    numb_double_paths = numb_double_paths + 1
                    double_idx.append(i)
                    double_idx.append(j)
    # single paths
    double_idx_list = np.sort(list(set(double_idx)))
    numb_single_paths = len(data_array) - len(double_idx_list)

    # --- Group the independent paths
    single_paths = np.zeros(numb_single_paths, dtype=object)
    for i in range(numb_single_paths):
        single_paths[i] = []
    double_paths = np.zeros(numb_double_paths, dtype=object)
    for i in range(numb_double_paths):
        double_paths[i] = []
    # double paths
    counter = 0
    for i, val in enumerate(data_array):
        if round(val[11]) == 0 and round(val[12]) == 0:  # origin
            for j, val2 in enumerate(data_array):
                if round(val[11] + val[2], 10) == round(val2[11], 10) and round(val[12] + val[3], 10) == round(val2[12], 10):
                    ntot = val[5] + val2[5]
                    double_paths[counter] = [val, val2, round(ntot)]
                    counter = counter + 1
    # single paths
    single_idx_list = np.setdiff1d(np.arange(len(data_array)), np.array(double_idx_list))
    for i, val in enumerate(single_idx_list):
        single_paths[i] = [data_array[val], round(data_array[val][9])]

    # print("single paths = ", single_paths)
    # print("double paths = ", double_paths)

    paths_list = [single_paths, double_paths]
    grouped_paths_list = []

    for paths in paths_list:

        # --- Count the number of total n
        n_tot = []
        for i, val in enumerate(paths):
            n_tot.append(val[-1])
        n_tot = np.sort(list(set(n_tot)))
        # print("n_tot = ", n_tot)
        numb_n = len(n_tot)
        # print("numb_n = ", numb_n)

        # --- Group the independent paths by total n
        grouped_paths = np.zeros(numb_n, dtype=object)
        for i in range(numb_n):
            grouped_paths[i] = []

        for i, nval in enumerate(n_tot):
            for j, val in enumerate(paths):
                if val[-1] == nval:
                    grouped_paths[i].append(val[:-1])
            grouped_paths[i].append(nval)

        grouped_paths_list.append(grouped_paths)

    # print("grouped single paths = ", grouped_paths_list[0])
    # print("grouped double paths = ", grouped_paths_list[1])

    return grouped_paths_list, backtrack_list


def peierls_factor(A_UC_val, nphi, xy_vec, delta_y, n):

    phase = - 2 * np.pi * nphi * xy_vec[0] * (n + delta_y/2) / A_UC_val
    factor = np.exp(1j * phase)

    return factor


def rammal_factor(nphi, nval, ninit, ntot, cos_theta, ay, k_val_val):

    phase = k_val_val[1]*ay*(ntot - ninit) - np.pi * nphi * cos_theta * ((nval + ntot)**2 - (nval + ninit)**2)
    factor = np.exp(1j * phase)

    return factor


def diag_func(A_UC_val, t_val, p_val, q_val, vec_list, n_val, k_val_val, i_val, cos_ang, backtrack_list, ay):
    nphi = p_val/q_val
    term = 0
    for idx, val in enumerate(vec_list):
        if val[-1] == i_val:  # extract rows with appropriate ntot
            ninits = []
            for k, val2 in enumerate(val[:-1]):  # for each path group
                factor = 1
                for val3 in val2:  # for each vector in path group (usually 2 for honeycomb)
                    ninits.append(val3[8])
                    NN_group = int(val3[6])
                    xy_vector = np.array([val3[2], val3[3]])
                    delta_n = val3[5]
                    factor *= - t_val[NN_group - 1] * (peierls_factor(A_UC_val, nphi, xy_vector, delta_n, n_val + val3[8])
                              * np.exp(1j * np.vdot(xy_vector, k_val_val)))
                term += factor
            abs_ninits = [abs(j) for j in ninits]
            term = term * rammal_factor(nphi, n_val, ninits[np.argmin(abs_ninits)], val[-1], cos_ang, ay, k_val_val)
    if i_val == 0:
        for i in backtrack_list:
            term = term + t_val[i-1]**2  # t^2 factor for every backtrack vector
    return term


def Hamiltonian(t_val, p_val, q_val, A_UC_val, vec_group_list_list, k_val, cos_angle, backtrack_list, ay):

    Hamiltonian_list = []

    for vec_group_list in vec_group_list_list:
        if vec_group_list != []:
            Hamiltonian = np.zeros((q_val, q_val), dtype=np.complex128)

            n_values = []
            for term in vec_group_list:
                n_values.append(term[-1])
            pos_n_vals_comb = [i for i in n_values if i >= 0]

            for i in pos_n_vals_comb:
                # upper_diag_array
                upper_diag_array = np.array([diag_func(A_UC_val, t_val, p_val, q_val, vec_group_list, n % q_val, k_val, i, cos_angle, backtrack_list, ay) for n in range(q_val)])
                Hamiltonian += np.roll(np.diag(upper_diag_array), i, axis=1)
                # lower_diag_array
                if i > 0:
                    Hamiltonian += np.roll(np.diag(np.conj(upper_diag_array)), i, axis=0)

            Hamiltonian_list.append(Hamiltonian)
        else:
            Hamiltonian_list.append(np.zeros((q_val, q_val)))

    return Hamiltonian_list


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
