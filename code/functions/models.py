import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from copy import deepcopy


def nearest_neighbor_finder(avec, acartvec, abasisvec, t_list, m_init, n_init):

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
            for k in abasisvec:
                vectors.append(r+k)

    # --- Shift the grid of vectors relative to an initial point (m_init, n_init)
    vectors = [np.subtract(i, np.array(m_init * acartvec[0] + n_init * acartvec[1])) for i in vectors]

    # --- Define data array with info on each vector
    data = np.zeros((len(vectors), 11), dtype=object)
    for i, r in enumerate(vectors):
        data[i, 0] = round(np.linalg.norm(r), 10)  # round so that we can use it for comparison
        data[i, 1] = np.angle(r[0]+1j*r[1])
        data[i, 2] = r[0]
        data[i, 3] = r[1]
        data[i, 4] = round(r[0] / acartvec[0][0])
        data[i, 5] = round(r[1] / acartvec[1][1])

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
                data[i, 6] = j+1
                data[i, 7] = m_init
                data[i, 8] = n_init
                data[i, 9] = m_init + data[i, 4]
    select_radii = [radii[i - 1] for i in numb_list]
    # delete rows with other radii
    rows_to_delete = []
    for i, row in enumerate(data):
        if row[0] not in select_radii:
            rows_to_delete.append(i)
    data = np.delete(data, rows_to_delete, axis=0)

    # --- Compute change sublattice flags
    if len(abasisvec) > 1:  # if there is a multi-particle basis
        # define unit cell in mn units
        x_unit = int(avec[0][0]/acartvec[0][0])
        y_unit = int(avec[1][1]/acartvec[1][1])
        # define avecs in mn units
        avec_mn = deepcopy(avec)
        acartx = np.linalg.norm(acartvec[0])
        acarty = np.linalg.norm(acartvec[1])
        for i, vec in enumerate(avec):
            avec_mn[i] = np.array([vec[0]/acartx, vec[1]/acarty])
        # define abasis in mn units
        abasis_mn = deepcopy(abasisvec)[1:]
        for i, vec in enumerate(abasisvec):
            if i > 0:
                abasis_mn[i-1] = np.array([vec[0]/acartx, vec[1]/acarty])
        # define x and y identifiers
        x_identifiers, y_identifiers = [], []
        for i in range(x_unit):
            if i not in avec_mn[:, 0] and i % abasis_mn[:, 0] == 0 and i != 0:
                x_identifiers.append(i)
        for j in range(y_unit):
            if j not in avec_mn[:, 1] and j % abasis_mn[:, 1] == 0 and j != 0:
                y_identifiers.append(j)
        # print("x, y = ", x_identifiers, y_identifiers)

        for i, val in enumerate(data):
            if ((np.abs(val[7] + val[4]) % x_unit in x_identifiers)
                    or (np.abs(val[8] + val[5]) % y_unit in y_identifiers)):
                data[i, 10] = 1
    # print("data = ", data)

    return data


def nearest_neighbor_sorter(data_array):

    # --- Delete backtrack vectors
    delete_list = []
    for i, val in enumerate(data_array):
        if val[7] + val[4] == 0 and val[8] + val[5] == 0:  # backtrack
            delete_list.append(i)
    data_array = np.delete(data_array, delete_list, axis=0)
    # print("filter data array = ", data_array)

    # --- Count the independent paths
    numb_paths = 0
    # count double paths
    double = True
    for i, val in enumerate(data_array):
        if val[7] == 0 and val[8] == 0:  # origin
            for j, val2 in enumerate(data_array):
                if val[7] + val[4] == val2[7] and val[8] + val[5] == val2[8]:
                    numb_paths = numb_paths + 1
    # --- Count single paths
    if numb_paths == 0:
        double = False
        for i, val in enumerate(data_array):
            if val[7] == 0 and val[8] == 0:  # origin
                numb_paths = numb_paths + 1
    # print("numb_paths = ", numb_paths)

    # --- Group the independent paths
    paths = np.zeros(numb_paths, dtype=object)
    for i in range(numb_paths):
        paths[i] = []
    counter = 0
    if double:
        for i, val in enumerate(data_array):
            if val[7] == 0 and val[8] == 0:  # origin
                for j, val2 in enumerate(data_array):
                    if val[7] + val[4] == val2[7] and val[8] + val[5] == val2[8]:
                        if np.max([np.abs(val[9]), np.abs(val2[9])]) == np.abs(val2[9]):
                            mtot = val2[9]
                        else:
                            mtot = val[9]
                        paths[counter] = [val, val2, mtot]
                        counter = counter + 1
    else:
        for i, val in enumerate(data_array):
            if val[7] == 0 and val[8] == 0:  # origin
                paths[counter] = [val, val[9]]
                counter = counter + 1
    # print("paths = ", paths)

    # --- Count the number of total m
    m_tot = []
    for i, val in enumerate(paths):
        m_tot.append(val[-1])
    m_tot = np.sort(list(set(m_tot)))
    # print("m_tot = ", m_tot)
    numb_m = len(m_tot)
    # print("numb_m = ", numb_m)

    # --- Group the independent paths by total m
    grouped_paths = np.zeros(numb_m, dtype=object)
    for i in range(numb_m):
        grouped_paths[i] = []

    for i, mval in enumerate(m_tot):
        for j, val in enumerate(paths):
            if val[-1] == mval:
                grouped_paths[i].append(val[:-1])
        grouped_paths[i].append(mval)
    # print("grouped_paths = ", grouped_paths)

    return grouped_paths


def peierls_factor(A_UC_val, nphi, vec, m):

    phase = 2 * np.pi * nphi * vec[1] * (m + vec[0]/2) / A_UC_val
    factor = np.exp(1j * phase)

    return factor


def diag_func(A_UC_val, t_val, p_val, q_val, vec_list, m_val, k_val_val, i_val):
    nphi = p_val/q_val
    term = 0
    for idx, val in enumerate(vec_list):
        if val[-1] == i_val:  # extract rows with appropriate mtot
            for k, val2 in enumerate(val[:-1]):  # for each path group
                factor = 1
                for val3 in val2:  # for each vector in path group (usually 2 for honeycomb)
                    NN_group = int(val3[6])
                    xy_vector = np.array([val3[2], val3[3]])
                    mn_vector = np.array([val3[4], val3[5]])
                    factor *= - t_val[NN_group - 1] * (peierls_factor(A_UC_val, nphi, mn_vector, (m_val + val3[7]) % (q_val+1))
                              * np.exp(1j * np.vdot(xy_vector, k_val_val)))
                term += factor
    return term


def Hamiltonian(t_val, p_val, q_val, A_UC_val, vec_group_list, k_val):

    Hamiltonian = np.zeros((q_val, q_val), dtype=np.complex128)

    m_values = []
    for term in vec_group_list:
        m_values.append(term[-1])
    pos_m_vals_comb = [i for i in m_values if i >= 0]

    for i in pos_m_vals_comb:
        # upper_diag_array
        upper_diag_array = np.array([diag_func(A_UC_val, t_val, p_val, q_val, vec_group_list, (m+i) % q_val, k_val, i) for m in range(q_val)])
        Hamiltonian += np.roll(np.diag(upper_diag_array), i, axis=1)
        # lower_diag_array
        if i > 0:
            Hamiltonian += np.roll(np.diag(np.conj(upper_diag_array)), i, axis=0)

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
