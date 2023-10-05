import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d


def nearest_neighbor_finder(avec, acartvec, t_list):

    # print("acart", acartvec)

    numb_list = []  # create list of NN to consider from t_list
    for i, t in enumerate(t_list):
        if t != 0:
            numb_list.append(i+1)

    vectors = []
    vectors_unit = []
    for i in range(-numb_list[-1], numb_list[-1]+1):
        for j in range(-numb_list[-1], numb_list[-1]+1):
            r_unit = np.array([i, j])
            vectors_unit.append(r_unit)
            r = np.matmul(r_unit, avec)
            vectors.append(r)

    data = np.zeros((len(vectors), 7), dtype=object)
    for i, r in enumerate(vectors):
        data[i, 0] = round(np.linalg.norm(r), 10)  # round so that we can use it for comparison
        data[i, 1] = np.angle(r[0]+1j*r[1])
        data[i, 2] = r[0]
        data[i, 3] = r[1]
        data[i, 4] = round(r[0] / acartvec[0][0])
        data[i, 5] = round(r[1] / acartvec[1][1])

    # --- 1) Extract the NN groups

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

    select_radii = [radii[i - 1] for i in numb_list]

    rows_to_delete = []
    for i, row in enumerate(data):
        if row[0] not in select_radii:
            rows_to_delete.append(i)

    data = np.delete(data, rows_to_delete, axis=0)

    # --- 2) Sort into vec groups

    data = data[data[:, 4].argsort()]  # sort by m offset

    m_values = np.sort(list(set(data[:, 4])))

    vec_group = np.zeros(len(m_values), dtype=object)
    for i in range(len(m_values)):
        vec_group[i] = []

    for i, vec_group_idx in enumerate(m_values):
        for row in data:
            if row[4] == vec_group_idx:
                vec_group[i].append((np.array([row[2], row[3]]), np.array([row[4], row[5]]), row[6]))

    #print(vec_group)

    return vec_group


def peierls_factor(nphi, vec, m):

    phase = 2 * np.pi * nphi * vec[1] * (m + vec[0]/2)
    factor = np.exp(1j * phase)

    return factor


def diag_func(t_val, nphi, vec_list, k, m):

    term = 0
    for vec_inf in vec_list:
        term += t_val[int(vec_inf[2])-1] * peierls_factor(nphi, vec_inf[1], m) * np.exp(1j * np.vdot(vec_inf[0], k))

    return term


def Hamiltonian(t_val, p_val, q_val, acartvec, vec_group_list, k_val):

    Hamiltonian = np.zeros((q_val, q_val), dtype=np.complex128)

    m_values = []
    for term in vec_group_list:
        m_values.append(term[0][1][0])

    pos_m_vals = [i for i in m_values if i >= 0]

    for i in pos_m_vals:
        if i == 0:
            Hamiltonian += np.diag(np.array([- diag_func(t_val, p_val / q_val, vec_group_list[len(vec_group_list)//2], k_val, m) for m in range(q_val)]))
        else:
            offset = 0 if 0 in pos_m_vals else -1
            # upper_diag_array
            upper_diag_array = np.array([- diag_func(t_val, p_val / q_val, vec_group_list[len(vec_group_list)//2+offset+i], k_val, (m+i) % q_val) for m in range(q_val)])
            Hamiltonian += np.roll(np.diag(upper_diag_array), i, axis=1)
            # lower_diag_array
            lower_diag_array = np.array([- diag_func(t_val, p_val / q_val, vec_group_list[len(vec_group_list)//2-i], k_val, (m+i) % q_val) for m in range(q_val)])
            Hamiltonian += np.roll(np.diag(lower_diag_array), i, axis=0)

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
