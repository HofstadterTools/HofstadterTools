import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d


def nearest_neighbor_finder(lattice, numb):

    if lattice == "square":
        vectors = []
        for i in range(-numb, numb+1):
            for j in range(-numb, numb+1):
                r = np.array([i, j])
                vectors.append(r)

    data = np.zeros((len(vectors), 4), dtype=object)
    for i, r in enumerate(vectors):
        data[i, 0] = np.linalg.norm(r)
        data[i, 1] = np.angle(r[0]+1j*r[1])
        data[i, 2] = r[0]
        data[i, 3] = r[1]

    # data = data[data[:, 0].argsort()]  # sort by increasing r

    # delete the first r=0 row
    mask = (data[:, 0] != 0)
    data = data[mask, :]
    # delete rows that exceed the NN number
    mask = (data[:, 0] <= numb)
    data = data[mask, :]
    # delete conjugate hoppings
    # mask = (0 <= data[:, 1])
    # data = data[mask, :]
    # mask = (data[:, 1] < np.pi)
    # data = data[mask, :]

    data = data[data[:, 2].argsort()]  # sort by m offset

    return data


def peierls_factor(nphi, vec, m):

    phase = 2 * np.pi * nphi * vec[1] * (m + vec[0]/2)
    factor = np.exp(1j * phase)

    return factor


def diag_func(nphi, vec_list, k, m):

    term = 0
    for vec in vec_list:
        term += peierls_factor(nphi, vec, m) * np.exp(1j * np.vdot(vec, k))

    return term


def Hamiltonian(p_val, q_val, vec_group_list, k_val):

    Hamiltonian = np.zeros((q_val, q_val), dtype=np.complex128)

    for m in range(q_val):
        Hamiltonian[m][m] = - diag_func(p_val / q_val, vec_group_list[1], k_val, m)

    for m in range(q_val-1):
        Hamiltonian[m][m+1] = - diag_func(p_val / q_val, vec_group_list[2], k_val, m)
        Hamiltonian[m+1][m] = - diag_func(p_val / q_val, vec_group_list[0], k_val, m+1)

    Hamiltonian[0][q_val-1] = - diag_func(p_val / q_val, vec_group_list[0], k_val, 0)
    Hamiltonian[q_val-1][0] = - diag_func(p_val / q_val, vec_group_list[2], k_val, q_val-1)

    # print(Hamiltonian)
    # 1/0

    return Hamiltonian


if __name__ == '__main__':

    data2 = nearest_neighbor_finder("square", 1)
    print(data2)

    min_term = np.min(data2[:, 2])
    max_term = np.max(data2[:, 2])

    vec_group = np.zeros(max_term - min_term + 1, dtype=object)
    for i, vec_group_idx in enumerate(range(min_term, max_term + 1)):
        vec_group[i] = []

    for i, vec_group_idx in enumerate(range(min_term, max_term + 1)):
        for row in data2:
            if row[2] == vec_group_idx:
                vec_group[i].append(np.array([row[2], row[3]]))

    print(vec_group[0])
    print(vec_group[1])
    print(vec_group[2])

    # print(vec_group)

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
                eigvals, eigvecs = np.linalg.eigh(Hamiltonian(1, num_bands, vec_group, k))
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
