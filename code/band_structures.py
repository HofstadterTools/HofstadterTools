# --- python imports
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

a0 = 1  # lattice constant
t = 1  # hopping amplitude


def define_unit_cell(model_val, q_val=4):

    if 'Squ' in model_val:
        if model_val is 'HofSqu1':
            num_bands_val = q_val
            # lattice vectors
            a1 = a0 * np.array([num_bands_val, 0])
            a2 = a0 * np.array([0, 1])
            avec_val = np.vstack((a1, a2))
            # reciprocal lattice vectors
            b1 = (2. * np.pi) / a0 * np.array([1 / num_bands_val, 0])
            b2 = (2. * np.pi) / a0 * np.array([0, 1])
            bvec_val = np.vstack((b1, b2))
            # symmetry points
            GA = np.array([0, 0])
            Y = np.array([0, 0.5])
            S = np.array([0.5, 0.5])
            X = np.array([0.5, 0])
            sym_points_val = [GA, Y, S, X]
        else:
            return ValueError("Requested model is not implemented in define_unit_cell function.")
    else:
        return ValueError("Requested lattice cannot be read from model name.")

    return num_bands_val, avec_val, bvec_val, sym_points_val


def hamiltonian(model_val, k_val, num_bands_val, p_val=1):

    # initialize the Hamiltonian
    Hamiltonian = np.zeros((num_bands_val, num_bands_val), dtype=np.complex128)

    if 'Squ' in model_val:

        # nearest neighbors
        delta = np.zeros((2, 2))
        delta[0, :] = a0 * np.array([1, 0])
        delta[1, :] = a0 * np.array([0, 1])

        if model_val is 'HofSqu1':

            q_val = num_bands_val
            nphi = p_val/q_val

            def h(k_val_val, m_val):
                return 2 * np.cos(2 * np.pi * nphi * m_val + k_val_val[1] * a0)

            for n in range(q_val):
                Hamiltonian[n][n] = t * h(k_val, n)

            for n in range(q_val-1):
                Hamiltonian[n][n+1] = t * np.exp(+1j*k_val[0]*a0)
                Hamiltonian[n+1][n] = t * np.exp(-1j*k_val[0]*a0)

            Hamiltonian[0][q_val-1] = t * np.exp(-1j*k_val[0]*a0)
            Hamiltonian[q_val-1][0] = t * np.exp(+1j*k_val[0]*a0)

        else:
            return ValueError("Requested model is not implemented in hamiltonian function.")

    else:
        return ValueError("Requested lattice cannot be read from model name.")

    return Hamiltonian


def berry_curv(eigvec, eigvec_alpha, eigvec_beta, eigvec_alpha_beta):

    Berry_curv = - np.imag(np.log(np.conj(eigvec).dot(eigvec_alpha) * np.conj(eigvec_alpha).dot(eigvec_alpha_beta)
                           * np.conj(eigvec_alpha_beta).dot(eigvec_beta) * np.conj(eigvec_beta).dot(eigvec)))

    return Berry_curv


if __name__ == '__main__':

    # initialization ###################################################################################################
    num_samples = 101
    model = 'HofSqu1'  # (HofSqu1, HalSquCN), (HalTriC3), (graphene, HalHexC1)

    flag_3D = False  # choose between 3D or 2D band structure
    p, q = 1, 5  # for Hofstadter model only
    ####################################################################################################################

    if flag_3D:
        # define unit cell
        num_bands, avec, bvec, sym_points = define_unit_cell(model, q_val=q)

        # construct bands
        eigenvalues = np.zeros((num_bands, num_samples, num_samples))  # real
        eigenvectors = np.zeros((num_bands, num_bands, num_samples, num_samples), dtype=np.complex128)  # complex
        for band in range(num_bands):
            for idx_x in range(num_samples):
                frac_kx = idx_x / (num_samples-1)
                for idx_y in range(num_samples):
                    frac_ky = idx_y / (num_samples-1)
                    k = np.matmul(np.array([frac_kx, frac_ky]), bvec)
                    eigvals, eigvecs = np.linalg.eig(hamiltonian(model, k, num_bands, p_val=p))
                    idx = np.argsort(eigvals)
                    eigenvalues[band][idx_x][idx_y] = np.real(eigvals[idx[band]])
                    eigenvectors[:, band, idx_x, idx_y] = eigvecs[:, idx[band]]

        # compute Chern numbers
        berry_fluxes = np.zeros((num_bands, num_samples - 1, num_samples - 1))  # real
        for band in range(num_bands):
            for idx_x in range(num_samples - 1):
                for idx_y in range(num_samples - 1):
                    berry_fluxes[band, idx_x, idx_y] = berry_curv(eigenvectors[:, band, idx_x, idx_y],
                                                                  eigenvectors[:, band, idx_x + 1, idx_y],
                                                                  eigenvectors[:, band, idx_x, idx_y + 1],
                                                                  eigenvectors[:, band, idx_x + 1, idx_y + 1])
        chern_numbers = np.zeros(num_bands)
        for band in range(num_bands):
            chern_numbers[band] = np.sum(berry_fluxes[band, :, :]) / (2 * np.pi)
            print(f"Chern number ({band}) = {chern_numbers[band]}")
            if band == 0:
                # print(f"Average Berry curvature ({band}) = {np.average(berry_fluxes[band, :, :])}")
                # print(f"Standard deviation Berry curvature ({band}) = {np.std(berry_fluxes[band, :, :])}")
                print(f"Stdev/|average| Berry curvature ({band}) = {np.std(berry_fluxes[band, :, :])/np.abs(np.average(berry_fluxes[band, :, :]))}")

        # construct figure
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        idx_x = np.linspace(0, num_samples - 1, num_samples, dtype=int)
        idx_y = np.linspace(0, num_samples - 1, num_samples, dtype=int)
        kx, ky = np.meshgrid(idx_x, idx_y)
        for i in range(num_bands):
            ax.plot_surface(kx, ky, eigenvalues[i, kx, ky])
        ax.set_xlabel('$k_x$')
        ax.set_ylabel('$k_y$')
        ax.set_zlabel('$E$')
    else:
        # define unit cell
        num_bands, avec, bvec, sym_points = define_unit_cell(model, q_val=q)

        # construct bands
        num_paths = len(sym_points)
        points_per_path = int(num_samples/num_paths)
        num_points = num_paths*points_per_path
        eigenvalues = np.zeros((num_bands, num_points))  # real
        count = 0
        for i in range(num_paths):
            for j in range(points_per_path):
                k = sym_points[i] + (sym_points[(i+1) % num_paths] - sym_points[i]) * float(j) / float(points_per_path - 1)
                k = np.matmul(k, bvec)
                eigvals = np.linalg.eigvals(hamiltonian(model, k, num_bands, p_val=p, tx_factor_val=1))
                idx = np.argsort(eigvals)
                for band in range(num_bands):
                    eigenvalues[band, count] = np.real(eigvals[idx[band]])
                count += 1

        # construct figure
        fig = plt.figure()
        ax = plt.subplot(111)
        for i in range(num_bands):
            ax.plot(eigenvalues[i])
        for i in range(1, num_paths):
            ax.axvline(i*points_per_path, color='k', linewidth=0.5, ls='--')
        ax.set_xlim([0, num_points])
        ax.set_xlabel('symmetry path')
        ax.set_ylabel('$E$')

    # band analysis
    band_gap = np.min(eigenvalues[1]) - np.max(eigenvalues[0])
    band_width = np.max(eigenvalues[0]) - np.min(eigenvalues[0])
    print(f"band width (0) = {band_width}")
    print(f"band gap (0-1) = {band_gap}")
    print(f"gap-to-width ratio (0-1) = {band_gap / band_width}")

    plt.show()
