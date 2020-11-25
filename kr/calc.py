"""
this file is responsible for the actual calculations
to solve the differential equations for N_i
"""

import numpy as np
import scipy.linalg as la

from numba import jit


class Calc:
    @staticmethod
    def I_f(I_0, k_I, t):
        return I_0 * np.exp(- k_I * t)

    @staticmethod
    def lamda_0_f(eta, k_I, I_, k_tc, k_td):
        return np.sqrt((2. * eta * k_I * I_) / (k_tc + k_td))

    @staticmethod
    def M_f(M_0, eta, k_td, k_I, t, k_tc, I_0, k_p, k_fm):
        M = M_0 * np.exp((2./k_I) * (k_p + k_fm) * np.power(
            (2. * eta * k_I * I_0 / (k_tc + k_td)), (1./2.)) *
            (np.exp(-(k_I/2.) * t) - 1.))
        return M

    @staticmethod
    def mu_1_f(M_0, M):
        return M_0 - M

    @staticmethod
    def X_f(M_0, M):
        return (M_0 - M) / M_0

    @staticmethod
    def calc_Ndot_vec_QSSA_triangular(N_vec, M, mu_1, k_fp, lamda_0, pivots, k_fm, eta, k_I, I_,
                                      k_tc, k_td, k_pstar, k_p,
                                      b_p_vec):
        y_vec = -k_fp * lamda_0 * pivots * N_vec  # d_tilde_vec
        y_vec[0] += -(k_fm * M * lamda_0 + 2. * eta * k_I * I_)  # e_tilde_vec

        # --- build QSSA matrix
        a_tilde_vec = -(k_fm * M +
                        (k_tc + k_td) * lamda_0 +
                        (k_fp + k_pstar) * mu_1)
        c_tilde_vec = -k_p * M * b_p_vec
        b_tilde_vec = (k_p * M * b_p_vec)[:-1]
        M_mat = np.diag(a_tilde_vec + c_tilde_vec) + np.diag(b_tilde_vec, -1)
        Ndot_vec = la.solve_triangular(M_mat, y_vec, lower=True)
        return Ndot_vec

    @staticmethod
    def calc_derivative(t, y,
                        pivots,  # params set by solver.set_f_params()
                        M_0, I_0, eta,
                        k_I, k_p,
                        k_tc, k_fm, k_td, k_fp, k_pstar,
                        a_p_vec, b_p_vec,
                        imax,
                        tc_method_handler,
                        slice_manager):

        # get_all_slices() returns a list of sub-vectors (slices)
        [N_vec] = slice_manager.get_all_slices(y)
        # N_vec, Ndot_vec = slice_manager.get_all_slices(y)
        # get analytical quantities
        I_ = Calc.I_f(I_0, k_I, t)
        lamda_0 = Calc.lamda_0_f(eta, k_I, I_, k_tc, k_td)
        M = Calc.M_f(M_0, eta, k_td, k_I, t, k_tc, I_0, k_p, k_fm)
        mu_1 = Calc.mu_1_f(M_0, M)

        # -- it would be possible to integrate M and I_ numerically
        # N_vec, M, I_ = slice_manager.get_slices_from_vec(y)
        # dM_dt_rhs = -(k_fm + k_p) * M * lamda_0
        # dmu_1_dt_rhs = k_p * M * lamda_0 (vgl. Butte Thesis, P.17)
        # dI_dt_rhs = - k_I * I_

        # bulid matrix and solve triangular system
        # with scipy.linalg.solve_triangular
        Ndot_vec = Calc.calc_Ndot_vec_QSSA_triangular(
            N_vec, M, mu_1,
            k_fp, lamda_0, pivots, k_fm, eta, k_I, I_,
            k_tc, k_td, k_pstar, k_p,
            b_p_vec)

        # ----- RHS of dNdt (27) ----

        rhs_term1_vec = -(k_fp + k_pstar) * lamda_0 * pivots * N_vec

        # N_vec at time t is provided by the ode solver
        # Ndot_vec at time t is provided by the QSSA system of equations
        rhs_term2_vec = (
            (k_fm * M +
             k_td * lamda_0 +
             k_fp * mu_1) * Ndot_vec)

        rhs_tc_term_vec = np.zeros(imax + 1)
        if k_tc != 0.0:
            rhs_tc_term_vec = (
                0.5 * k_tc * tc_method_handler.calc_tc_vec(Ndot_vec, imax))

        dNdt_vec = rhs_term1_vec + rhs_term2_vec + rhs_tc_term_vec
        return dNdt_vec

    # --- Termination by combination (tc) calculation functions:
    # -- not numpy-idiomatic (explicit python loops)

    @staticmethod
    @jit
    def calc_tc_CM_LOOPS(m_k, m_a, Ndot_vec, imax):
        """ Makes use of the connection matrix produced by gen_tc_CM_LOOPS;
            runs fast because of @jit;
            nomenclature: (i, j) -> (k, k+1)
            parameters: m_k contains k_ij, and -1 at invalid points
                        m_a contains corresponding a_ij, and 0 at invalid points

            This procedure was originally written by M. Rommel and adapted by C. Schnur
        """
        sums = np.zeros(len(Ndot_vec))
        for i in range(imax + 1):
            for j in range(0, imax + 1):
                k = m_k[i, j]
                if k == -1:  # if sum condition is not matched (see gen_tc)
                    continue
                a = m_a[i, j]
                b = 1.0 - a
                if i == j:
                    a *= 0.5
                    b *= 0.5
                sums[k] += a * Ndot_vec[i] * Ndot_vec[j]
                assert (k + 1 <= imax)
                sums[k + 1] += b * Ndot_vec[i] * Ndot_vec[j]
        return sums

    # --- quantities to plot
    @staticmethod
    def dispersity(N_vec_slice, pivots_slice):
        """ calculates the dispersity according to the text """
        return ((np.sum(N_vec_slice * np.power(pivots_slice, 2.)) / np.sum(N_vec_slice * pivots_slice)) /
                (np.sum(N_vec_slice * pivots_slice) / np.sum(N_vec_slice)))

    @staticmethod
    def gelfraction(mu_1, pivots, N_vec, imax):
        mu_1_sol = np.sum(
            pivots[:imax + 1] * N_vec[:imax + 1])
        return (mu_1 - mu_1_sol) / mu_1

    @staticmethod
    @jit
    def reconstruct(N_vectors, pivots, imax):
        P_vectors = []
        for N_vec in N_vectors:
            P_vec = np.zeros(imax + 1)

            for i in range(imax - 1):
                P_vec[i] = N_vec[i] / (pivots[i + 1] - pivots[i])
            P_vectors.append(P_vec)

        return P_vectors

    # --- misc

    @staticmethod
    def mu_1_f_naive_sum(pivots, N_vec):
        # only usable if a large number of chains is simulated
        # when chains get infinite, it is unusable
        return np.sum(pivots * N_vec)

    @staticmethod
    def sum_weight(N_vec_slice, pivots_slice):
        return np.sum(N_vec_slice * pivots_slice)

    @staticmethod
    def cartesian_product(x, y):
        """
        cartesian prodct {x} x {y}
        credit: https://stackoverflow.com/a/11144716
        """
        return np.transpose([np.tile(x, len(y)), np.repeat(y, len(x))])
