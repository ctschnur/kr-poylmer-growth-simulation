from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
from __future__ import unicode_literals

import numpy as np

from numba import jit


class Gen:
    """
    Contains functions used to generate resources (e.g. lists of indices or
    connection matrices) that are used in later calculation to reduce computation time
    """
    # --- not numpy-idiomatic
    @staticmethod 
    @jit
    def gen_tc_CM_LOOPS(pivots, i_max):
        """ Generates Connection Matrix as described in the text for TbC 
            nomenclature: (i, j) -> (k)
            return: m_k contains k_ij, m_a contains corresponding a_ij
        """
        m_k = np.zeros((i_max + 1, i_max + 1), dtype=np.int32)
        m_a = np.zeros((i_max + 1, i_max + 1), dtype=np.float64)
        m_k += -1  # k=-1 means: condition not matched, continue
        for i in range(i_max + 1):
            for j in range(i_max + 1):
                for k in range(j, i_max - 1):
                    if ((pivots[k] < pivots[i] + pivots[j]) and
                            (pivots[i] + pivots[j] <= pivots[k + 1])):
                        m_k[i, j] = k
                        m_a[i, j] = ((pivots[k + 1] - (pivots[i] + pivots[j])) /
                                     (pivots[k + 1] - pivots[k]))
                        break
        return m_k, m_a


    # --- pivot generation


    @staticmethod 
    def gen_pivots_method_1(imax, nbar_imax):
        """ primary method for generating pivots;
            takes the index as pivot for s^i < i """
        s = np.power(nbar_imax, 1. / imax)
        pivots = s ** np.arange(1, imax + 1)
        # import ipdb; ipdb.set_trace()  # noqa BREAKPOINT
        # print(pivots[:5], ", ...", pivots[-3:])
        
        pivots = np.append(0., pivots)  # prepend 0 as first pivot
        pivots = np.rint(pivots)  # round to nearest integer
        for i in np.arange(0, imax + 1):  # clamping the values
            if pivots[i] < i:
                pivots[i] = i
        
        return pivots


    # --- alternative pivot generation method

    @staticmethod
    def gen_pivots_method_2(imax, nbar_imax):
        """ alternative method for generating pivots;
            if s_m doesn't genertate a unique list large enough, 
            iteratively decrease s_m and generate larger lists, until
            the unique list is large enough """  
        s_m = None
        m = 0  # iteratively determine s_m by increasing counter
        pivots_m = np.array([])

        while True:
            s_m = nbar_imax ** (1. / (imax + m))
            pivots_m = np.unique(np.rint(s_m ** np.arange(imax + m + 1)))

            if pivots_m.size >= imax:
                break
            
            m += 1

        pivots = pivots_m[-imax:]  # get imax of the biggest values
        pivots = np.append(0., pivots_m)  # prepend zero as the first pivot

        assert False not in (pivots == np.unique(pivots))  # check for doubles
        print("pivots generated: s: ", s_m, ", m: ", m)
        # import ipdb; ipdb.set_trace()  # noqa BREAKPOINT
        print(pivots[:5], ", ...", pivots[-3:])
        return pivots


    # --- misc


    @staticmethod 
    def gen_one_more_pivot_for_method12(pivots):
        """ compatible with gen_pivots_method_1 and gen_pivots_method_2
            (because they generate exponentially increasing values);
            generates one more pivot (after main pivot generation has ended)
            in order to calculate the last propagation sharing coefficient 
            b_i_max. The calculation of b_i_max needs pivots[i_max + 1], which 
            would be out of bounds if not explicitly calculated. Not using 
            b_i_max (effectively setting it to 0) would result in an artificial 
            accumulation of mass in the generation imax of active chains """
        max_idx = pivots.size - 1  # i starts at 0
        x = np.power(pivots[-1], 1. / np.float64(max_idx))  # recalculate factor

        cur_i = max_idx + 1
        cur_nbar = np.power(x, cur_i)

        # making sure that all generated values are unique after
        # rounding to the *nearest integer* with numpy.rint
        if(cur_nbar - cur_i > 0.5):
            new_pivot = cur_nbar
        else:
            new_pivot = cur_i

        pivots_extended = np.rint(np.append(pivots, new_pivot))

        assert (np.unique(pivots_extended).size == pivots_extended.size)

        return np.rint(new_pivot)


    @staticmethod 
    def gen_ab_p_vecs(pivots, i_p_max, pivot_at_imaxplus1):
        """ generates sharing coefficients a and b for propagation 
            Parameter pivot_at_imaxplus1: Explanation see gen_one_more_pivot() """
        # generate b_p_vec of all available pivots
        ks = np.arange(i_p_max + 1, dtype=np.int32)
        pivots_ = np.append(pivots, pivot_at_imaxplus1)

        a_p_vec = np.zeros(i_p_max + 1)
        a_p_vec[ks] += ((pivots_[(ks + 1)] - (pivots_[ks] + 1.)) /
                        (pivots_[(ks + 1)] - pivots_[ks]))

        b_p_vec = np.zeros(i_p_max + 1)
        b_p_vec[ks] += 1. / (pivots_[(ks + 1)] - pivots_[ks])

        # Debugging: check that all a, b are all in [0., 1.] and a + b = 1
        assert (False not in (np.logical_and((a_p_vec >= 0.0), (a_p_vec <= 1.0))))
        assert (False not in (np.logical_and((b_p_vec >= 0.0), (b_p_vec <= 1.0))))
        assert (False not in (a_p_vec + b_p_vec == 1.0))

        return a_p_vec, b_p_vec
