""" Kr class: Management of importing, calculation and exporting """

import numpy as np
import scipy.integrate as scint
import json
import time

from kr.calc import Calc
from kr.gen import Gen


class Kr_constants:
    def __init__(self, Pdict):
        self._original_Pdict = Pdict
        self.import_params(Pdict)

    def import_params(self, Pdict):
        """ for each variable present in Json, manually create a
            member variable with the appropriate value (JSON automatically 
            imports the right data types (if correctly specified in 
            JSON file) """
        self._original_Pdict = Pdict

        self.M_0 = Pdict["M_0"]
        self.I_0 = Pdict["I_0"]
        self.k_p = Pdict["k_p"]
        self.k_tc = Pdict["k_tc"]
        self.k_fm = Pdict["k_fm"]
        self.k_pstar = Pdict["k_pstar"]
        self.k_I = Pdict["k_I"]
        self.k_td = Pdict["k_td"]
        self.k_fp = Pdict["k_fp"]
        self.eta = Pdict["eta"]

        self.imax = Pdict["imax"]
        self.nbar_imax = Pdict["nbar_imax"]

        self.t_step = Pdict["t_step"]
        self.conversions_to_save_clds_at = Pdict[
            "conversions_to_save_clds_at"]

    def reset(self, Pdict=None):
        if Pdict is None:
            self.import_params(self._original_Pdict)
        else:
            self.Pdict = Pdict

    @staticmethod
    def complete_dict(default_dict, custom_dict):
        """ take items that are missing from custom dict out of the default dict 
        """
        for dkey, dval in default_dict.items():
            if dkey not in custom_dict.keys():
                custom_dict[dkey] = dval

        return custom_dict

    @staticmethod
    def import_dict_from_file(json_file_path, custom_dict_name="DEFAULT"):
        """ import parameters as a dictionary from a json file """
        with open(json_file_path) as f:
            data = json.load(f)

        if "DEFAULT" not in data.keys():
            print("warning: no dict named DEFAULT found in ", json_file_path)

        if custom_dict_name not in data.keys():
            print("warning: no dict named ", custom_dict_name,
                  " found in ", json_file_path)
            exit(1)

        return Kr_constants.complete_dict(
            data["DEFAULT"], data[custom_dict_name])

    @staticmethod
    def compare_dicts(default_dict, dict2):
        dict1 = default_dict.copy()
        # first make all lists tuples

        def make_lists_tuples(dict_to_convert):
            ret_dict = dict_to_convert.copy()
            for key, value in ret_dict.items():
                if isinstance(value, list):
                    ret_dict[key] = tuple(value)
            return ret_dict

        set1 = set(make_lists_tuples(dict1).items())
        set2 = set(make_lists_tuples(dict2).items())

        # return (set1 ^ set2)  # symmetric, shows all differences
        return (set1 - set2)  # asymmetric, shows elemenst of set1 \ set2

    def get_str(self):
        """ report all currently initialized values """
        s = ""
        s += "\tM_0 " + str(self.M_0) + "\tI_0 " + str(self.I_0) + "\n"
        s += "\tk_I " + str(self.k_I) + "\tk_p " + str(self.k_p) + "\n"
        s += "\tk_tc " + str(self.k_tc) + "\tk_td " + str(self.k_td) + "\n"
        s += "\tk_fm " + str(self.k_fm) + "\tk_fp " + str(self.k_fp) + "\n"
        s += "\tk_pstar " + str(self.k_pstar) + "\n"
        s += "\timax " + str(self.imax) + "\tnbar_imax " + \
            str(self.nbar_imax) + "\n"
        s += "\tt_step " + str(self.t_step) + "\n"
        s += "\tconversions_to_save_clds_at " + \
            str(self.conversions_to_save_clds_at) + "\n"
        s += "\teta " + str(self.eta)
        return s


class Kr:
    def __init__(self):
        print("Kr object created")

    def gen_resources(self):
        """ calculates: 
            - pivots
            - propagation sharing coefficients
            - termination by combination (tc) resources (e.g. connection
              matrix) """

        print("generating pivots \t\t\t(", self.gen_pivots_func.__name__, ")")
        # if self.gen_pivots_func is not None:
        self.pivots = self.gen_pivots_func(self.d.imax, self.d.nbar_imax)
        # else:
        #    print("Err: gen_pivots_func is not assigned")

        if self.d.k_p != 0.:
            print("k_p != 0. -> generating resources \t(",
                  Gen.gen_ab_p_vecs.__name__, ")")
            self.a_p_vec, self.b_p_vec = Gen.gen_ab_p_vecs(
                self.pivots, self.d.imax,
                Gen.gen_one_more_pivot_for_method12(self.pivots))
        if self.d.k_tc != 0.:
            print("k_tc != 0. -> generating resources \t(",
                  self.tc_method_handler.gen_tc_function_obj.__name__, ")")
            if self.tc_method_handler is not None:
                self.tc_method_handler.gen_tc_resources(
                    self.pivots, self.d.imax)
            else:
                print("no tc_method_handler given")
                exit()

    def run(self, imported_params, gen_pivots_func=Gen.gen_pivots_method_1,
            tc_method_handler=None):
        """ manages each run of the simulation 
            parameters: imported_params: Kr_constant containing a dictionary 
                                         with e.g. imported reaction constants
                        tc_method_handler: new instance specifying what methods
                                           gen_tc... and calc_tc... to use 
        """
        # --- assign some function pointers to use
        self.tc_method_handler = tc_method_handler
        self.gen_pivots_func = gen_pivots_func

        # --- handle imported constants
        self.d = Kr_constants(imported_params)
        # optionally print the imported params
        # print("Loaded constants: \n", self.d.get_str())

        # --- generate resources
        self.gen_resources()

        start = time.time()  # measure execution time
        # --- call the actual numerical calculation
        self.solve_fpbes()
        end = time.time()
        print("execution time of solve_fpbes(): {0:.2f} seconds".format(
            end - start))

        # --- reconstruct all N_vectors that were saved during calculation
        self.P_vectors = Calc.reconstruct(self.N_vectors,
                                          self.pivots, self.d.imax)
        # --- calculate distributions based on the reconstructed values
        # e.g. reaction rate distribution for r_fp

        # plot rate of chain transfer to polymer at specific conversions
        # for different terminated chain lenghts
        # we have m * P_m
        # we need sum_n Rdot_n  (= lamda_0)
        # r_fp = k_fp * m * P_m * Rdot_n
        # sum_n=0^size(P_vec) r_fp = k_fp * m * P_m * sum_n=0^size(P_vec) Rdot_n
        #                          = k_fp * m * P_m * lamda_0

        def r_fp_dead_dist_vec(k_fp, pivots, P_vec, lamda_0):
            return k_fp * pivots * P_vec * lamda_0

        self.r_fp_dead_dist_vec = []
        # iterate over saved conversions
        for idx, P_vec in enumerate(self.P_vectors):
            self.r_fp_dead_dist_vec.append(
                r_fp_dead_dist_vec(
                    self.d.k_fp, self.pivots, P_vec,
                    self.lamda_0_at_cld_conv_vec[idx]))

        # --- save away clds as dictionarys adding printing information
        clds = []
        for idx in range(len(self.d.conversions_to_save_clds_at)):
            clds.append({"conversion": self.d.conversions_to_save_clds_at[idx],
                         "imax": self.d.imax,
                         "pivots": self.pivots,
                         "P_n": self.P_vectors[idx],
                         "r_fp_dead_dist": self.r_fp_dead_dist_vec[idx]
                         })

        # --- this simulation run gets it's own hash, so that all figure files
        #     generated (hash in filenames) can be associated with this run
        print("hash: ", self.get_run_hash_str())

        return self.evolution_data_unit_dict, clds

    def solve_fpbes(self):
        """ runs the actual integration and collects plotting information """
        print("solve_fpbes()")

        # init Data_units to collect data
        self.collect_data_init()

        # init the ode solver
        solver = scint.ode(Calc.calc_derivative)  # assign a function
        solver.set_integrator('lsoda', ixpr=True)

        # init the slice manager (insert order is important)
        self.slice_manager = Slice_manager()
        # N_vec
        self.slice_manager.push_slice(np.zeros(self.d.imax + 1))
        # # Ndot_vec
        # self.slice_manager.push_slice(np.zeros(self.d.imax + 1))
        # # mu_1
        # self.slice_manager.push_slice(np.zeros(1))

        # # I
        # self.slice_manager.push_slice(np.array([self.d.I_0]))

        # -- get the vector from Slice_manager
        solver.set_initial_value(self.slice_manager.get_initial_value_vec())

        # pass all necessary parameters to the numerical solver in addition to
        # tc_method_handler and slice_manager
        solver.set_f_params(  # parameters for Calc.calc_derivative
            self.pivots,
            self.d.M_0, self.d.I_0, self.d.eta,
            self.d.k_I, self.d.k_p,
            self.d.k_tc, self.d.k_fm, self.d.k_td, self.d.k_fp, self.d.k_pstar,
            self.a_p_vec, self.b_p_vec,
            self.d.imax,
            self.tc_method_handler,
            self.slice_manager)

        X = 0.  # initial conversion value
        count = 0  # keep track of how many loops (outside runs) are necessary
        # integrate up until a maximum conversion
        while (solver.successful() and
                X < max(self.d.conversions_to_save_clds_at)):
            solver.integrate(solver.t + self.d.t_step)

            # get_all_slices() returns a list of sub-vectors (slices)
            [N_vec] = self.slice_manager.get_all_slices(solver.y)  # with QSSA
            # N_vec, Ndot_vec = slice_manager.get_all_slices(solver.y)

            # here is where the slice manager can come in handy
            # (note that insertion order must be kept)
            # N_vec, M, mu_1, I_ = self.slice_manager.get_all_slices(solver.y)

            # re-calculate some quantities analytically

            I_ = Calc.I_f(self.d.I_0, self.d.k_I, solver.t)
            lamda_0 = Calc.lamda_0_f(
                self.d.eta, self.d.k_I, I_, self.d.k_tc, self.d.k_td)
            M = Calc.M_f(self.d.M_0, self.d.eta, self.d.k_td, self.d.k_I,
                         solver.t, self.d.k_tc, self.d.I_0, self.d.k_p, self.d.k_fm)
            mu_1 = Calc.mu_1_f(self.d.M_0, M)
            X = Calc.X_f(self.d.M_0, M)

            # collect data at conversion X (time t)
            self.evolution_data_unit_dict["dispersity_until_i_max"].append_xy_values(
                X,
                Calc.dispersity(
                    N_vec, self.pivots))
            self.evolution_data_unit_dict["dispersity"].append_xy_values(
                X, Calc.dispersity(N_vec, self.pivots))
            self.evolution_data_unit_dict["conversion"].append_xy_values(
                solver.t, X)
            self.evolution_data_unit_dict["mu_1"].append_xy_values(
                X, mu_1)
            self.evolution_data_unit_dict["lamda_0"].append_xy_values(
                X, lamda_0)
            self.evolution_data_unit_dict["gelfraction"].append_xy_values(
                X, Calc.gelfraction(mu_1, self.pivots, N_vec, self.d.imax))
            self.collect_distribution_data(
                self.evolution_data_unit_dict["conversion"].y_vector,
                N_vec, lamda_0)

            tc_method_string = ""
            if self.tc_method_handler is not None:
                tc_method_string = self.tc_method_handler.get_str()

            print(("outside state: \tsolver.t = {0:.2f} " +
                   "conversion = {1:.2f} ").format(
                solver.t, X),
                tc_method_string,
                end='\r')
            count += 1

        print("\nsolve_fpbes(): while loop ran ", count, " times")

    def collect_data_init(self):
        """ sets up dictionary where this run's data is collected """

        self.N_vectors = []  # save clds at specific conversions
        self.lamda_0_at_cld_conv_vec = []  # save lamda_0 there
        # index needed in collect_distribution_data to count through the list of preset
        # conversions_to_save_clds_at
        self.cur_conv_save_index = 0

        # --- bundle up data for plotting:
        # make a dictionary that saves some evolution data
        # in so called Data_units (simplify plotting)
        self.evolution_data_unit_dict = {  # save misc quantities
            "gelfraction": Data_unit(
                x_vector=np.array([]),
                y_vector=np.array([]),
                x_name="conversion",
                y_name="gelfraction"),
            "dispersity": Data_unit(
                x_vector=np.array([]),
                y_vector=np.array([]),
                x_name="conversion",
                y_name="dispersity"),
            "conversion": Data_unit(
                x_vector=np.array([]),
                y_vector=np.array([]),
                x_name="time",
                y_name="conversion"),
            "lamda_0": Data_unit(
                x_vector=np.array([]),
                y_vector=np.array([]),
                x_name="conversion",
                y_name="lamda_0"),
            "mu_1": Data_unit(
                x_vector=np.array([]),
                y_vector=np.array([]),
                x_name="conversion",
                y_name="mu_1"),
            "dispersity_until_i_max": Data_unit(
                x_vector=np.array([]),
                y_vector=np.array([]),
                x_name="conversion",
                y_name="dispersity_until_i_max"),
        }

    def collect_distribution_data(self, conversions, N_vec, lamda_0):
        """ collects cld at specific preset conversions """
        if len(conversions) >= 2:
            if ((conversions[-2] < self.d.conversions_to_save_clds_at[
                    self.cur_conv_save_index]) and
                (conversions[-1] >= self.d.conversions_to_save_clds_at[
                    self.cur_conv_save_index])):
                print("\nsaving N_vec at conversion = ",
                      "{0:.2f}".format(conversions[-1]))
                self.N_vectors.append(N_vec)
                self.lamda_0_at_cld_conv_vec.append(lamda_0)
                self.cur_conv_save_index += 1

    def get_oneline_str(self):
        """ provides information about the run (without linebreaks 
            -> for filenames or headings) """
        string = ""
        if self.tc_method_handler is not None:
            string += self.tc_method_handler.get_str() + ";"

        return string

    def get_linebreaks_str(self):
        """ provides information about the run (with linebreaks) """
        return self.get_oneline_str().replace(";", "\n")

    def get_run_hash_str(self):
        """ gets the hash (or at least the last 5 letters) of the kr object
            (to append in front of filenames, making the association between
            the run and the generated plots/csv files) easier """
        return str(self.__hash__())[-5:]


class Tc_method_handler:
    """ Kr.run() is passed an instance of Tc_method_handler, which handles 
        the selection of different gen_tc and calc_tc methods automatically. 
        Also, different gen_tc methods can return a different number of 
        resources data types, that are then passed along to calc_tc during 
        integration. 
    """

    def __init__(self,
                 gen_tc_function_obj=Gen.gen_tc_CM_LOOPS,
                 calc_tc_function_obj=Calc.calc_tc_CM_LOOPS):
        """ it has both generation and calculation functions assigned 
            to it """
        self.gen_tc_function_obj = gen_tc_function_obj
        self.calc_tc_function_obj = calc_tc_function_obj

    def gen_tc_resources(self, *args):
        """ *args: the generating functions usually only 
            need the pivots[] and i_max. resource_tuple is passed to 
            the specified calc function in calc_tc_vec """
        self.resource_tuple = self.gen_tc_function_obj(*args)

    def calc_tc_vec(self, *args):
        """ calculates in vector form all contributions to dNidt from 
            termination by combination. *args are usually only 
            Ndot_vec, i_max. The resource_tuple is passed before 
            the additional *args """
        return self.calc_tc_function_obj(*(self.resource_tuple), *args)

    def get_str(self):
        """ get state information as a string """
        return str(
            "(" + self.gen_tc_function_obj.__name__ + "," +
            self.calc_tc_function_obj.__name__ + ")")


class Slice_manager:
    """ In order to integrate other quantities than only dNidt, they 
        need to be packed all into the same vector [N, Ndot, mu, ...]. 
        If one wants to make experimental changes, one must keep track 
        of what is stored where in this vector. This class tries to 
        simplify that by storing starting index and size of 
        *sub*-vectors (slices) as they are inserted. """

    def __init__(self):
        """ A triplet consists of 
            (start index, size of slice, slice_initial_values). """
        self.list_triplets = []

    def push_slice(self, np_array_init_val=np.array([])):
        """ insert a vector (becomes a slice of the whole vector) """
        # get top index from slices
        size_sum = 0
        for start_idx, size, init_val in self.list_triplets:
            size_sum += size

        self.list_triplets.append(  # size_sum is top index
            (size_sum, np_array_init_val.size, np_array_init_val))

    def get_all_slices(self, np_array):
        """ returns slices in the order they were inserted as a list 
            of vectors, given the vector to be sliced as a parameter
            e.g.: N_vec, Ndot_vec, ... = sm.get_all_slices(whole_vec)
        """
        list_of_slices = []
        for s_idx, size, slice_init_vals in self.list_triplets:
            c_slice = np_array[s_idx:s_idx + size]
            if c_slice.size == 1:
                c_slice = np.asscalar(c_slice)
            list_of_slices.append(c_slice)

        return list_of_slices

    def get_initial_value_vec(self):
        """ returns the whole initial value vector """
        i_vec = np.array([], dtype=np.float64)
        for s_idx, size, slice_init_vals in self.list_triplets:
            i_vec = np.append(i_vec, slice_init_vals)

        return i_vec


class Data_unit:
    """ collects data and labels during integration of fPBEs for easier plotting and exporting """

    def __init__(self, x_vector, y_vector, x_name="x_name", y_name="y_name"):
        self.x_vector = x_vector
        self.y_vector = y_vector
        self.x_name = x_name
        self.y_name = y_name

    def append_xy_values(self, x_values=np.array([]), y_values=np.array([])):
        self.x_vector = np.append(self.x_vector, x_values)
        self.y_vector = np.append(self.y_vector, y_values)

    def get_xy_names(self):
        return self.x_name, self.y_name

    def get_xy_vectors(self):
        return self.x_vector, self.y_vector
