"""
helps plotting clds
and misc other quantities over conversion
"""

import numpy as np
import matplotlib.pyplot as plt


class Plot_utility:
    @staticmethod
    def plot_clds(CLDs, refdatas, kr_obj,
                  labels={"config_and_bundle_str": "", "ref": "", "own": ""},
                  export=False, showtitle=False, mpl_figure=plt.figure(), savefig=True, savefig_pgf=False, labelquantity2=None, linestyle="-", custom_xlim_2=None, ylabel="$n P_n \quad / (\si{.mol.L^{-1}})$", plot_r_fp_dead_dist=False):
        """ Parameters: CLDs: list of dictionaries each containing one bundled-up
                            cld with additional printing information.
                        refdatas: dictionaries with reference clds
                        kr_obj: for the filename (hash) and axvline markings
                                indicating certain pivots
                        export: string; if export is not None, it 
                                     will write n and P_n to the specified path
                                     with the corresponding conversion and hash
                                     attached to the filename
                        plot_r_fp_dead_dist_scale_factor: found out through numerical experiment
                                                          for FIG06
        """

        # activate figure
        plt.figure(mpl_figure.number)

        config_and_bundle_str = labels["config_and_bundle_str"]

        if showtitle is True:
            plt.title(config_and_bundle_str)

        plt.xlabel("$\log_{10}n$")
        plt.ylabel(ylabel)

        from itertools import cycle
        marker_list = [".", "+", "x", "*"]
        marker_cycle = cycle(marker_list)

        for cur_refcld in refdatas:
            xcoords = cur_refcld["data"][:, 0]
            ycoords = cur_refcld["data"][:, 1]

            label_ref = ""
            if labels["ref"] is not "":
                label_ref += labels["ref"] + ", "

            plt.plot(xcoords, ycoords, next(marker_cycle), markersize=5, alpha=0.8, color="k",
                     label=label_ref + cur_refcld["label"])

        max_x_coord = 0.
        for cur_cld in CLDs:
            xcoords = np.log10(
                Plot_utility.fix_log_plot_args(cur_cld["pivots"]))
            max_x_coord = max(max(xcoords), max_x_coord)
            label = ""
            if labels["own"] is not "":
                label += labels["own"] + ", "

            label += ("$X = " + str(cur_cld["conversion"]) + " $")
            if labelquantity2 is not None:
                label += (", " + labelquantity2)

            nP_n = cur_cld["pivots"] * cur_cld["P_n"]
            base_line, = plt.plot(
                xcoords, nP_n, "-", label=label,
                linestyle=linestyle)

            if export is not False:
                Plot_utility.export_xy(
                    kr_obj.get_run_hash_str() + config_and_bundle_str + "_" +
                    str(cur_cld["conversion"]) + ".csv",
                    cur_cld["pivots"], cur_cld["P_n"])

            if custom_xlim_2 is not None:
                plt.xlim(0, custom_xlim_2)
            else:
                plt.xlim(0., max(xcoords))

            if plot_r_fp_dead_dist is not False:
                # scale r_fp_dead_dist, so that max(r_fp_dead_dist) is max(m
                # P_m)
                dist_unscaled = cur_cld["pivots"] * cur_cld["r_fp_dead_dist"]
                # ycoords = dist_unscaled * (np.max(nP_n) /
                # np.max(dist_unscaled)) # this scaling would be non-uniform if
                # several distributions (at different conversions) are plotted together

                ycoords = dist_unscaled

                print("(np.max(nP_n) / np.max(dist_unscaled)): ",
                      (np.max(nP_n) / np.max(dist_unscaled)),
                      ", np.max(dist_unscaled): ",
                      np.max(dist_unscaled))
                plt.plot(
                    xcoords, ycoords,
                    "-", label=label + r"$\sum_{m=0}^{\infty} r_{fp}$",
                    linestyle=linestyle)

        # plt.axvline(np.log10(kr_obj.pivots[kr_obj.imax]),
        #             label="nbar_imax",
        #             linestyle="-", alpha=0.5,
        #             color=base_line.get_color())

        plt.grid(True)
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        plt.legend(loc='best', frameon=True)
        if savefig is True:
            plt.tight_layout()
            plt.savefig(kr_obj.get_run_hash_str() +
                        "_" + config_and_bundle_str)
        if savefig_pgf is True:
            plt.tight_layout()
            plt.savefig(kr_obj.get_run_hash_str() + "_" +
                        config_and_bundle_str + ".pgf", format="pgf")

    @staticmethod
    def fix_log_plot_args(vec):
        """ log of numbers <= 0 is undefined, so set 0 to 0.1, 
            which prevents numpy from complaining """
        vec_fixed = np.array(vec.copy()).astype(np.float64)
        vec_fixed[0] = 0.1
        return vec_fixed

    @staticmethod
    def export_xy(filename, x_values, y_values):
        np.savetxt(filename, (x_values, y_values))


class Snu:
    def __init__(self, sym_code, num_code, stay_code=False):
        self.sym_code = sym_code
        self.sym_latex = str(sym_code)

        self.num_code = num_code
        self.num_latex = str(num_code)

        self.unit = ""

        self.stay_code = stay_code

        self.format_num_latex("{0}")  # automatically format if it's an array

    def set_sym_latex_stay_code(self, stay_code=True):
        if stay_code is True:
            self.sym_latex = "\\texttt{{{0}}}".format(
                self.sym_code).replace("_", "\_")
        return self

    def set_sym_latex(self, sym_latex):
        if sym_latex is not None:
            self.sym_latex = sym_latex
        return self

    def format_num_latex(self, fstr):
        self.fstr = fstr
        if fstr is not None:
            self.num_latex = fstr.format(self.num_code)
        if isinstance(self.num_code, list):
            self.num_latex = "\\text{"
            for nb in self.num_code[:-1]:  # format array of numbers
                self.num_latex += (fstr + ", ").format(nb)
            self.num_latex += fstr.format(self.num_code[-1])
            self.num_latex += "}"

        return self

    def set_unit(self, unit=None):
        if unit is not None:
            self.unit = unit
        return self

    def set_num(self, num, reformat=True):
        self.num_code = num
        self.num_latex = num
        if reformat is True:
            self.format_num_latex(self.fstr)
        return self

    def get_latex_senu(self, dollars=False):  # returns string "sym = num unit"
        s = ""
        s += self.sym_latex
        s += " = "
        s += self.num_latex
        s += self.unit
        if dollars is True:
            s = "$ " + s + " $"
        return s

    @staticmethod
    def get_snus_dict_from_normal_dict(normal_dict):
        snus_dict = {}
        for key, value in normal_dict.items():
            snus_dict[key] = Snu(key, value)

        return snus_dict

    @staticmethod
    def get_normal_dict_from_snus_dict(snus_dict):
        """ snus_dict is a dictionary 
            with (key = the python variable's name, 
                  value = snu(python variable's name, value) 
            the returned normal dict contains only 
            strings for keys and values """
        normal_dict = {}
        for key, value in snus_dict.items():
            normal_dict[value.sym_latex] = r"{0:<10} & {1:<10}".format(
                value.num_latex, value.unit)

        return normal_dict

    @staticmethod
    def latex_table_from_dictionary(sym_val_dict):
        """ this function puts the pre-formatted dict into a latex table, sorted alphabetically"""

        table_str = ""
        # sym_val_dict_copy = dict(sym_val_dict)
        import collections
        sym_val_dict_copy = collections.OrderedDict(
            sorted(sym_val_dict.items()))
        for key, value in sym_val_dict_copy.items():
            table_str += " {:<10} & ".format(str(key))
            table_str += " {:<10}".format(str(value))
            table_str += "\\\\\n"
        return table_str
