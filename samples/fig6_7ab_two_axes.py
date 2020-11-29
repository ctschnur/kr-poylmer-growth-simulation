from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# -- import most important classes
from kr.kr import Kr, Kr_constants, Tc_method_handler
# -- import functions that run fast, wich are later passed to a Kr object
from kr.gen import Gen
from kr.calc import Calc
from kr.plot_utils import Plot_utility, Snu


def mpl_settings():
    # -- plotting settings
    pgf_with_custom_preamble = {
        "figure.figsize": (6.4*0.8, 4.8*0.7),
        "savefig.format": "png",  # change this to pgf or png
        "pgf.rcfonts": False,  # don't use pre-specified fonts (in rcParmams)
        "text.usetex": True,  # use latex backend not just built-in MathText
        "text.latex.unicode": True,  # to be able to pass unicode characters to mpl
        "text.latex.preamble": [  # required for actual rendering to png
            r"\usepackage{amsmath, siunitx}",
        ],
        "pgf.preamble": [  # when exporting pgf code, mpl checks it for compilablity
            r"\usepackage{amsmath, siunitx}",
        ]}
    mpl.rcParams.update(pgf_with_custom_preamble)

    from cycler import cycler
    mpl.rcParams['axes.prop_cycle'] = cycler(color='bgrcmyk')


def compare_rfp_cld_fig6(my_kr):

    def plot_together(CLDs, refdatas, kr_obj, mpl_axes, r_fp_axes=None, plot_m_k_fp=False,
                      labels={"config_and_bundle_str": "",
                              "ref": "", "own": ""},
                      labelquantity2=None, linestyle="-", custom_xlim_2=None, ylabel="$m P_m \quad / (\si{.mol.L^{-1}})$", plot_r_fp_dead_dist=False, plot_r_fp_dead_dist_scale_factor=2e4):

        mpl_axes.set_xlabel("$\log_{10}m$")
        mpl_axes.set_ylabel(ylabel)

        from itertools import cycle
        marker_list = [".", "+", "x", "*"]
        marker_cycle = cycle(marker_list)

        for cur_refcld in refdatas:
            xcoords = cur_refcld["data"][:, 0]
            ycoords = cur_refcld["data"][:, 1]

            label_ref = ""
            if labels["ref"] is not "":
                label_ref += labels["ref"] + ", "

            mpl_axes.plot(xcoords, ycoords, next(marker_cycle), markersize=6, alpha=0.8, color="k",
                          label=label_ref + cur_refcld["label"])


        if plot_m_k_fp is True:
            m_vec = Plot_utility.fix_log_plot_args(CLDs[0]["pivots"])
            xcoords = np.log10(m_vec)
            ycoords = m_vec * kr_obj.d.k_fp
            r_fp_axes.plot(
                xcoords, ycoords,
                "-", label="$m \cdot k_{fp}$",
                linestyle="--")


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
            mpl_axes.plot(
                xcoords, nP_n, "-", label=label,
                linestyle=linestyle)

            if custom_xlim_2 is not None:
                mpl_axes.set_xlim(0, custom_xlim_2)
            else:
                mpl_axes.set_xlim(0., max(xcoords))

            if ((plot_r_fp_dead_dist is not False) and
                    (r_fp_axes is not None)):
                # dist_unscaled = cur_cld["pivots"] * cur_cld["r_fp_dead_dist"]
                dist_unscaled = cur_cld["r_fp_dead_dist"] / cur_cld["P_n"]
                plot_r_fp_dead_dist_scale_factor = 1.
                ycoords = dist_unscaled * plot_r_fp_dead_dist_scale_factor
                r_fp_axes.plot(
                    xcoords, ycoords,
                    "-", label=label,
                    linestyle="--")




    config_title = "FIG06"
    print(" -------- Comparing different imax for --------- ")
    print(" --- ", config_title)
    importdict = Kr_constants.import_dict_from_file("conf.json", config_title)
    importdict["conversions_to_save_clds_at"] = [0.4, 0.8, 0.81]
    # importdict["conversions_to_save_clds_at"].append(0.81)

    # 1. data acquisition
    # 2. plotting that data
    # --- clds
    refdatas_cld = [
        {
            "data": np.loadtxt(
                "paper_ref_data/butteFig6-KR-solid-40PercentConversion.csv",
                skiprows=1, delimiter=','),
            "label": r"Butté, $X = 0.4 $"  # , $ i_{\mathrm{max}} = 400$"
        },
        {
            "data": np.loadtxt(
                "paper_ref_data/butteFig6-KR-solid-80PercentConversion.csv",
                skiprows=1, delimiter=','),
            "label": r"Butté, $X = 0.8 $"  # , $ i_{\mathrm{max}} = 400$"
        }
    ]

    imaxlist = [400]
    linestyleslist = ["-", ":", "-."]

    evolutions_list = []
    clds_list = []
    for c_imax in imaxlist:
        print("c_imax: ", c_imax)
        importdict["imax"] = c_imax

        # show first, which simulation parameters are different
        # compared with DEFAULT configuration
        print("custom simulation parameters: \n",
              Kr_constants.compare_dicts(
                  importdict,
                  Kr_constants.import_dict_from_file("conf.json", "DEFAULT")))

        evolution_data_unit_dict, clds = my_kr.run(
            importdict,
            tc_method_handler=Tc_method_handler(Gen.gen_tc_CM_LOOPS,
                                                Calc.calc_tc_CM_LOOPS))
        evolutions_list.append(evolution_data_unit_dict)
        clds_list.append(clds)

    # 2. plotting
    # 2.1. plot clds (reference data + 3 with varying imax ) in
    # into one standalone figure

    fig_clds, ax1 = plt.subplots()
    default_size = fig_clds.get_size_inches()
    fig_clds.set_size_inches(default_size[0] * 1.25, default_size[1] * 1.25)
    # plot reference clds

    imax_snu = Snu("imax", importdict["imax"])
    imax_snu.set_sym_latex(r"i_{\mathrm{max}}")

    plot_together(
        CLDs=[], refdatas=refdatas_cld,
        labels={
            "config_and_bundle_str": "reference_cld",
            "ref": "",
            "own": ""
        },
        kr_obj=my_kr, mpl_axes=ax1,
        labelquantity2=imax_snu.get_latex_senu(dollars=True))

    ax2 = ax1.twinx()
    ax2.set_ylabel(
        r"$ m \cdot k_{fp} \quad / (\si{.L.mol^{-1}.s^{-1}})$")

    for index, c_clds in enumerate(clds_list):
        imax_snu.set_num(imaxlist[index])
        first2_clds = c_clds[:-1]
        plot_together(
            CLDs=first2_clds, refdatas=[],
            labels={
                "config_and_bundle_str": "Fig06-own-clds-increasing-imax",
                "ref": "reference simulated curve",
                "own": ""
            },
            kr_obj=my_kr, mpl_axes=ax1, 
            r_fp_axes=ax2,
            plot_m_k_fp=True,
            linestyle=linestyleslist[index])


    ax1.grid(True)
    ax2.grid(False)

    cur_ylim = ax1.get_ylim()
    ax1.set_ylim(0, cur_ylim[1])
    cur_ylim = ax2.get_ylim()
    ax2.set_ylim(0, cur_ylim[1])

    ax1.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    ax2.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    ax1.legend(loc='upper left', frameon=True, title="$m P_m$ - Achse:")
    ax2.legend(loc='upper right', bbox_to_anchor=(1.-0.1, 1),
               frameon=True, title="$ m \cdot k_{fp}$ - Achse: ")
    plt.tight_layout()
    plt.savefig(my_kr.get_run_hash_str() + "_" + "plot_together")
    plt.savefig(my_kr.get_run_hash_str() + "_" +
                "plot_together" + ".pgf", format="pgf")

    # plot rate of chain transfer to polymer at specific conversions
    # for different terminated chain lenghts
    # we have m * P_m
    # we need sum_n Rdot_n  (= lamda_0)
    # r_fp = k_fp * m * P_m * Rdot_n
    # sum_n=0^size(P_vec) r_fp = k_fp * m * P_m * Rdot_n
def main():
    mpl_settings()  # matplotlib settings

    my_kr = Kr()

    compare_rfp_cld_fig6(my_kr)


if __name__ == "__main__":
    main()


