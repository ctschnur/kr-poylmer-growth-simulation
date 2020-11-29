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


def compare_imax_fig6(my_kr):
    config_title = "FIG06"
    print(" -------- Comparing different imax for --------- ")
    print(" --- ", config_title)
    importdict = Kr_constants.import_dict_from_file("conf.json", config_title)
    # importdict["conversions_to_save_clds_at"] = [0.2, 0.4, 0.4]
    importdict["conversions_to_save_clds_at"].append(0.99)

    # 1. data acquisition
    # 2. plotting that data
    # --- clds
    refdatas_cld = [
        {
            "data": np.loadtxt(
                "paper_ref_data/butteFig6-KR-solid-40PercentConversion.csv",
                skiprows=1, delimiter=','),
            "label": r"$X = 0.4 $, $ i_{\mathrm{max}} = 400$"
        },
        {
            "data": np.loadtxt(
                "paper_ref_data/butteFig6-KR-solid-80PercentConversion.csv",
                skiprows=1, delimiter=','),
            "label": r"$X = 0.8 $, $ i_{\mathrm{max}} = 400$"
        }
    ]

    imaxlist = [400, 100, 1000]
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

    fig_clds = plt.figure()
    default_size = fig_clds.get_size_inches()
    fig_clds.set_size_inches(default_size[0] * 1.25, default_size[1] * 1.25)
    # plot reference clds

    imax_snu = Snu("imax", importdict["imax"])
    imax_snu.set_sym_latex(r"i_{\mathrm{max}}")

    Plot_utility.plot_clds(
        CLDs=[], refdatas=refdatas_cld,
        labels={
            "config_and_bundle_str": "reference_cld",
            "ref": "Butté",
            "own": ""
        },
        kr_obj=my_kr, savefig=False, mpl_figure=fig_clds,
        labelquantity2=imax_snu.get_latex_senu(dollars=True))

    for index, c_clds in enumerate(clds_list):
        imax_snu.set_num(imaxlist[index])
        first2_clds = c_clds[:-1]
        Plot_utility.plot_clds(
            CLDs=first2_clds, refdatas=[],
            labels={
                "config_and_bundle_str": "Fig06-own-clds-increasing-imax",
                "ref": "reference simulated curve",
                "own": ""
            },
            kr_obj=my_kr, savefig_pgf=True, mpl_figure=fig_clds,
            labelquantity2=imax_snu.get_latex_senu(dollars=True),
            linestyle=linestyleslist[index])

    # exit()

    # --- evolutions

    from matplotlib import gridspec
    fig_evols = plt.figure()
    xd, yd = fig_evols.get_size_inches()  # get default size
    fig_evols.set_size_inches([xd * 1.25 * 1., yd * 1.35 * 2. * 0.75])
    gs = gridspec.GridSpec(2, 1, height_ratios=[0.5, 1])
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])

    ax1.text(-0.1, 1.0, r"\textbf{a)}", transform=ax1.transAxes)
    ax2.text(-0.1, 1.0, r"\textbf{b)}", transform=ax2.transAxes)

    # --- dispersity
    bundle_name = "dispersity"
    x_name, y_name = evolution_data_unit_dict[bundle_name].get_xy_names()
    ax1.set_xlabel("$X$")
    ax1.set_ylabel("$\mathrm{PDI}$")
    ax1.set_ylim(0., 35)

    config_and_bundle_str = "FIG07ab" + "(" + config_title + ")" + bundle_name

    # reference data
    ref_x, ref_y = np.loadtxt("paper_ref_data/butteFig7a-KR-solid.csv",
                              skiprows=1, delimiter=',', unpack=True)
    ax1.plot(ref_x, ref_y, ".", markersize=5, alpha=0.8, color="k",
             label=("reference simulated curve" + ", " +
                    imax_snu.set_num(400).get_latex_senu(dollars=True)))

    for index, c_evolution in enumerate(evolutions_list):
        imax_snu.set_num(imaxlist[index])
        x, y = c_evolution[bundle_name].get_xy_vectors()
        ax1.plot(x, y, linestyle=linestyleslist[index],
                 label=imax_snu.get_latex_senu(dollars=True))

    # --- gelfraction
    bundle_name = "gelfraction"
    x_name, y_name = evolution_data_unit_dict[bundle_name].get_xy_names()
    ax2.set_xlabel("$X$")
    ax2.set_ylabel("$g$")

    config_and_bundle_str += bundle_name

    # reference data
    ref_x, ref_y = np.loadtxt("paper_ref_data/butteFig7b-KR-solid.csv",
                              skiprows=1, delimiter=',', unpack=True)
    ax2.plot(ref_x, ref_y, ".", markersize=5, alpha=0.8, color="k",
             label=("reference simulated curve" + ", " +
                    imax_snu.set_num(400).get_latex_senu(dollars=True)))

    for index, c_evolution in enumerate(evolutions_list):
        imax_snu.set_num(imaxlist[index])
        x, y = c_evolution[bundle_name].get_xy_vectors()
        ax2.plot(x, y, linestyle=linestyleslist[index],
                 label=imax_snu.get_latex_senu(dollars=True))

    ax1.grid()
    ax2.grid()
    ax1.legend(loc='best')
    ax2.legend(loc='best')
    plt.tight_layout()
    plt.savefig(my_kr.get_run_hash_str() + "_" + config_and_bundle_str)
    plt.savefig((my_kr.get_run_hash_str() +
                 "_" + config_and_bundle_str + ".pgf"), format="pgf")

    # exit()

    # --- generate plots of details with different pivot generation method
    # use all the curves above, in addition to the same curves

    importdict["conversions_to_save_clds_at"][-1] = 0.83

    evolutions_list_m2 = []
    clds_list_m2 = []
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
            gen_pivots_func=Gen.gen_pivots_method_2,
            tc_method_handler=Tc_method_handler(Gen.gen_tc_CM_LOOPS,
                                                Calc.calc_tc_CM_LOOPS))
        evolutions_list_m2.append(evolution_data_unit_dict)
        clds_list_m2.append(clds)

    # --- evolutions

    from matplotlib import gridspec
    fig_evols_2 = plt.figure()
    xd, yd = fig_evols_2.get_size_inches()  # get default size
    fig_evols_2.set_size_inches([xd * 1.25 * 1., yd * 1.35 * 2. * 0.75])
    gs = gridspec.GridSpec(3, 1, height_ratios=[0.5, 0.5, 0.5])
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    ax3 = plt.subplot(gs[2])

    ax1.text(-0.1, 1.0, r"\textbf{a)}", transform=ax1.transAxes)
    ax2.text(-0.1, 1.0, r"\textbf{b)}", transform=ax2.transAxes)
    ax3.text(-0.1, 1.0, r"\textbf{c)}", transform=ax3.transAxes)

    # --- dispersity
    bundle_name = "dispersity"
    x_name, y_name = evolution_data_unit_dict[bundle_name].get_xy_names()
    ax1.set_xlabel("$X$")
    ax1.set_ylabel("$\mathrm{PDI}$")

    ax1.set_xlim(0.37, 0.83)
    ax1.set_ylim(15., 29.)

    config_and_bundle_str = "FIG07ab" + "(" + config_title + ")" + bundle_name

    # reference data
    ref_x, ref_y = np.loadtxt("paper_ref_data/butteFig7a-KR-solid.csv",
                              skiprows=1, delimiter=',', unpack=True)
    ax1.plot(ref_x, ref_y, ".", markersize=5, alpha=0.8, color="k",
             label=("reference simulated curve" + ", " +
                    imax_snu.set_num(400).get_latex_senu(dollars=True)))

    for index, c_evolution in enumerate(evolutions_list[:]):
        imax_snu.set_num(imaxlist[index])
        x, y = c_evolution[bundle_name].get_xy_vectors()
        ax1.plot(x, y, linestyle=linestyleslist[index],
                 label=imax_snu.get_latex_senu(dollars=True) + ", PM 1")
        # for each imax, plot the other method
        x, y = evolutions_list_m2[index][bundle_name].get_xy_vectors()
        ax1.plot(x, y, linestyle=linestyleslist[index],
                 label=imax_snu.get_latex_senu(dollars=True) + ", PM 2")

    # --- gelfraction
    bundle_name = "gelfraction"
    x_name, y_name = evolution_data_unit_dict[bundle_name].get_xy_names()
    ax2.set_xlabel("$X$")
    ax2.set_ylabel("$g$")

    ax2.set_ylim(0., 0.19)
    ax2.set_xlim(0.25, 0.77)

    config_and_bundle_str += bundle_name

    # reference data
    ref_x, ref_y = np.loadtxt("paper_ref_data/butteFig7b-KR-solid.csv",
                              skiprows=1, delimiter=',', unpack=True)
    ax2.plot(ref_x, ref_y, ".", markersize=5, alpha=0.8, color="k",
             label=("reference simulated curve" + ", " +
                    imax_snu.set_num(400).get_latex_senu(dollars=True)))

    for index, c_evolution in enumerate(evolutions_list[:]):
        imax_snu.set_num(imaxlist[index])
        x, y = c_evolution[bundle_name].get_xy_vectors()
        ax2.plot(x, y, linestyle=linestyleslist[index],
                 label=imax_snu.get_latex_senu(dollars=True) + ", PM 1")
        # for each imax, plot the other method
        x, y = evolutions_list_m2[index][bundle_name].get_xy_vectors()
        ax2.plot(x, y, linestyle=linestyleslist[index],
                 label=imax_snu.get_latex_senu(dollars=True) + ", PM 2")

    # clds reference
    ax3.set_xlabel("$\log_{10}n$")
    nPn_label = "$n P_n \quad / (\si{.mol.L^{-1}})$"
    ax3.set_ylabel(nPn_label)

    ax3.set_xlim(1.5-0.2, 3.5-0.2)
    ax3.set_ylim(2.75*1e-4, 3.35*1e-4)

    ref_x, ref_y = np.loadtxt("paper_ref_data/butteFig6-KR-solid-80PercentConversion.csv",
                              skiprows=1, delimiter=',', unpack=True)
    ax3.plot(ref_x, ref_y, ".", markersize=5, alpha=0.8, color="k",
             label=("reference simulated curve" + ", " + "$X = 0.8$, " +
                    imax_snu.set_num(400).get_latex_senu(dollars=True)))

    # clds_list (only 100 and 400)
    for index, c_clds in enumerate(clds_list[:]):
        imax_snu.set_num(imaxlist[index])
        # iterate over clds

        # Pivot Method 1
        for cur_cld in c_clds[1:2]:  # only X=0.4 and 0.8, not 0.9
            xcoords = np.log10(
                Plot_utility.fix_log_plot_args(cur_cld["pivots"]))

            label = "$X = " + str(cur_cld["conversion"]) + " $"
            label += ", " + \
                imax_snu.set_num(cur_cld["imax"]).get_latex_senu(dollars=True)

            ax3.plot(
                xcoords, cur_cld["pivots"] * cur_cld["P_n"], linestyle=linestyleslist[index],
                label=label + ", PM 1")

        # Pivot Method 2
        for cur_cld in clds_list_m2[index][1:2]:  # only 0.4 and 0.8
            xcoords = np.log10(
                Plot_utility.fix_log_plot_args(cur_cld["pivots"]))

            label = "$X = " + str(cur_cld["conversion"]) + " $"
            label += ", " + \
                imax_snu.set_num(cur_cld["imax"]).get_latex_senu(dollars=True)

            ax3.plot(
                xcoords, cur_cld["pivots"] * cur_cld["P_n"], linestyle=linestyleslist[index],
                label=label + ", PM 2")

    ax3.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

    # new data ---

    ax1.grid()
    ax2.grid()
    ax3.grid()
    ax1.legend(loc='best')
    ax2.legend(loc='best')
    ax3.legend(loc='best')
    plt.tight_layout()
    plt.savefig(my_kr.get_run_hash_str() + "_" + config_and_bundle_str + "_3")
    plt.savefig((my_kr.get_run_hash_str() +
                 "_" + config_and_bundle_str + "_3" + ".pgf"), format="pgf")


def mpl_settings():
    # -- plotting settings
    pgf_with_custom_preamble = {
        "figure.figsize": (6.4*0.8, 4.8*0.7),
        "legend.fontsize": 10,
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


def main():
    mpl_settings()  # matplotlib settings

    my_kr = Kr()

    compare_imax_fig6(my_kr)


if __name__ == "__main__":
    main()
