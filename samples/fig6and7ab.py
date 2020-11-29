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
from kr.plot_utils import Plot_utility


def fig6_7ab(my_kr):
    config_title = "FIG06"
    print(" --- ", config_title)
    # import a configuration of initial parameters
    importdict = Kr_constants.import_dict_from_file("conf.json", config_title)
    importdict["conversions_to_save_clds_at"] = [0.4, 0.8, 0.99]

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

    # --- clds
    bundle_name = "clds"
    config_and_bundle_str = config_title + bundle_name
    refdatas_cld = [
        {
            "data": np.loadtxt(
                "paper_ref_data/butteFig6-KR-solid-40PercentConversion.csv",
                skiprows=1, delimiter=','),
            "label": "$X = 0.4$"
        },
        {
            "data": np.loadtxt(
                "paper_ref_data/butteFig6-KR-solid-80PercentConversion.csv",
                skiprows=1, delimiter=','),
            "label": "$X = 0.8$"
        }
    ]

    Plot_utility.plot_clds(
        clds, refdatas=refdatas_cld,
        labels={
            "config_and_bundle_str": config_and_bundle_str,
            "ref": "reference simulated curve",
            "own": "own simulated curve"
        },
        kr_obj=my_kr, savefig=True, savefig_pgf=True)

    plt.clf()
    from matplotlib import gridspec
    fig = plt.figure()
    xd, yd = fig.get_size_inches()  # get default size
    fig.set_size_inches([xd * 1., yd * 2. * 0.75])
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
             label="reference simulated curve")

    x, y = evolution_data_unit_dict[bundle_name].get_xy_vectors()
    ax1.plot(x, y, linestyle='-', label="own simulated curve")

    # --- gelfraction
    bundle_name = "gelfraction"
    x_name, y_name = evolution_data_unit_dict[bundle_name].get_xy_names()
    ax2.set_xlabel("$X$")
    ax2.set_ylabel("$g$")

    config_and_bundle_str += bundle_name

    # reference data
    ref_x, ref_y = np.loadtxt("paper_ref_data/butteFig7b-KR-solid.csv",
                              skiprows=1, delimiter=',', unpack=True)
    plt.plot(ref_x, ref_y, ".", markersize=5, alpha=0.8, color="k",
             label="reference simulated curve")

    # bundle data
    x, y = evolution_data_unit_dict[bundle_name].get_xy_vectors()
    plt.plot(x, y, linestyle='-', label="own simulated curve")

    ax1.grid()
    ax2.grid()
    ax1.legend(loc='best')
    ax2.legend(loc='best')
    plt.tight_layout()
    plt.savefig(my_kr.get_run_hash_str() + config_and_bundle_str)
    plt.savefig(my_kr.get_run_hash_str() +
                config_and_bundle_str + ".pgf", format="pgf")
    plt.gcf().clear()


# --- general functions needed for plotting
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


def main():
    mpl_settings()  # matplotlib settings

    my_kr = Kr()

    fig6_7ab(my_kr)


if __name__ == "__main__":
    main()
