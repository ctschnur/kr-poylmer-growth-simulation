from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# -- import most important classes
from kr.kr import Kr, Kr_constants
# -- import functions that run fast, wich are later passed to a Kr object
from kr.plot_utils import Plot_utility


def fig4(my_kr):
    config_title = "FIG04"
    print(" --- ", config_title)
    # import a configuration of initial parameters
    importdict = Kr_constants.import_dict_from_file("conf.json", config_title)

    # show first, which simulation parameters are different
    # compared with DEFAULT configuration
    print("custom simulation parameters: \n",
          Kr_constants.compare_dicts(
              importdict,
              Kr_constants.import_dict_from_file("conf.json", "DEFAULT")))

    evolution_data_unit_dict, clds = my_kr.run(importdict)

    # --- dispersity, until X = 0.8
    bundle_name = "dispersity"
    config_and_bundle_str = config_title + "_" + bundle_name
    x_name, y_name = evolution_data_unit_dict[bundle_name].get_xy_names()

    fig, ax = plt.subplots()
    ax.set_xlabel("$X$")
    ax.set_ylabel("$\mathrm{PDI}$")

    # reference data
    x_ref, y_ref = np.loadtxt(
        "paper_ref_data/butteFig4b-solid-01moments.csv",
        skiprows=1, delimiter=',', unpack=True)

    ax.plot(x_ref, y_ref, ".", markersize=5, alpha=0.8, color="k",
            label="reference simulated curve")

    # bundle data
    x, y = evolution_data_unit_dict[bundle_name].get_xy_vectors()
    ax.plot(x, y, linestyle='-', label="simulated")

    ax.grid()
    ax.legend(loc="best")
    plt.tight_layout()
    plt.savefig(my_kr.get_run_hash_str() + "_" +
                config_and_bundle_str)
    plt.savefig(my_kr.get_run_hash_str() + "_" +
                config_and_bundle_str + ".pgf", format="pgf")

    fig_clds = plt.figure()
    # show only cld for X = 0.6
    clds = [clds[0]]
    Plot_utility.plot_clds(
        clds,
        refdatas=[
            {
                "data": np.loadtxt(
                    "paper_ref_data/butteFig4a.csv", skiprows=1, delimiter=','),
                "label": "$X=0.6$"
            }
        ],
        labels={"config_and_bundle_str": config_title + "_clds",
                "ref": "reference simulated curve",
                "own": "simulated"},
        kr_obj=my_kr, mpl_figure=fig_clds, savefig=True, savefig_pgf=True)


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

    fig4(my_kr)


if __name__ == "__main__":
    main()
