from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# -- import most important classes
from kr.kr import Kr, Kr_constants

from kr.plot_utils import Plot_utility


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


def fig5(my_kr):
    config_title = "FIG05-1000"
    print(" --- ", config_title)
    # import a configuration of initial parameters
    importdict = Kr_constants.import_dict_from_file("conf.json", config_title)

    # show first, which simulation parameters are different
    # compared with DEFAULT configuration
    print("custom simulation parameters: \n",
          Kr_constants.compare_dicts(
              importdict,
              Kr_constants.import_dict_from_file("conf.json", "DEFAULT")))

    evolution_data_unit_dict, clds = my_kr.run(
        importdict)

    # --- dispersity
    bundle_name = "dispersity"
    config_and_bundle_str = config_title + bundle_name
    x_name, y_name = evolution_data_unit_dict[bundle_name].get_xy_names()

    fig, ax = plt.subplots()
    ax.set_xlabel("$X$")
    ax.set_ylabel("$\mathrm{PDI}$")

    # reference data
    ref_x, ref_y = np.loadtxt("paper_ref_data/butteFig5b-KR1000pointssolid.csv",
                              skiprows=1, delimiter=',', unpack=True)
    ax.plot(ref_x, ref_y, ".", markersize=5, alpha=0.8, color="k",
            label=r"Butté et al.")

    # bundle data
    x, y = evolution_data_unit_dict[bundle_name].get_xy_vectors()
    ax.plot(x, y, linestyle='-', label="eigene Daten")

    ax.grid()
    ax.legend(loc='upper left')
    plt.tight_layout()
    plt.savefig(my_kr.get_run_hash_str() + config_and_bundle_str)
    plt.savefig(my_kr.get_run_hash_str() + "_" +
                config_and_bundle_str + ".pgf", format="pgf")

    fig_clds = plt.figure()
    Plot_utility.plot_clds(
        clds,
        refdatas=[
            {
                "data": np.loadtxt(
                    "paper_ref_data/butteFig5a-KR1000pointssolid.csv",
                    skiprows=1, delimiter=','),
                "label": "$X=0.6$"
            }
        ],
        labels={"config_and_bundle_str": "FIG05-1000",
                "ref": r"Butté et al.", "own": "eigene Daten"},
        kr_obj=my_kr, savefig=True, savefig_pgf=True, mpl_figure=fig_clds)


def main():
    mpl_settings()  # matplotlib settings

    my_kr = Kr()

    fig5(my_kr)


if __name__ == "__main__":
    main()
