from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import numpy as np
# import matplotlib.pyplot as plt
import matplotlib as mpl

# -- import most important classes
from kr.kr import Kr, Kr_constants, Tc_method_handler
# -- import functions that run fast, wich are later passed to a Kr object
from kr.gen import Gen
from kr.calc import Calc
from kr.plot_utils import Plot_utility


def fig9(my_kr):
    config_title = "FIG09"
    config_and_bundle_str = "FIG09" + "(" + config_title + ")"
    print(" --- ", config_title)
    importdict = Kr_constants.import_dict_from_file("conf.json", config_title)

    # show first, which simulation parameters are different
    # compared with DEFAULT configuration
    print("custom simulation parameters: \n",
          Kr_constants.compare_dicts(
              importdict,
              Kr_constants.import_dict_from_file("conf.json", "DEFAULT")))

    evolution_data_unit_dict, clds = my_kr.run(
        importdict,
        tc_method_handler=Tc_method_handler(
            Gen.gen_tc_CM_LOOPS,
            Calc.calc_tc_CM_LOOPS))

    # --- clds
    bundle_name = "clds"
    config_and_bundle_str += "_" + bundle_name
    refdatas_cld = [
        {
            "data": np.loadtxt(
                "paper_ref_data/butteFig9-KR-solid-40percentConversion.csv",
                skiprows=1, delimiter=','),
            "label": "$X=0.4$"
        },
        {
            "data": np.loadtxt(
                "paper_ref_data/butteFig9-KR-solid-80percentConversion.csv",
                skiprows=1, delimiter=','),
            "label": "$X=0.8$"
        }
    ]

    Plot_utility.plot_clds(
        clds, refdatas=refdatas_cld,
        labels={
            "config_and_bundle_str": config_and_bundle_str,
            "ref": "Butté et al.",
            "own": "eigene Daten"
        },
        kr_obj=my_kr)


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

    fig9(my_kr)


if __name__ == "__main__":
    main()
