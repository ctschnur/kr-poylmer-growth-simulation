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


def fig8(my_kr):
    config_title = "FIG06"
    config_and_bundle_str = "FIG08" + "(" + config_title + ")"
    print(" --- ", config_and_bundle_str)
    # import a configuration of initial parameters
    importdict = Kr_constants.import_dict_from_file("conf.json", config_title)
    importdict["conversions_to_save_clds_at"] = [0.99]

    # show first, which simulation parameters are different
    # compared with DEFAULT configuration
    print("custom simulation parameters: \n",
          Kr_constants.compare_dicts(
              importdict,
              Kr_constants.import_dict_from_file("conf.json", "DEFAULT")))

    # importdict["imax"] = 800

    nbar_imax_list = np.array([2.e6, 5.e6, 1.e7])
    linestyle_list = ["-", "--", "-."]

    fig1, ax1 = plt.subplots(nrows=1, ncols=1)
    fig2, ax2 = plt.subplots(nrows=1, ncols=1)

    # --- dispersity
    ax1_bundle_name = "dispersity"
    ax1.set_xlabel("$X$")
    ax1.set_ylabel("$\mathrm{PDI}$")
    ax1.set_ylim(0., 70.)

    # --- gelfraction
    ax2.set_xlabel("$X$")
    ax2.set_ylabel("$g$")

    # reference data
    ref_x, ref_y = np.loadtxt("paper_ref_data/butteFig8a-DashDotted-nbarimax2e6.csv",
                              skiprows=1, delimiter=',', unpack=True)
    ax1.plot(ref_x, ref_y,
             label=r"Butté et al., $\bar{n}_{i_\mathrm{max}} = \num{2e6}$", linestyle='-', color='k')

    ref_x, ref_y = np.loadtxt("paper_ref_data/butteFig8a-DashDotted-nbarimax5e6.csv",
                              skiprows=1, delimiter=',', unpack=True)
    ax1.plot(ref_x, ref_y,
             label=r"Butté et al., $\bar{n}_{i_\mathrm{max}} = \num{5e6}$", linestyle='--', color='k')

    ref_x, ref_y = np.loadtxt("paper_ref_data/butteFig8a-DashDotted-nbarimax1e7.csv",
                              skiprows=1, delimiter=',', unpack=True)
    ax1.plot(ref_x, ref_y,
             label=r"Butté et al., $\bar{n}_{i_\mathrm{max}} = \num{1e7}$", linestyle='-.', color='k')

    # for gelfraction, collect line data from plots,
    # so a subplot containing a zoomed area an be created

    ref_x, ref_y = np.loadtxt("paper_ref_data/butteFig8b-DashDotted-nbarimax2e6.csv",
                              skiprows=1, delimiter=',', unpack=True)
    line1 = plt.Line2D(
        ref_x, ref_y, label=r"Butté et al., $\bar{n}_{i_\mathrm{max}} = \num{2e6}$", linestyle='-', color='k')

    ref_x, ref_y = np.loadtxt("paper_ref_data/butteFig8b-DashDotted-nbarimax5e6.csv",
                              skiprows=1, delimiter=',', unpack=True)
    line2 = plt.Line2D(
        ref_x, ref_y, label=r"Butté et al., $\bar{n}_{i_\mathrm{max}} = \num{5e6}$", linestyle='--', color='k')

    ref_x, ref_y = np.loadtxt("paper_ref_data/butteFig8b-DashDotted-nbarimax1e7.csv",
                              skiprows=1, delimiter=',', unpack=True)
    line3 = plt.Line2D(
        ref_x, ref_y, label=r"Butté et al., $\bar{n}_{i_\mathrm{max}} = \num{1e7}$", linestyle='-.', color='k')

    lines_gelfraction = [line1, line2, line3]
    for index, c_nbar_imax in enumerate(nbar_imax_list):
        importdict["nbar_imax"] = c_nbar_imax
        evolution_data_unit_dict, clds = my_kr.run(
            importdict,
            tc_method_handler=Tc_method_handler(
                Gen.gen_tc_CM_LOOPS,
                Calc.calc_tc_CM_LOOPS))

        # dispersity -- actually plot (creates Line2D and assigns to axes, fig.)
        x, y = evolution_data_unit_dict[ax1_bundle_name].get_xy_vectors()
        ax1.plot(x, y, linestyle=linestyle_list[index],
                 label=r"eigene Daten, $\bar{n}_{i_\mathrm{max}} = " + r"\num{{{0:.0E}}}$".format(c_nbar_imax), color="C1")

        # gel fraction -- don't assign Line2D with axes or figure, just create it
        x, y = evolution_data_unit_dict["gelfraction"].get_xy_vectors()
        line = plt.Line2D(x, y, linestyle=linestyle_list[index],
                          label=r"eigene Daten, $\bar{n}_{i_\mathrm{max}} = " + r"\num{{{0:.0E}}}$".format(c_nbar_imax))

        lines_gelfraction.append(line)

    # set colors for gelfraction
    lines_gelfraction[-1].set_color("C1")
    lines_gelfraction[-2].set_color("C1")
    lines_gelfraction[-3].set_color("C1")
    # now add Line2D's to ax2, which is essentially the same as plt.plot()
    import copy  # I'm not sure how much is actually copied, deepcopy didn't work
    [ax2.add_line(copy.copy(line)) for line in lines_gelfraction]

    ax1.grid()
    ax2.grid()
    ax1.legend(loc='best')
    ax2.legend(loc='best')
    ax2.set_ylim(0.0, 0.8)

    plt.figure(fig1.number)
    plt.tight_layout()
    plt.savefig(my_kr.get_run_hash_str() +
                config_and_bundle_str + "_dispersities")
    plt.savefig(my_kr.get_run_hash_str() + config_and_bundle_str +
                "_dispersities.pgf", format="pgf")
    plt.gcf().clear()

    plt.figure(fig2.number)
    plt.tight_layout()
    plt.savefig(my_kr.get_run_hash_str() +
                config_and_bundle_str + "_gelfractions")
    plt.savefig(my_kr.get_run_hash_str() + config_and_bundle_str +
                "_gelfractions.pgf", format="pgf")
    plt.gcf().clear()

    # subplots of original version (ax2) and zoomed in version (in new figure)
    # save all lines of ax2 (axes can't be exchanged between figures)
    # create new subplot-figure with gridspec
    from matplotlib import gridspec
    fig = plt.figure()
    xd, yd = fig.get_size_inches()  # get default size
    fig.set_size_inches([xd * 1., yd * 2. * 0.75])
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 0.5])
    ax_new1 = plt.subplot(gs[0])
    ax_new2 = plt.subplot(gs[1])
    ax_new2.add_artist(mpl.offsetbox.AnchoredText(
        "vergrößerter Ausschnitt", loc=2, frameon=False))

    # add all the lines to the 1st and 2nd new axes
    for line in lines_gelfraction:
        ax_new1.add_line(copy.copy(line))
        ax_new2.add_line(copy.copy(line))

    # zoom in (show only certain range and keep figure size constant)
    ax_new2.set_xlim(0.4, 0.8)
    ax_new2.set_ylim(0., 0.15)
    ax_new2.grid()

    # unfortunately, the labels must be set manually
    ax_new1.set_xlabel("$X$")
    ax_new1.set_ylabel("$g$")
    ax_new2.set_xlabel("$X$")
    ax_new2.set_ylabel("$g$")

    # lines are created with Line2D, no plot() necessary
    ax_new1.set_ylim(0.0, 0.8)
    ax_new1.grid()
    ax_new2.grid()
    ax_new1.legend(loc="best")
    plt.tight_layout()
    plt.savefig(my_kr.get_run_hash_str() +
                config_and_bundle_str + "_gelfractions_zoomed")
    plt.savefig(my_kr.get_run_hash_str() + config_and_bundle_str +
                "_gelfractions_zoomed.pgf", format="pgf")
    plt.gcf().clear()


def main():
    mpl_settings()  # matplotlib settings

    my_kr = Kr()

    fig8(my_kr)


if __name__ == "__main__":
    main()
