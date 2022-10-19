#!/usr/bin/env python3

# Standard includes
import pickle
import numpy as np
import boost_histogram as bh
import matplotlib as mpl
import matplotlib.pyplot as plt

# Style setup
# import seaborn as sns

# sns.set_palette("muted")
# sns.set_color_codes()
# sns.set_style("ticks")
# sns.set_style({"xtick.direction": "in", "ytick.direction": "in"})
# sns.set_style({"axes.grid": "True", "grid.color": "0.95"})

# plt.rcParams["figure.figsize"] = [6, 6]
# plt.rcParams["figure.dpi"] = 100
# plt.rcParams["axes.formatter.min_exponent"] = 0

import mplhep as hep

hep.style.use("LHCb2")

plt.rcParams["font.size"] = 16
plt.rcParams["figure.dpi"] = 50  # Comment out/set to 300 for production plots
plt.rcParams["axes.formatter.min_exponent"] = 0

# Config
particles = ["K", "Pi", "P", "Mu_nopt"]
particles_label = {"K": "$K$", "Pi": "$\\pi$", "Mu_nopt": "$\\mu$", "P": "$p$"}
pidcuts = ["UBDT>0.25"]
mags = ["up"]
vars = ["Brunel_P", "Brunel_PT"]
vars_label = {"Brunel_P": "$p$", "Brunel_PT": "$p_T$"}
dirs = ["pidcalib_ubdt_eff"]
output = "plots"

# To be filled
hists = {}
hists2 = {}
cuts2 = []
cuts3 = []
#  for cut in np.linspace(0, 1, 21):
for cut in [0.25]:
    cut = format(cut, '.2f')
    cuts2.append(f"UBDT>{cut}")
    cuts3.append(f"Brunel_MC15TuneV1_ProbNNmu>{cut}")

# Open files and obtain all the histograms
for mag in mags:
    for particle in particles:
        for pidcut in pidcuts:
            for var in vars:
                for dir in dirs:
                    with open(
                            f"../efficiency_gen/{dir}/effhists-Turbo16-{mag}-{particle}-{pidcut}-{var}.pkl", "rb"
                    ) as f:
                        hists[f"{particle}_{mag}_{pidcut}_{var}_{dir}"] = pickle.load(f)

                    for cut in cuts3:
                        with open(
                                f"../efficiency_gen/{dir}/effhists-Turbo16-{mag}-{particle}-{cut}-{var}.pkl", "rb"
                        ) as f:
                            hists2[f"eff_{particle}_{mag}_{cut}_{var}"] = pickle.load(f)
                            hists2[f"passing_{particle}_{mag}_{cut}_{var}"] = pickle.load(f)
                            hists2[f"total_{particle}_{mag}_{cut}_{var}"] = pickle.load(f)
                    for cut in cuts2:
                        with open(
                                f"../efficiency_gen/{dir}/effhists-Turbo16-{mag}-{particle}-{cut}-{var}.pkl", "rb"
                        ) as f:
                            hists2[f"eff_{particle}_{mag}_{cut}_{var}"] = pickle.load(f)
                            hists2[f"passing_{particle}_{mag}_{cut}_{var}"] = pickle.load(f)
                            hists2[f"total_{particle}_{mag}_{cut}_{var}"] = pickle.load(f)
plots_save = True
plots_format = ".pdf"

colors = {
    "Mu_nopt_UBDT>0.25": "xkcd:blue",
    "K_UBDT>0.25": "xkcd:red",
    "Pi_UBDT>0.25": "xkcd:green",
    "P_UBDT>0.25": "xkcd:purple",

}

particle_label = {
    'Mu': r'$\mu$',
    'Mu_nopt': r'$\mu$',
    'K': '$K$',
    'Pi': r'$\pi$',
    'P': '$p$',
}

# Efficiency as a function of a variable, with binning
i=0
for dir in dirs:
    for mag in mags:
        for var in vars:
            plt.figure(i)
            i+=1
            for particle in particles:
                for pidcut in pidcuts:
                    name=particle+"_"+mag+"_"+pidcut+"_"+var+"_"+dir
                    pidcut_lbl = pidcut.replace("UBDT", r"\mu$BDT$")
                    name_label = particle_label[particle]+fr", ${pidcut_lbl}$"
                    plt.hist(
                        hists[name].axes[0].edges[:-1],
                        bins=hists[name].axes[0].edges,
                        weights=hists[name].values(),
                        histtype="stepfilled",
                        label=name_label.replace("_", " ").replace("Pi", r"$\pi$").replace("Mu", r"$\mu$"),
                        color=colors[particle+"_"+pidcut],
                        edgecolor=colors[particle+"_"+pidcut],
                        linewidth=1.5,
                        fc=(*mpl.colors.to_rgb(colors[particle+"_"+pidcut]), 0.03),
                    )
                    plt.ylim(0, 1.6)
                    plt.margins(x=-0.01)
                    plt.legend()
                    plt.legend(fontsize=20)
                    plt.xlabel(vars_label[var]+" [MeV]")
                    plt.ylabel("Efficiency")
                    plt.figtext(0.15, 0.85, "LHCb $\\sqrt{s} = 13$ TeV\n2016 MagUp")
                    plt.figtext(0.15, 0.8, "IsMuon & MuonUnbiased & DLL$\\mu$ > 2")
            if plots_save:
                plt.savefig(output+"/eff_"+var+"_"+mag+"_"+dir+plots_format)

# Efficiency curves
for var in vars:
    plt.figure(i)
    i+=1
    for particle in particles:
        #  if particle == "Mu":
        #      continue
        Mu_eff_up = [
            hists2[f"passing_Mu_nopt_up_{cut}_{var}"].sum().value / hists2[f"total_Mu_nopt_up_{cut}_{var}"].sum().value
            for cut in cuts2
        ]
        bg_eff_up = [
            1-hists2[f"passing_{particle}_up_{cut}_{var}"].sum().value / hists2[f"total_{particle}_up_{cut}_{var}"].sum().value
            for cut in cuts2
        ]

        Mu_eff_up_pnnmu = [
            hists2[f"passing_Mu_nopt_up_{cut}_{var}"].sum().value / hists2[f"total_Mu_nopt_up_{cut}_{var}"].sum().value
            for cut in cuts3
        ]
        bg_eff_up_pnnmu = [
            1-hists2[f"passing_{particle}_up_{cut}_{var}"].sum().value / hists2[f"total_{particle}_up_{cut}_{var}"].sum().value
            for cut in cuts3
        ]
        plt.plot(Mu_eff_up, bg_eff_up,
                 "s-", markersize=8,
                 color=colors[particle+"_UBDT>0.25"],
                 label=particle_label[particle])

        plt.xlim(0, 1.05)
        plt.ylim(0, 1.05)
        plt.xlabel("Signal efficiency")
        plt.ylabel("Background rejection efficiency")
        plt.figtext(0.15, 0.25, "LHCb $\\sqrt{s} = 13$ TeV\n2016 MagUp")
        plt.figtext(0.15, 0.2, "IsMuon & MuonUnbiased & DLL$\\mu$ > 2")
        plt.legend(bbox_to_anchor=(0.02, 0.96), loc="upper left", fontsize = 18)
    if plots_save:
        plt.savefig(output+"/rej_v_eff_unbiased_"+var+plots_format)
