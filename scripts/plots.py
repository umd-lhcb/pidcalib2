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

plt.rcParams["font.size"] = 20
plt.rcParams["figure.dpi"] = 50  # Comment out/set to 300 for production plots
plt.rcParams["axes.formatter.min_exponent"] = 0

# Config
particles = ["K", "Pi", "Mu"]
pidcuts = ["UBDT>0.25", "UBDT>0.65"]
mags = ["up"]
vars = ["Brunel_P", "Brunel_PT"]
dirs = ["pidcalib_output_many_8_8"]
output = "plots_9_20_old"

# To be filled
hists = {}
hists2 = {}
cuts2 = []
cuts3 = []
for cut in np.linspace(0, 1, 21): 
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
    "Mu_UBDT>0.25": "xkcd:blue",
    "Mu_UBDT>0.65": "xkcd:pastel blue",
    "K_UBDT>0.25": "xkcd:red",
    "K_UBDT>0.65": "xkcd:pink",
    "Pi_UBDT>0.25": "xkcd:green",
    "Pi_UBDT>0.65": "xkcd:pastel green",
    "P_UBDT>0.25": "xkcd:purple",
    "P_UBDT>0.65": "xkcd:pastel purple",
    
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
                    name_label = particle+"_"+mag+"_"+pidcut+"_"+var
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
                    plt.xlabel(var+" [MeV/c]")
                    plt.ylabel("Efficiency")
                    plt.figtext(0.2, 0.8, "LHCb\n $\\sqrt{s}$=13 TeV 2016 Validation")
            if plots_save:
                plt.savefig(output+"/eff_"+var+"_"+mag+"_"+dir+plots_format)

# Efficiency curves                
for var in vars:        
    plt.figure(i)
    i+=1
    for particle in particles:
        if particle == "Mu":
            continue
        Mu_eff_up = [
            hists2[f"passing_Mu_up_{cut}_{var}"].sum().value / hists2[f"total_Mu_up_{cut}_{var}"].sum().value
            for cut in cuts2
        ]
        bg_eff_up = [
            1-hists2[f"passing_{particle}_up_{cut}_{var}"].sum().value / hists2[f"total_{particle}_up_{cut}_{var}"].sum().value
            for cut in cuts2
        ]

        Mu_eff_up_pnnmu = [
            hists2[f"passing_Mu_up_{cut}_{var}"].sum().value / hists2[f"total_Mu_up_{cut}_{var}"].sum().value
            for cut in cuts3
        ]
        bg_eff_up_pnnmu = [
            1-hists2[f"passing_{particle}_up_{cut}_{var}"].sum().value / hists2[f"total_{particle}_up_{cut}_{var}"].sum().value
            for cut in cuts3
        ]
        plt.plot(Mu_eff_up, bg_eff_up,
                 "s-", markersize=8,
                 color=colors[particle+"_UBDT>0.25"],
                 label=particle+" Up UBDT")
        plt.plot(Mu_eff_up_pnnmu, bg_eff_up_pnnmu,
                 "s-", markersize=8,
                 color=colors[particle+"_UBDT>0.65"],
                 label=particle+" Up ProbNNmu")
        # Mu_eff_down = [
        #     hists2[f"passing_Mu_down_{cut}_{var}"].sum().value / hists2[f"total_Mu_down_{cut}_{var}"].sum().value
        #     for cut in cuts2
        # ]
        # bg_eff_down = [
        #     1-hists2[f"passing_{particle}_down_{cut}_{var}"].sum().value / hists2[f"total_{particle}_down_{cut}_{var}"].sum().value
        #     for cut in cuts2
        # ]
        # Mu_eff_down_pnnmu = [
        #     hists2[f"passing_Mu_down_{cut}_{var}"].sum().value / hists2[f"total_Mu_down_{cut}_{var}"].sum().value
        #     for cut in cuts3
        # ]
        # bg_eff_down_pnnmu = [
        #     1-hists2[f"passing_{particle}_down_{cut}_{var}"].sum().value / hists2[f"total_{particle}_down_{cut}_{var}"].sum().value
        #     for cut in cuts3
        # ]
        # plt.plot(Mu_eff_down, bg_eff_down,
        #          ".-",
        #          color=colors[particle+"_UBDT>0.25"],
        #          label=particle+" Down UBDT")
        # plt.plot(Mu_eff_down_pnnmu, bg_eff_down_pnnmu,
        #          ".-",
        #          color=colors[particle+"_UBDT>0.65"],
        #          label=particle+" Down ProbNNmu")
        plt.xlim(0, 1.05)
        plt.ylim(0, 1.05)
        plt.xlabel("Signal Efficiency")
        plt.ylabel("Background Rejection Efficiency")
        plt.figtext(0.2, 0.2, "IsMuon==1 & MuonUnbiased==1 & DLLmu>2 Online")
        plt.legend(bbox_to_anchor=(0.02, 0.8), loc="upper left", fontsize = 20)
    if plots_save:
        plt.savefig(output+"/rej_v_eff_unbiased_"+var+plots_format)
                
# for var in vars:        
#     plt.figure(i)
#     i+=1
#     Mu_eff_up = [
#         hists2[f"passing_Mu_up_{cut}_{var}"].sum().value / hists2[f"total_Mu_up_{cut}_{var}"].sum().value
#         for cut in cuts2
#     ]
#     bg_eff_up = [
#         1-0.8*(hists2[f"passing_Pi_up_{cut}_{var}"].sum().value / hists2[f"total_Pi_up_{cut}_{var}"].sum().value)-0.15*(hists2[f"passing_K_up_{cut}_{var}"].sum().value / hists2[f"total_K_up_{cut}_{var}"].sum().value)-0.05*(hists2[f"passing_P_up_{cut}_{var}"].sum().value / hists2[f"total_P_up_{cut}_{var}"].sum().value)
#         for cut in cuts2
#     ]
    
#     Mu_eff_up_pnnmu = [
#         hists2[f"passing_Mu_up_{cut}_{var}"].sum().value / hists2[f"total_Mu_up_{cut}_{var}"].sum().value
#         for cut in cuts3
#     ]
#     bg_eff_up_pnnmu = [
#         1-0.8*(hists2[f"passing_Pi_up_{cut}_{var}"].sum().value / hists2[f"total_Pi_up_{cut}_{var}"].sum().value)-0.15*(hists2[f"passing_K_up_{cut}_{var}"].sum().value / hists2[f"total_K_up_{cut}_{var}"].sum().value)-0.05*(hists2[f"passing_P_up_{cut}_{var}"].sum().value / hists2[f"total_P_up_{cut}_{var}"].sum().value)
#         for cut in cuts3
#     ]
#     plt.plot(Mu_eff_up, bg_eff_up,
#              "s-", markersize=8,
#              color=colors["Mu_UBDT>0.25"],
#              label="Bkg Up UBDT")
#     plt.plot(Mu_eff_up_pnnmu, bg_eff_up_pnnmu,
#              "s-", markersize=8,
#              color=colors["Mu_UBDT>0.65"],
#              label="Bkg Up ProbNNmu")
#     plt.xlim(0.695, 1.05)
#     plt.ylim(0.295, 1.05)
#     plt.xlabel("Signal Efficiency")
#     plt.ylabel("Background Rejection Efficiency")
#     plt.figtext(0.2, 0.2, "IsMuon==1 & MuonUnbiased==1 & DLLmu>2 Offline")
#     plt.legend(bbox_to_anchor=(0.02, 0.8), loc="upper left", fontsize = 20)
#     if plots_save:
#         plt.savefig("plots/rej_v_eff_weighted_"+var+plots_format)
                
