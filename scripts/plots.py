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

hists = {}
hists2 = {}
particles = ["K", "Pi", "Mu", "P"]
pidcuts = ["UBDT>0.25", "UBDT>0.65"]
cuts2 = ["UBDT>0.00", "UBDT>0.05", "UBDT>0.10", "UBDT>0.15", "UBDT>0.20", "UBDT>0.25", "UBDT>0.30", "UBDT>0.35", "UBDT>0.40", "UBDT>0.45", "UBDT>0.50", "UBDT>0.55", "UBDT>0.60", "UBDT>0.65", "UBDT>0.70", "UBDT>0.75", "UBDT>0.80", "UBDT>0.85", "UBDT>0.90", "UBDT>0.95", "UBDT>1.00"]
mags = ["up", "down"]
vars = ["Brunel_P", "Brunel_PT"]
dirs = ["pidcalib_output", "pidcalib_output_precut", "pidcalib_output_sy"]

for mag in mags:
    for particle in particles:
        for pidcut in pidcuts:
            for var in vars:
                for dir in dirs:
                    with open(
                            f"../efficiency_gen/{dir}/effhists-Turbo16-{mag}-{particle}-{pidcut}-{var}.pkl", "rb"
                    ) as f:
                        hists[f"{particle}_{mag}_{pidcut}_{var}_{dir}"] = pickle.load(f)

for mag in mags:
    for particle in particles:
        for cut in cuts2:
            for var in vars:
                with open(
                        f"../efficiency_gen/pidcalib_output_many/effhists-Turbo16-{mag}-{particle}-{cut}-{var}.pkl", "rb"
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
i=0
for dir in dirs:
    for mag in mags:
        for var in vars:        
            plt.figure(i)
            i+=1
            for particle in particles:
                for pidcut in pidcuts:
                    name=particle+"_"+mag+"_"+pidcut+"_"+var+"_"+dir
                    plt.hist(
                        hists[name].axes[0].edges[:-1],
                        bins=hists[name].axes[0].edges,
                        weights=hists[name].values(),
                        histtype="stepfilled",
                        label=name.replace("_", " ").replace("Pi", r"$\pi$").replace("Mu", r"$\mu$"),
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
                plt.savefig("eff_"+var+"_"+mag+"_"+dir+plots_format)
            
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
        Mu_eff_down = [
            hists2[f"passing_Mu_down_{cut}_{var}"].sum().value / hists2[f"total_Mu_down_{cut}_{var}"].sum().value
            for cut in cuts2
        ]
        bg_eff_up = [
            1-hists2[f"passing_{particle}_up_{cut}_{var}"].sum().value / hists2[f"total_{particle}_up_{cut}_{var}"].sum().value
            for cut in cuts2
        ]
        bg_eff_down = [
            1-hists2[f"passing_{particle}_down_{cut}_{var}"].sum().value / hists2[f"total_{particle}_down_{cut}_{var}"].sum().value
            for cut in cuts2
        ]
        plt.plot(Mu_eff_up, bg_eff_up,
                 "s-", markersize=8,
                 color=colors[particle+"_UBDT>0.25"],
                 label=particle+" Up")
        plt.plot(Mu_eff_down, bg_eff_down,
                 ".-",
                 color=colors[particle+"_UBDT>0.25"],
                 label=particle+" Down")
        plt.xlim(0, 1.05)
        plt.ylim(0, 1.05)
        plt.xlabel("Signal Efficiency")
        plt.ylabel("Background Rejection Efficiency")
        plt.figtext(0.2, 0.8, "LHCb\n $\sqrt{s}$=13 TeV 2016 Validation")
        plt.legend(bbox_to_anchor=(0.02, 0.8), loc="upper left")
    if plots_save:
        plt.savefig("rej_v_eff_"+var+plots_format)
                
