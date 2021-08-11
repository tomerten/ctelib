from glob import glob

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt


def main():
    file = "./cpp/tests/bin/Distribution_bucket_0_11_08_2021.dat"
    df = pd.read_csv(
        file, delim_whitespace=True, header=None, names=["id", "x", "px", "y", "py", "t", "delta"]
    )
    df["turn"] = df.groupby("id").cumcount() + 1
    dfg = df.groupby("turn")

    emitfile = "./cpp/tests/bin/CTE_Emittances_0_11_08_2021.dat"
    dfe = pd.read_csv(
        emitfile, delim_whitespace=True, header=None, names=["id", "ex", "ey", "sigs", "sige"]
    )
    plt.ion()
    fig = plt.figure(constrained_layout=True, figsize=(14, 10))

    # create plot grid
    gs = fig.add_gridspec(2, 4)

    # create subplots and set titles
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[0, 2:])
    ax4 = fig.add_subplot(gs[1, 0])
    ax5 = fig.add_subplot(gs[1, 1])
    ax6 = fig.add_subplot(gs[1, 2])
    ax7 = fig.add_subplot(gs[1, 3])

    ax1.set_title(r"$x-p_x$ init")
    ax2.set_title(r"$y-p_y$ init")
    ax3.set_title(r"$t-d\gamma$ init")
    ax4.set_title(r"L ip1")
    ax5.set_title(r"L ip2")
    ax6.set_title(r"L ip5")
    ax7.set_title(r"L ip8")

    ex_plot = ax4.plot([], [], marker=".", linestyle="-")[0]
    ey_plot = ax5.plot([], [], marker=".", linestyle="-")[0]
    sigs_plot = ax6.plot([], [], marker=".", linestyle="-")[0]

    for turn, g in dfg:
        ax1.scatter(g["x"], g["px"], s=3)
        ax2.scatter(g["y"], g["py"], s=3)
        ax3.scatter(g["t"].values, g["delta"].values, s=3)
        ex_plot.set_xdata(np.append(ex_plot.get_xdata(), [turn]))
        ex_plot.set_ydata(np.append(ex_plot.get_ydata(), [dfe.loc[dfe["id"] == turn - 1, "ex"]]))

        ey_plot.set_xdata(np.append(ey_plot.get_xdata(), [turn]))
        ey_plot.set_ydata(np.append(ey_plot.get_ydata(), [dfe.loc[dfe["id"] == turn - 1, "ey"]]))

        sigs_plot.set_xdata(np.append(sigs_plot.get_xdata(), [turn]))
        sigs_plot.set_ydata(
            np.append(sigs_plot.get_ydata(), [dfe.loc[dfe["id"] == turn - 1, "sigs"]])
        )
        plt.draw()
        ax1.set_xlim(-0.001, 0.001)
        ax2.set_xlim(-0.001, 0.001)
        ax3.set_xlim(0, 1.5e-9)
        ax1.set_ylim(-0.001, 0.001)
        ax2.set_ylim(-0.001, 0.001)
        ax3.set_ylim(-1e-2, 1e-2)
        ax4.set_xlim(0, 100)
        ax4.set_ylim(0, 1e-8)
        ax5.set_xlim(0, 100)
        ax5.set_ylim(0, 1e-9)
        ax6.set_xlim(0, 100)
        ax6.set_ylim(0, 1e-2)
        ax7.set_xlim(0, 100)
        plt.pause(0.01)
        ax1.cla()
        ax2.cla()
        ax3.cla()

        # plt.show()
    plt.ioff()


if __name__ == "__main__":
    main()
