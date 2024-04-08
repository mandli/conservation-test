#!/usr/bin/env python

import os

import numpy
import matplotlib.pyplot as plt

import clawpack.geoclaw.surge.plot as surgeplot

def parse_amr_log(path=None):

    if path is None:
        path = os.path.join(os.getcwd(), "_output", "fort.amr")

    t = []
    mass = []
    momentum = []
    u2 = []
    mass_diff = []
    momentum_diff = []
    u2_diff = []
    with open(path, 'r') as amr_file:
        for line in amr_file:
            if "time t =" in line:
                if line.split()[5] == "mass":
                    t.append(float(line.split()[3].strip(',')))
                    mass.append(float(line.split()[7]))
                    mass_diff.append(float(line.split()[10]))
                elif line.split()[5] == 'mom':
                    momentum.append(float(line.split()[7]))
                    momentum_diff.append(float(line.split()[10]))
                elif line.split()[5] == 'u**2':
                    u2.append(float(line.split()[7]))
                    u2_diff.append(float(line.split()[10]))
                else:
                    raise ValueError("Invalid type of conservation found.")

    t = numpy.array(t)
    mass = numpy.array(mass)
    momentum = numpy.array(momentum)
    u2 = numpy.array(u2)
    mass_diff = numpy.array(mass_diff)
    momentum_diff = numpy.array(momentum_diff)
    u2_diff = numpy.array(u2_diff)

    return t, mass, momentum, u2, mass_diff, momentum_diff, u2_diff

def plot_conservation():

    data = parse_amr_log()
    t = surgeplot.sec2days(data[0])

    figs = []

    # Plot mass
    figs.append(plt.figure())
    axes = figs[-1].add_subplot(1, 1, 1)
    axes.plot(t, data[1] / data[1][0] - 1.0, 'b-', label='Mass')
    axes.plot(t, data[4] / data[1][0], 'r--', label='Difference')
    axes.ticklabel_format(style="plain", useMathText=True)
    axes.set_title("Mass")
    axes.legend()

    axes.set_xlim([-1, 5])
    axes.set_xlabel('Days relative to landfall')
    axes.set_xticks([-1, 0, 1, 2, 3, 4, 5])
    axes.set_xticklabels([r"$-1$", r"$0$", r"$1$", r"$2$", r"$3$", r"$4$", r"$5$"])
    axes.grid(True)

    figs.append(plt.figure())
    axes = figs[-1].add_subplot(1, 1, 1)
    axes.plot(t, data[2], 'b-', label='Momentum')
    axes.plot(t, data[5], 'r--', label='Difference')
    axes.ticklabel_format(style="plain", useMathText=True)
    axes.set_title("Momentum")
    axes.legend()

    axes.set_xlim([-1, 5])
    axes.set_xlabel('Days relative to landfall')
    axes.set_xticks([-1, 0, 1, 2, 3, 4, 5])
    axes.set_xticklabels([r"$-1$", r"$0$", r"$1$", r"$2$", r"$3$", r"$4$", r"$5$"])
    axes.grid(True)

    figs.append(plt.figure())
    axes = figs[-1].add_subplot(1, 1, 1)
    axes.plot(t, data[3], 'b-', label=r'$|\vec{u}|^2$')
    axes.plot(t, data[6], 'r--', label='Difference')
    axes.ticklabel_format(style="plain", useMathText=True)
    axes.set_title(r"$|\vec{u}|^2$")
    axes.legend()
    
    axes.set_xlim([-1, 5])
    axes.set_xlabel('Days relative to landfall')
    axes.set_xticks([-1, 0, 1, 2, 3, 4, 5])
    axes.set_xticklabels([r"$-1$", r"$0$", r"$1$", r"$2$", r"$3$", r"$4$", r"$5$"])
    axes.grid(True)

    print("Mass Max, Diff = (%s, %s)" % (numpy.max(data[1]), numpy.max(data[4])))
    print("Momentum Max, Diff = (%s, %s)" % (numpy.max(data[2]), numpy.max(data[5])))
    print("u^2 Max, Diff = (%s, %s)" % (numpy.max(data[3]), numpy.max(data[6])))

    # # Mass difference
    # figs.append(plt.figure())
    # axes = figs[1].add_subplot(1, 1, 1)
    # axes.plot(data[0], data[3] / data[1][0])
    # axes.ticklabel_format(style="plain", useMathText=True)
    # axes.set_title("Mass Difference")

    # # Plot momentum
    # figs.append(plt.figure())
    # axes = figs[2].add_subplot(1, 1, 1)
    # axes.plot(data[0], data[2] / data[2][0])
    # axes.ticklabel_format(style="plain", useMathText=True)
    # axes.set_title("Total Momentum")

    # # Momentum difference
    # figs.append(plt.figure())
    # axes = figs[3].add_subplot(1, 1, 1)
    # axes.plot(data[0], data[4] / data[2][0])
    # axes.ticklabel_format(style="plain", useMathText=True)
    # axes.set_title("Momentum Difference")

    return figs


if __name__ == '__main__':
    plot_conservation()
    plt.show()