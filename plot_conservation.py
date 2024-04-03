#!/usr/bin/env python

import os

import numpy
import matplotlib.pyplot as plt

def parse_amr_log(path=None):

    if path is None:
        path = os.path.join(os.getcwd(), "_output", "fort.amr")

    t = []
    mass = []
    momentum = []
    mass_diff = []
    momentum_diff = []
    with open(path, 'r') as amr_file:
        for line in amr_file:
            if "time t =" in line:
                if line.split()[5] == "mass":
                    t.append(float(line.split()[3].strip(',')))
                    mass.append(float(line.split()[7]))
                    mass_diff.append(float(line.split()[10]))
                else:
                    momentum.append(float(line.split()[7]))
                    momentum_diff.append(float(line.split()[10]))

    t = numpy.array(t)
    mass = numpy.array(mass)
    momentum = numpy.array(momentum)
    mass_diff = numpy.array(mass_diff)
    momentum_diff = numpy.array(momentum_diff)

    return t, mass, momentum, mass_diff, momentum_diff

def plot_conservation():

    data = parse_amr_log()

    figs = []

    # Plot mass
    figs.append(plt.figure())
    axes = figs[0].add_subplot(1, 1, 1)
    axes.plot(data[0], data[1] / data[1][0] - 1.0, 'b-', label='Mass')
    axes.plot(data[0], data[3] / data[1][0], 'r--', label='Difference')
    axes.ticklabel_format(style="plain", useMathText=True)
    axes.set_title("Mass")
    axes.legend()

    figs.append(plt.figure())
    axes = figs[1].add_subplot(1, 1, 1)
    axes.plot(data[0], data[2] / data[2][0] - 1.0, 'b-', label='Momentum')
    axes.plot(data[0], data[4] / data[2][0], 'r--', label='Difference')
    axes.ticklabel_format(style="plain", useMathText=True)
    axes.set_title("Momentum")
    axes.legend()

    print("Mass Max, Diff = (%s, %s)" % (numpy.max(data[1]), numpy.max(data[3])))
    print("Momentum Max, Diff = (%s, %s)" % (numpy.max(data[2]), numpy.max(data[4])))

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