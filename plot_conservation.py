#!/usr/bin/env python

import os
import sys

import numpy as np
import matplotlib.pyplot as plt

import clawpack.geoclaw.surge.plot as surgeplot

def parse_amr_log(path=None):

    if path is None:
        path = os.path.join(os.getcwd(), "_output", "fort.amr")
    else:
        path = os.path.join(path, "fort.amr")

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
                else:
                    raise ValueError("Invalid type of conservation found.")

    t = np.array(t)
    mass = np.array(mass)
    momentum = np.array(momentum)
    mass_diff = np.array(mass_diff)
    momentum_diff = np.array(momentum_diff)

    return [t, mass, momentum, mass_diff, momentum_diff]

def plot_conservation(data):
    t = data[0]
    # data[1] = data[1] / data[1][0] - 1.0

    fig, axs = plt.subplots(2, sharex=True)
    fig.suptitle("Differences")
    titles = ['Mass', 'Momentum']
    for (i, field_title) in enumerate(titles):
        axs[i].plot(t, data[i + 3], 'b-', label=field_title)
        axs[i].set_title(field_title)
        axs[i].grid(True)
        axs[i].set_xlim([t.min(), t.max()])

    return fig


if __name__ == '__main__':
    # Load data
    if len(sys.argv) <= 1:
        base_path = os.getcwd()
    else:
        base_path = sys.argv[1]
    data = parse_amr_log(os.path.join(os.getcwd(), "_output"))
    
    # Print out max diffs
    index = data[3].argmax()
    print(f"max mass difference = {data[3][index]} at t = {data[0][index]}")
    index = data[4].argmax()
    print(f"max momentum difference = {data[4][index]} at t = {data[0][index]}")

    # Plot conservation
    plot_conservation(data)

    plt.show()