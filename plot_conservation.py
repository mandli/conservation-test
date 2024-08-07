#!/usr/bin/env python

import os

import numpy
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

    t = numpy.array(t)
    mass = numpy.array(mass)
    momentum = numpy.array(momentum)
    mass_diff = numpy.array(mass_diff)
    momentum_diff = numpy.array(momentum_diff)

    return [t, mass, momentum, mass_diff, momentum_diff]


# def plot_conservation(experiments, base_path=None):

#     if base_path is None:
#         base_path = os.getcwd()

#     t2boundary = 100e3 / numpy.sqrt(9.8 * 3000)

#     fig, axes = plt.subplots(3, len(experiments))
#     if len(experiments) == 1:
#         print(axes.shape)
#     titles = ['Mass', 'Momentum', r"$|\vec{u}|^2$"]
#     for (j, experiment) in enumerate(experiments):

#         # experiment_string = f"{experiment[0]}_{experiment[1]}_{experiment[2]}"
#         # path = os.path.join(base_path, experiment_string + "_output")
#         path = os.path.join(base_path, "_output")
#         data  = parse_amr_log(path)
#         t = data[0]

#         # Modify mass plot to do difference
#         data[1] = data[1] / data[1][0] - 1.0

#         for (i, field_title) in enumerate(titles):
#             axes[i][j].plot(t, data[i + 1], 'b-', label=field_title)
#             axes[i][j].ticklabel_format(style="plain", useMathText=True)
#             axes[i][j].set_title(field_title + " - " + experiment_string)
#             axes[i][j].grid(True)

#     return fig


def plot_conservation(base_path=None):

    if base_path is None:
        base_path = os.getcwd()

    t2boundary = 100e3 / numpy.sqrt(9.8 * 3000)

    fig, axs = plt.subplots(2, sharex=True)
    titles = ['Mass', 'Momentum']
    path = os.path.join(base_path, "_output")
    data  = parse_amr_log(path)
    t = data[0]

    # Modify mass plot to do difference
    data[1] = data[1] / data[1][0] - 1.0

    for (i, field_title) in enumerate(titles):
        axs[i].plot(t, data[i + 1], 'b-', label=field_title)
        axs[i].ticklabel_format(style="plain", useMathText=True)
        axs[i].set_title(field_title)
        axs[i].grid(True)
        axs[i].set_xlim([0, 20])

    return fig


if __name__ == '__main__':
    # base_path = os.path.join(os.environ.get('DATA_PATH', os.getcwd()),
    #                          "conservation_tests")
    experiments = []
    # Experiment 1
    # for init_condition in ['hump']:
    #     for order in [1, 2]:
    #         for transverse_waves in [0, 1, 2]:
    #             experiments.append([order, transverse_waves, init_condition])

    # Experiment 2
    # for init_condition in ['hump', 'pressure']:
    #     for order in [1, 2]:
    #         for transverse_waves in [0]:
    #             experiments.append([order, transverse_waves, init_condition])

    # Experiment 3
    # for init_condition in ['step',]:
    #     for order in [1]:
    #         for transverse_waves in [0]:
    #             experiments.append([order, transverse_waves, init_condition])

    # plot_conservation(experiments, base_path=os.getcwd())

    # Direct Local
    plot_conservation()

    plt.show()