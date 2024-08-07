import os

import numpy
import matplotlib.pyplot as plt
import datetime

import clawpack.visclaw.colormaps as colormap
import clawpack.visclaw.gaugetools as gaugetools
import clawpack.clawutil.data as clawutil
import clawpack.amrclaw.data as amrclaw
import clawpack.geoclaw.data as geodata


import clawpack.geoclaw.surge.plot as surgeplot

# Set indices

surgeplot.wind_field = 2
surgeplot.pressure_field = 4

try:
    from setplotfg import setplotfg
except:
    setplotfg = None


def setplot(plotdata=None):
    """"""

    if plotdata is None:
        from clawpack.visclaw.data import ClawPlotData
        plotdata = ClawPlotData()

    # clear any old figures,axes,items data
    plotdata.clearfigures()
    plotdata.format = 'ascii'

    # Load data from output
    clawdata = clawutil.ClawInputData(2)
    clawdata.read(os.path.join(plotdata.outdir, 'claw.data'))
    physics = geodata.GeoClawData()
    physics.read(os.path.join(plotdata.outdir, 'geoclaw.data'))
    surge_data = geodata.SurgeData()
    surge_data.read(os.path.join(plotdata.outdir, 'surge.data'))
    friction_data = geodata.FrictionData()
    friction_data.read(os.path.join(plotdata.outdir, 'friction.data'))

    # Load storm track
    # track = surgeplot.track_data(os.path.join(plotdata.outdir, 'fort.track'))

    # Set afteraxes function
    # def surge_afteraxes(cd):
    #     surgeplot.surge_afteraxes(cd, track, plot_direction=False,
    #                                          kwargs={"markersize": 4})

    # Color limits
    surface_limits = [-0.1, 0.1]
    speed_limits = [0.0, 0.55]

    def friction_after_axes(cd):
        plt.title(r"Manning's $n$ Coefficient")

    # ==========================================================================
    #   Plot specifications
    # ==========================================================================
    regions = {"Full Domain": {"xlimits": (clawdata.lower[0], clawdata.upper[0]),
                               "ylimits": (clawdata.lower[1], clawdata.upper[1]),
                               "figsize": (6.4, 4.8)}, 
              }

    for (name, region_dict) in regions.items():

        # Surface Figure
        plotfigure = plotdata.new_plotfigure(name="Surface - %s" % name)
        # plotfigure.kwargs = {"figsize": region_dict['figsize']}
        plotaxes = plotfigure.new_plotaxes()
        plotaxes.title = "Surface"
        plotaxes.xlimits = region_dict["xlimits"]
        plotaxes.ylimits = region_dict["ylimits"]
        # plotaxes.afteraxes = surge_afteraxes

        surgeplot.add_surface_elevation(plotaxes, bounds=surface_limits)
        surgeplot.add_land(plotaxes, bounds=[0.0, 20.0])
        plotaxes.plotitem_dict['surface'].amr_patchedges_show = [0] * 10
        plotaxes.plotitem_dict['land'].amr_patchedges_show = [0] * 10

        # Speed Figure
        plotfigure = plotdata.new_plotfigure(name="Currents - %s" % name)
        # plotfigure.kwargs = {"figsize": region_dict['figsize']}
        plotaxes = plotfigure.new_plotaxes()
        plotaxes.title = "Currents"
        plotaxes.xlimits = region_dict["xlimits"]
        plotaxes.ylimits = region_dict["ylimits"]
        # plotaxes.afteraxes = surge_afteraxes

        surgeplot.add_speed(plotaxes, bounds=speed_limits)
        surgeplot.add_land(plotaxes, bounds=[0.0, 20.0])
        plotaxes.plotitem_dict['speed'].amr_patchedges_show = [0] * 10
        plotaxes.plotitem_dict['land'].amr_patchedges_show = [0] * 10

    # Transect
    def plot_transect(current_data, y0=0.0):
        plt.gca().grid(True) 
        # y = current_data.y
        # dy = current_data.dy
        # index = numpy.where(abs(y - y0) <= dy / 2.0)[1][0]
        # x = current_data.x[:, index]
        # q = current_data.q
        # eta = q[3, :, index]
        # u = numpy.where(q[0, :] > 0, q[1, :] / q[0, :], 0.0)

        # axs = [None, None]
        # axs[0] = plt.gca()

        # eta_line = axs[0].plot(x, eta, 'b', label=r"$\eta$")
        # axs[0].set_ylabel('Surface (m)')
        # axs[0].set_ylim(surface_limits)

        # axs[1] = axs[0].twinx()
        # u_line = axs[1].plot(x, u, 'r', label=r"$u$")
        # axs[1].set_ylabel("Velocity (m/s)")
        # axs[1].set_ylim([-speed_limits[1], speed_limits[1]])
        
        # axs[0].grid(True)
        # lines = eta_line + u_line
        # axs[0].legend(lines, [line.get_label() for line in lines])


    def transect(current_data, field=3, y0=0.0):
        y = current_data.y
        dy = current_data.dy
        index = numpy.where(abs(y - y0) <= dy / 2.0)[1][0]
        return current_data.x[:, index], current_data.q[field, :, index]

    plotfigure = plotdata.new_plotfigure(name="Surface Transect")
    plotfigure.show = True
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = "Surface Transect"
    plotaxes.xlabel = "x (m)"
    plotaxes.ylabel = r"$\eta$"
    plotaxes.xlimits = [clawdata.lower[0], clawdata.upper[0]]
    plotaxes.ylimits = surface_limits
    plotaxes.afteraxes = plot_transect

    plotitem = plotaxes.new_plotitem(plot_type="1d_from_2d_data")
    plotitem.map_2d_to_1d = transect
    plotitem.plotstyle = 'ko-'
    plotitem.kwargs = {"markersize": 3}

    # ========================================================================
    #  Figures for gauges
    # ========================================================================
    plotfigure = plotdata.new_plotfigure(name='Gauge Surface and Velocity', figno=300,
                                         type='each_gauge')
    plotfigure.show = True
    plotfigure.clf_each_gauge = True

    def gauge_afteraxes(cd):

        axes = plt.gca()

        gauge = cd.gaugesoln
        t = gauge.t
        eta_line = axes.plot(t, gauge.q[3, :], 'b', label=r"$\eta$")
        axes.plot(t, numpy.zeros(t.shape), 'b-.')
        axes.set_ylabel('Surface (m)')
        axes.set_ylim([-0.2, 0.2])

        axes2 = axes.twinx()
        u = numpy.where(gauge.q[0, :] > 0, gauge.q[1, :] / gauge.q[0, :], 0.0)
        u_line = axes2.plot(t, u, 'r', label=r"$u$")
        axes2.plot(t, numpy.zeros(t.shape), 'r-.')
        axes2.set_ylabel("Velocity (m/s)")
        axes2.set_ylim([-0.5, 0.5])

        # Fix up plot - in particular fix time labels
        axes.set_title('Station %s' % cd.gaugeno)
        axes.set_xlim([clawdata.t0, clawdata.tfinal])
        # axes.set_xlabel('Days relative to landfall')
        # axes.set_xticks([-1, 0, 1, 2, 3, 4, 5])
        # axes.set_xticklabels([r"$-1$", r"$0$", r"$1$", r"$2$", r"$3$", r"$4$", r"$5$"])
        axes.grid(True)
        lines = eta_line + u_line
        axes.legend(lines, [line.get_label() for line in lines])

        plt.gcf().tight_layout()

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.afteraxes = gauge_afteraxes
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.show = False

    #
    #  Gauge Location Plot
    #
    def gauge_location_afteraxes(cd):
        plt.subplots_adjust(left=0.12, bottom=0.06, right=0.97, top=0.97)
        gaugetools.plot_gauge_locations(cd.plotdata, gaugenos='all',
                                        format_string='ko', add_labels=False)
        # gaugetools.plot_gauge_locations(cd.plotdata, gaugenos=[1, 5, 9],
        #                                 format_string='ko', add_labels=True)

    plotfigure = plotdata.new_plotfigure(name="Gauge Locations")
    plotfigure.show = True

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Gauge Locations'
    plotaxes.scaled = True
    plotaxes.xlimits = regions['Full Domain']['xlimits']
    plotaxes.ylimits = regions['Full Domain']['ylimits']
    plotaxes.afteraxes = gauge_location_afteraxes
    surgeplot.add_surface_elevation(plotaxes, bounds=surface_limits)
    surgeplot.add_land(plotaxes, bounds=[0.0, 20.0])
    plotaxes.plotitem_dict['surface'].amr_patchedges_show = [0] * 10
    plotaxes.plotitem_dict['land'].amr_patchedges_show = [0] * 10

    # -----------------------------------------
    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_gaugenos = 'all'          # list of gauges to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?
    plotdata.parallel = True                 # parallel plotting

    return plotdata
