#!/usr/bin/env python

import os
import numpy
import datetime

import batch.batch

import clawpack.geoclaw.topotools as topotools

days2seconds = lambda days: days * 60.0**2 * 24.0

# def eta(x, y, A=1.0, sigma=50e3, x0=0.0):
#     return A * numpy.exp(-(x - x0)**2 / sigma**2)


class ConservationJob(batch.batch.Job):
    r""""""

    def __init__(self, order, transverse_waves, init_condition="pressure", 
                       base_path='./'):

        super(ConservationJob, self).__init__()

        self.type = "conservation_tests"
        self.name = ""
        self.prefix = "%s_%s_%s" % (order, transverse_waves, init_condition)
        self.executable = "xgeoclaw"

        self.init_condition = init_condition
        self.qinit_path = None

        # Create base data object
        import setrun
        self.rundata = setrun.setrun()

        self.rundata.clawdata.order = order
        self.rundata.clawdata.transverse_waves = transverse_waves

        if self.init_condition.lower() == "pressure":
            # This is the default implementation so nothing right now needs
            # to be done
            self.rundata.surge_data.pressure_forcing = True
        elif self.init_condition.lower() == "hump":
            self.rundata.qinit_data.qinit_type = 4
            self.qinit_path = os.path.abspath(os.path.join(base_path, 'hump.xyz'))
            self.rundata.qinit_data.qinitfiles = [[self.qinit_path]]
        else:
            pass # No perturbation


    def __str__(self):
        output = super(ConservationJob, self).__str__()
        output += "\n  Order: %s" % self.rundata.clawdata.order
        output += "\n  Trans: %s" % self.rundata.clawdata.transverse_waves
        output += "\n  init: %s" % self.init_condition
        return output


    def write_data_objects(self):
        r""""""

        # Write out all data files
        super(ConservationJob, self).write_data_objects()

        # Write out qinit data if need be
        if self.init_condition == "hump":
            N = 400
            qinit = topotools.Topography(topo_func=self.rundata.eta)
            qinit.x = numpy.linspace(self.rundata.clawdata.lower[0], 
                                     self.rundata.clawdata.upper[0], N)
            qinit.y = numpy.linspace(self.rundata.clawdata.lower[1], 
                                     self.rundata.clawdata.upper[1], N)
            qinit.write(self.qinit_path, topo_type=1)
        elif self.init_condition == "step":
            pass


if __name__ == '__main__':

    jobs = []
    for init_condition in ['step']:
        for order in [1, 2]:
            for transverse_waves in [0]:
                jobs.append(ConservationJob(order, transverse_waves, 
                                            init_condition=init_condition))

    controller = batch.batch.BatchController(jobs)
    print(controller)
    controller.run()
