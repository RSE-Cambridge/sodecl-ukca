# ---------------------------------------------------------------------------//
#  Copyright (c) 2017 Eleftherios Avramidis <el.avramidis@gmail.com>
#
#  Distributed under The MIT License (MIT)
#  See accompanying file LICENSE
# ---------------------------------------------------------------------------//

import os
import subprocess
import numpy
import numpy.matlib
import platform
import sodecl_interface

# Calls the SODECL executable.
#
#  Executes the SODECL executable for specific parameters.
#
#  @param openclplatform       	is the number of the OpenCL platform.
#  @param opencldevice       	is the is the number of the OpenCL device for the selected platform.
#  @param openclkernel       	is the path to the OpenCL kernel with the definition of the SODE system.
#  @param initx       			is the array with the initial state of the SODE system.
#  @param params       			is the array with the parameters of the SODE system.
#  @param sodesolver       		is the integration solver for the SODE system.
#  @param orbits       			is the number of orbits to integrate the SODE system.
#  @param nequat       			is the number of equations of the SODE system.
#  @param nnoiseprocesses       is the number of noise processes needed by the SODE system.
#  @param dt       				is the time step for the SODE solver.
#  @param tspan       			is the time span for the integration of the SODE system.
#  @param ksteps       			is the number of SODE solver step for each call to the OpenCL device.
#  @param localgroupsize       	is the size of the OpenCL local group size.
#
#  @return              		the integration results.

def sodecl(openclplatform, opencldevice, openclkernel,
           initx, params, sodesolver,
           orbits, nequat, nnoiseprocesses,
           dt, tspan, ksteps, localgroupsize):

    t0 = numpy.array([0])
    t0 = numpy.matlib.repmat(t0, orbits, 1)
    nparams = params[0].size

    t0=t0.flatten()
    initx=initx.flatten()
    params=params.flatten()

    if nnoiseprocesses == 0:
        # solvers = {'se': 0,
        #            'e': 1,
        #            'r': 2,
        #            'ie': 3,
        #            'im': 4}

        if sodesolver>0 and sodesolver<5:
            solver2user = sodesolver
        else:
            print("Invalid SODE solver selected! Using Euler!")
            solver2user = 1
    else:
        solver2user = 0

    results = sodecl_interface.sodeclcall(  t0,
                        initx,
                        params,
                        openclplatform,
                        opencldevice,
                        openclkernel,
                        solver2user,
                        orbits,
                        nequat,
                        nparams,
                        nnoiseprocesses,
                        dt,
                        tspan,
                        ksteps,
                        localgroupsize)

    results = results.reshape(orbits*nequat, int(results.shape[0] / (orbits*nequat)), order='F')

    t=numpy.arange(0, tspan, dt)

    return t, results
