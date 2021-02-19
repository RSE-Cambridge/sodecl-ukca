# ---------------------------------------------------------------------------//
#  Copyright (c) 2018 Eleftherios Avramidis <el.avramidis@gmail.com>
#
#  Distributed under The MIT License (MIT)
#  See accompanying file LICENSE
# ---------------------------------------------------------------------------//

import random
import numpy
import sodecl
import matplotlib.pyplot as plt

if __name__ == '__main__':

    openclplatform = 0
    opencldevice = 0
    openclkernel = 'ukca.cl'
    solver = 3
    dt = 600
    tspan = 5*86400
    ksteps = 1
    localgroupsize = 0
    J1a = 1e-30
    J4a = 1e-30
    M = 2.5e19
    O2 = 0.2*M
    H20 = 1e17

    orbits = [1000, 10000, 100000]
    nequat = 8
    nparams = 11
    nnoi = 0

    for orbit in orbits:
        # Initial conditions
        initx = numpy.ndarray((orbit, nequat))
        for o in range(orbit):
            initx[o][0] = 1e11 #NO2
            initx[o][1] = 1e9 #NO
            initx[o][2] = 1e12 #O3
            initx[o][3] = 1e3 #O
            initx[o][4] = 1e3 #OH
            initx[o][5] = 1e5 #HO2
            initx[o][6] = 1e13 #CO
            initx[o][7] = 100 #O1D


        # Parameters values
        params = numpy.ndarray((orbit, nparams))
        for o in range(orbit):
            params[o][0] = J1a # should be J1a + 1e-2*photons
            params[o][1] = 6e-34*M*O2 # k2
            params[o][2] = 1e-11 # k3
            params[o][3] = J4a # should be J4a + 1e-5*photons 
            params[o][4] = 1e-11*M # k5
            params[o][5] = 1e-10*H20 # k6
            params[o][6] = 1e-13 # k7
            params[o][7] = 8.5e-12 # k8
            params[o][8] = 1.4e-11 # k9
            params[o][9] = 1e6 # should be 1e6*photons
            params[o][10] = H20 # H20

        t, results = sodecl.sodecl(openclplatform, opencldevice, openclkernel,
                                initx, params, solver,
                                orbit, nequat, nnoi,
                                dt, tspan, ksteps, localgroupsize)

    #print(results)
    #print(t)
    fig,axs = plt.subplots(3,3)
    axs[0,0].plot(t,results[0*orbit,:])
    axs[0,1].plot(t,results[1*orbit,:])
    axs[0,2].plot(t,results[2*orbit,:])
    axs[1,0].plot(t,results[3*orbit,:])
    axs[1,1].plot(t,results[4*orbit,:])
    axs[1,2].plot(t,results[5*orbit,:])
    axs[2,0].plot(t,results[6*orbit,:])
    axs[2,1].plot(t,results[7*orbit,:])


    sunlight = [0,0,0,0,0,0.5,1,1.5,2,2.5,3,4,5,4,3,2.5,2,1.5,1,0.5,0,0,0,0]
    photons_plot = []
    for ti in t:
        sun_index = int((ti/3600) % 24)
        photons_now = sunlight[sun_index]
        photons_plot.append(photons_now)

    axs[2,2].plot(t,photons_plot)
    plt.show()
