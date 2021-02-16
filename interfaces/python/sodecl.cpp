//---------------------------------------------------------------------------//
// Copyright (c) 2015 Eleftherios Avramidis <el.avramidis@gmail.com>
//
// Distributed under The MIT License (MIT)
// See accompanying file LICENSE
//---------------------------------------------------------------------------//

#include "iostream"
#include "sodecl.hpp"
#include <cstdio>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

using namespace std;

pybind11::array_t<double> sodeclcall( std::vector<double> &a_t0,
                std::vector<double> &a_y0,
                std::vector<double> &a_params,
                int a_platform, 
                int a_device, 
                string a_system, 
                int a_solver_type, 
                int a_orbits, 
                int a_equats, 
                int a_nparams, 
                int a_nnoi, 
                double a_dt,
                double a_tspan, 
                int a_ksteps, 
                int a_local_group_size)
{
	py::array_t<double> py_array;
	sodecl::solver_Type a_solver;

    switch (a_solver_type)
    {
    case 0:
        a_solver = sodecl::solver_Type::StochasticEuler;
        break;
    case 1:
        a_solver = sodecl::solver_Type::Euler;
        break;
    case 2:
        a_solver = sodecl::solver_Type::RungeKutta;
        break;
    case 3:
        a_solver = sodecl::solver_Type::ImplicitEuler;
        break;
    case 4:
        a_solver = sodecl::solver_Type::ImplicitMidpoint;
        break;
    default:
        std::cout << "Unknown SODE solver selection." << std::endl;
        std::vector<double> myvector(0);
				py_array = py::array_t<double>({myvector.size()},{8}, &myvector[0]);
				return py_array;
    }

	sodecl::sodeclmgr *mysodeclmgr = new sodecl::sodeclmgr("kernels", 
															&a_system[0], 
															a_solver, 
															a_dt, 
															a_tspan, 
															a_ksteps, 
															a_equats, 
															a_nparams, 
															a_nnoi, 
															a_orbits, 
															sodecl::output_Type::Array);

	int success;

	// Choose device
	success = mysodeclmgr->choose_device(a_platform, sodecl::device_Type::ALL, a_device);
	if (success == 0)
	{
		cerr << "Error selecting OpenCL device!" << endl;
	}

	mysodeclmgr->set_t0(a_t0.data());
	mysodeclmgr->set_y0(a_y0.data());
	mysodeclmgr->set_params(a_params.data());

	// Set the local group size.
	mysodeclmgr->set_local_group_size(a_local_group_size);

	// Setup and run the SODE solver
	int ret = mysodeclmgr->setup_sode_solver();
		if (ret == 0)
		{
				std::vector<double> myvector(0);
				py_array = py::array_t<double>({myvector.size()},{8}, &myvector[0]);
				return py_array;
		}
    
	mysodeclmgr->run_sode_solver();

	//return mysodeclmgr->m_output;

	// Explicitly returning a numpy array here somehow prevents a lot of
	// needless copying/memory usage versus letting Python take care of it on
	// its own.
	// Not sure if we can avoid *all* copying without ditching the Py.
	py_array = py::array_t<double>({mysodeclmgr->m_output.size()},{8}, &mysodeclmgr->m_output[0]);
	return py_array;
}

PYBIND11_MODULE(sodecl_interface, m) {
    m.doc() = "sodecl plugin";

    m.def("sodeclcall", &sodeclcall, "A function that integrates a SDE/ODE system.");
}
