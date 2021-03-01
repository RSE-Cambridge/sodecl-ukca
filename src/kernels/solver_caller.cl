//---------------------------------------------------------------------------//
// Copyright (c) 2015 Eleftherios Avramidis <el.avramidis@gmail.com>
//
// Distributed under The MIT License (MIT)
// See accompanying file LICENSE
//---------------------------------------------------------------------------//


// First dimension -> simulated time
// Second dimension -> equations
//
// Example: time is 720 10-min steps, 8 equations
//0    7
//8    15
//8*719 8*719+7

// i = simulation step index
// j = equation index
// k = thread id
int get_y0_index(int i, int j, int k) {

// within a timestep, this is your position - offset thread_id * numeq, plus equation index
// to compute the timestep, you need to add step * (numeq * numeq)
    return  (i * _listsize_  +  k) * _numeq_  + j;

}
__kernel void solver_caller(__global double *t0,
	__global double *y0,
	__global double *params_g)

{
	int i = get_global_id(0);

	double t;
	double y[_numeq_];
	double detterm[_numeq_];
	double params[_numpar_];


	t = t0[i]; 

	int k = i * _numeq_;

	k = i * _numpar_;
	for (int m = 0; m < _numpar_; ++m)
	{
		params[m] = params_g[k + m];
	}

    for (int ieq = 0; ieq < _numeq_; ++ieq)
     {
         y[ieq] = y0[get_y0_index(0, ieq, i)];
     }
	for (int it = 0; it < _numsteps_; ++it)
	{
		sode_solver(_m_dt_, t, y, detterm, params);

		t = t + _m_dt_;
		for (int ieq = 0; ieq < _numeq_; ++ieq)
		{
			y[ieq] = detterm[ieq];
            y0[get_y0_index(it, ieq, i)] = y[ieq];
		}

        t0[i] = t;

    }
}
