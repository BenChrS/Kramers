#pragma once

#include <stdio.h>
#include <iostream>
#include <stdio.h>

#include <tr1/functional>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>

#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>

#include "noiseDiss.h"
#include "simulationOptions.h"
#include "results.h"

using namespace std;
using namespace std::tr1;

// Adams-Bashforth scheme for white noise
//solves the Langevin equation with white noise and friction coefficient gamma
int adamsBWhiteNoise(const vector<double>& tVec, const double& dt, const double& x0, const double& v0, const double& mass,const double& potw,
		     const double& k_b, const double& temperature,const vector<double>& whiteNoise, const double gamma, 
		     const double D,vector<double>& xVec, vector<double>& velocVec,function<double(const double&, const double&)> forceExtFunc,
		     const bool& KramersRate,const double rxborder,const double lxborder,gsl_rng* res_rng,gsl_rng* noise_rng,const int& initCond,
		     const NoiseDiss&);

// Euler scheme for white noise
//solves the Langevin equation with white noise and friction coefficient gamma with the euler method
int eulerWhiteNoise(const vector<double>& tVec, const double& dt, const double& x0, const double& v0, const double& mass,
		    const double& k_b, const double& temperature,const vector<double>& whiteNoise, const double& gamma, 
		    const double& D,vector<double>& xVec, vector<double>& velocVec,function<double(const double&, const double&)> forceExtFunc,
		    const bool& KramersRate,const double rxborder,const double lxborder,gsl_rng* res_rng,gsl_rng* noise_rng,const int& initCond,
		    const NoiseDiss&);

// Euler scheme for colored noise
// solves the Langevin-differential equation with colored Noise and dissipation kernel with the euler method
// still testing
int eulerColNoise(const vector<double>& tVec, const double& dt, const double& x0, const double& v0, const double& mass,
		  const double& k_b, const double& temperature, const vector<double>& colNoise, 
		  const vector<double>& dissipationKernel, vector<double>& xVec, vector<double>& velocVec,
		  function<double(const double&, const double&)> forceExtFunc,const bool& KramersRate,const double xborder,
		  const double lxborder,gsl_rng* res_rng,gsl_rng* noise_rng,const int& intCond,const NoiseDiss&);

// Euler-Cromer scheme for colored noise, more stable for harmonic oscillations
// solves the Langevin-differential equation with colored Noise and dissipation kernel with the Euler-Cromer method
int eulerCromerColNoise(const vector<double>& tVec, const double& dt, const double& x0, const double& v0, const double& mass,
			const double& k_b, const double& temperature, const vector<double>& colNoise, const vector<double>& dissipationKernel,
		        vector<double>& xVec, vector<double>& velocVec,function<double(const double&, const double&)> forceExtFunc,
			const bool& KramersRate,const double xborder,const double lxborder,gsl_rng* res_rng,gsl_rng* noise_rng,
 		        const int& intCond,const NoiseDiss&);

// // solves the Langevin-differential equation with colored Noise and dissipation kernel with the Leap-Frog method
// // not working - difficult to apply as a=a(v,x)
// int leapFrogColNoise(const vector<double>& tVec, const double dt, const vector<double>& y0, const double mass,
// 		     const vector<double>& colNoise, const vector<double>& dissipationKernel,
// 		     vector<double>& yVec, vector<double>& velocVec,
// 		     function<double(double, double)> forceExtFunc);

// Adams-Bashforth scheme for colored noise
// solves the Langevin-differential equation with colored Noise and dissipation kernel with the Adams-Bashforth-Method
// without implicit integral solving
// energy is not perfectly conserved
// if absolute position of particle becomes larger than xCutoff, all following positions will be set to the same value
// set xCutoff to negative value for turning that feature off
int adamsBColNoise(const vector<double>& tVec, const double& dt, const double& x0, const double& v0, const double& mass,const double& potw,
		     const double& k_b, const double& temperature,const vector<double>& colNoise, const vector<double>& dissipationKernel,
		     vector<double>& xVec, vector<double>& velocVec,
		     function<double(const double&, const double&)> forceExtFunc, const double& xCutoff, 
		     const bool& KramersRate, const double& xborder,const double lxborder,gsl_rng* res_rng,gsl_rng* noise_rng,
		     const int& intCond,const NoiseDiss&);

class SdeSolver{
public:
	function<int(const vector<double>&, const double&, const double&, const double&, const vector<double>&, vector<double>&, vector<double>&, function<double(const double&, const double&)>)> sdeSolverFunc;
	SdeSolver(const SimulationOptions&,const NoiseDiss&);
};
