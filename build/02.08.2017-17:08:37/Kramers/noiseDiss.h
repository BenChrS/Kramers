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

#include "simulationOptions.h"

using namespace std;
using namespace std::tr1;

// no noise
void generateNoNoise(vector<double>& output);

// generate a sequence of white noise
void generateWhiteNoise (const double& dt, gsl_rng* r,  vector<double>& output);

// generate a sequence of colored noise
void generateColoredNoise(const double& dt, gsl_rng* r, const vector<double>& G, vector<double>& colNoise);

// generate exponential correlation function
//I(t)=0.5*D/tau*exp(-t/tau)
vector<double> generateCorrFunc0(const double& t0, const double& dt, const int& nSteps, const double& D, const double& tau);

// generate Gaussian correlation function
//I(t)=1/(a*sqrt(Pi))*e^(-(t/a)^2)
vector<double> generateCorrFunc1(const double& t0, const double& dt, const int& nSteps, const double& D, const double& a);

// generate fancy correlation function
vector<double> generateCorrFunc2(const double& t0, const double& dt, const int& nSteps, const double& D, const double& chi);

// generate correlation function for massless theory
vector<double> generateCorrFunc3(const double& t0, const double& dt, const int& nSteps, const double& D, const double& alpha);

// generate dissipation kernel
//Gamma(t)=0.5/(k*T)*I(t) -- dissipation-fluctuation theorem
vector<double> generateDissipationKernel(const int& nMax, const vector<double>& corrFunction,
					 const double& temperature, const double& k, const bool& optimisationBool, const double& dPrecision);

// calculate the function G(t)
void calcG(const double& dt, const vector<double>& corrFunc, const int& corrFuncNr,
		     const bool& timeSaving = false, const double& gPrecision = pow(10,-3));

// initializes random generator
void initRandGen(const long& seed, vector<gsl_rng *>& randGen);

class NoiseDiss{
public:
	function< vector<double>(const double&, const double&, const int&) > corrFunc;
	function<void(const double&, gsl_rng*, vector<double>&)> noiseFunc;
	  // ----prepare simulation functions and dependent variables----
	  vector<double> tFourierVec;  // all times for nFourier steps

	  vector<double> corrFVec; // I(t) for t=[0, nFourier*dt]
	  vector<double> dispKernel; //Gamma(t) for t=[0, nSteps*dt]
	  vector<double> noise;
	  vector< vector<double> > allNoise;
	  vector<double>G;
	  vector<double> tGVec;
	  long seed; // seed for random generator
	  vector<gsl_rng *> randGen;


	NoiseDiss(SimulationOptions& so);
	//int storeAsXML(const string& foldername);
	//int loadFromXML(const string& filename);
};
