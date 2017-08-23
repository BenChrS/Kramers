#pragma once

#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <time.h>
#include <iterator>
#include <tr1/functional>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_histogram.h>
#include <fstream>
#include <omp.h>
#include <sstream>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <string.h>
#include <numeric>
#include <complex>

#include "simulationOptions.h"
#include "vectorCalculations.h"
#include "filesFolder.h"
#include "results.h"
#include "potential.h"
#include "noiseDiss.h"
//#include "kappasimulation.h"

using namespace std;
using namespace std::tr1;

void calcVDistribution(vector<double>& vec,vector<double>& vec1, const vector<double>& data, const double& rangeBeg, const double& rangeEnd, const int& nBins);

//averaged probability density Distribution of velocity over range from rowBeg to rowEnd
void calcvDistributionMean(vector<double>& vec,vector<double>& vec1, const vector< vector<double> >& data, const int& rowBeg, const int& rowEnd, const double& rangeBeg,const double& rangeEnd, const int& nBins);

// returns the Maxwell Boltzmann probability density distribution for one dimension v in the range [rangeBeg, rangeEnd]
void calcVDistTheo(vector<double>& vec, const double& temperature, const double& mass, const double& k, const double& rangeBeg, const double& rangeEnd, const int& nBins);

// computes the average position of the particles, note: already implemented by buildaverage in vectorCalculations
void calcxAv(vector<double>& vec,const vector< vector<double> >& data);

//computes the average position of the particles over avnum simulations
void calcxAv1(vector<double>& vec,const vector< vector<double> >& xAv);

//Markovian case
void calcxAvTheo(vector<double>& vec,const double dt,const double w, const double gamma, const double x0, const double v0,const int nSteps);

void calcxSqrAvTheo(vector<double>& vec,const double dt,const double w, const double gamma,
		    const double mass,const double D,const double x0, const double v0,const int nSteps);

void calcxAvTheo1(vector<double>& vec,const double dt,const double x0,const double v0,const int nSteps);

void calcxAvTheo2(vector<double>& vec,const double dt,const double x0,const double v0,const int nSteps);
//Non-Markovian case
/*double u(double t, const double l1,const double l2, const double l3, const double tau, const double gamma);

double v(double t, const double l1,const double l2, const double l3, const double tau, const double gamma);

void calcxAvTheo1(vector<double>& vec,const double dt,const double l1, const double l2, const double l3,const double tau, const double gamma, const double x0, const double v0,const int nSteps);
*/
// computes the value of kappa (ratio between E_kin and B_eff) for given v0
void kappa(vector<double>& vec,const double v0,const double& mass,const double& B_eff,const int& k);

//gnuplot file für Potenzial
// erstmal nur für zeitunabhängiges Potenzial
void GnuPotential(vector <double>& xVec, vector <double>& potVec, function<double(const double&, const double&)> potExtFunc,
	       const int& n, const double& x0, const double& xn);

//total number of particles to the left of the barrier (->n_a) for comparison with analytical value
void countLeft(vector<double>& vec,const vector< vector<double> >& data, const double& border);

// makes statistic of how many particles are on the right of border over all times
void countRight(vector<double>& vec,vector<double>& vec1, vector<double>& vec2, vector<double>& vec3, const vector< vector<double> >& data, const double& border,const double& Kborder,const bool& Kramers,const int& k,const double& dt);

void calcFlux(vector<double>& vec,const vector<double>& vec1,const vector<double>& rightDist,const double& dt);

//counts number of the particles crossing the potential well for every time step
void fluxLeftRight(vector <double>& vec,vector <double>& vec1,const double border,const vector< vector<double> >& data,const double& dt);

//average of kramersrate
void averageFluxLeftRight(vector<double>& vec, double& averageRate);

// returns the averaged potential energy over all times
void calcPotEnergyAv(vector<double>& vec, const vector< vector<double> >& allX, const vector<double>& tVec, function<double(double, double)> potentialFunc);

//calculates kinetic Energy out of a velocity Vector
void calcKinEnergy(vector<double>& vec, const vector<double>& velocVec, const double& mass);

void calcKinEnergyAv(vector<double>& vec, const vector< vector<double> >& allVeloc, const double& mass);

// calculates potential energy of all particles at one specific time t
// written for prepareEnergyAnimation function
void calcPotEnergy(vector<double>& vec, const vector<double>& xVec, const double& t, function<double(double, double)> potentialFunc);

//calculates totale Energy of all particles at one specific time t
//written for prepareEnergyAnimation function
void calcTotalEnergy(vector<double>& vec, const vector<double>& xVec, const vector<double>& velocVec, const double& t,
			       const double& mass, function<double(double, double)> potentialFunc);

//I(0,t') = Xi(0)*Xi(t')
void calcCorr(vector<double>& vec, const vector<double>& data, const int& startN, const int& endN);

//calculates correlation function I(t)=xi(0)*xi(t) and averages it
void buildCorrAv(vector<double>& vec, const vector< vector<double> >& allNoise, const int& startN, const int& endN);

//----animation functions----

//writes for each time one file for gnuplot with position as x-axis and total energy at y-axis
int prepareEnergyAnimation(const vector< vector<double> >& allX, const vector< vector<double> >& allVeloc, const vector<double>& tVec,
			       const double& mass, function<double(double, double)> potentialFunc, const int& stride, const string& headerLine, const string& folderName);

//same as prepareEnergyAnimation& & writes potential function in file for each time step
int prepareEnergyAnimationTimePot(const vector< vector<double> >& allX, const vector< vector<double> >& allVeloc, const vector<double>& tVec,
			       const double& mass, function<double(double, double)> potentialFunc, const int& potentialN, const double& potBorderL,
				  const double& potBorderR, const int& skipN, const string& headerLine, const string& folderName);

//writes one file for gnuplot of the velocity distribution for each time step
int prepareVda(const vector< vector<double> >& allVeloc, const vector<double>& vDistT, vector<vector<double> >& DistTotal,
	  const double& rangeBeg, const double& rangeEnd, const int& nBins, const int& stride,
	  const string& headerLine, const string& folderName);

//writes one file for gnuplot of the position distribution for each time step
int preparePda(const vector< vector<double> >& allx, const vector<double>& pDistT,vector<vector<double> >& DistTotal,
	  const double& rangeBeg, const double& rangeEnd, const int& nBins, const int& stride,
	  const string& headerLine, const string& folderName);

//dicretized bessel Function K0, used for comparison with numerical results
void besselFunction(vector<double>& vec, const vector<double>& t_vec);

//discretized theoretical fourier transform of the exponential correlation function
void fourier1(vector<double>& vec, const double& samplingRate, const int& sampleN, const int& nSteps);


// analytically expected kinetic energy for exponential correlation function as described in Xus paper
void calcKinEnergyTheo(vector<double>& vec, const vector<double>& tVec, const double& temperature,
		     const double& mass, const double& D, const double& tau, const double& k);

// theoretical spectral density of velocity for colored noise and no potential
void calcSpecVTheo0(vector<double>& vec, const vector<double>& omegaVec, const double& k_b, const double& temperature, const double& mass, const double& D, const double& tau);

// theoretical spectral density of velocity for colored noise and harmonic potential k*x^2
void calcSpecXTheo1(vector<double>& vec, const vector<double>& omegaVec, const double& k_b, const double& temperature, const double& mass, const double& D, const double& tau, const double& potK);

// spectral density of given vector corrVec which is the correlation of some function
void calcSpec(vector<double>& vec, const vector<double>& corrVec, const double& dt);

// <Delta x^2> for white Noise without potential and v0 = 0
void calcXSqrTheo0(vector<double>& vec, const vector<double>& tVec, const double& mass, const double& gamma, const double& k_b, const double& temperature);

// <1/2*mass* v^2> for white Noise without potential
void calcKinEnergyTheo0(vector<double>& vec, const vector<double>& tVec, const double& mass, const double& gamma, const double& k_b, const double& temperature, const double& v0);

//does all calculations after simulating particles, like statistics...
void doAfterMath(const Filenames& filenames,const Foldernames& foldernames, const SimulationOptions& so, const Potential& potential, const Results& results, const NoiseDiss& noiseDiss, const string& currentDateString, KappaSimulation& ksim);
