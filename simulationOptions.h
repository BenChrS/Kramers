#pragma once

#include <stdio.h>
#include <iostream>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <boost/concept_check.hpp>

using namespace std;

class SimulationOptions{
private:


public:

	string commentLine;
	//physical settings
	double k_b; // k_boltzmann constant
	double mass; //mass of particle
	double temperature;
	double D; // see paper
	double tau; // see paper /only important for first correlation function
	double a; // only important for second correlation function
	double chi; //correlation time for third correlation function
	double alpha; // correlation time for the massless theory
	double minAlpha; //Minimum des Kernels
	double aWhite; // eigenvalue for white noise
	double aColoured; // eigenvalue for coloured noise
        double bColoured;
	double cColoured;
	//statistical/program settings
	double t0; //time interval [t0, t1]
	double tEnd;
	double tSettling; // time needed for I(t) to be approximately 0
	double timeSettled; //approximate time particles need to be in equilibrium - only important for kinetic Energy Average - not yet in external call
	int nStepsFactor;
	int nStepsTwo;
	int nSteps; // number of final datapoints (stored), must be devidable by 2
	int npTen;
	int npTwo;
	int np; //number of averaged simulations (number of particles)
	bool optimisationBool;
	double gPrecisionTen;
	double gPrecision; //if opitmisationBool is set to true gPrecision will be the smallest value in G(t)
	double dPrecisionTen;
	double dPrecision; //if optimisationBool is set to true dPrecision will be the smallest value in dispKernel(t)
	int maxThreadNumber; //must be greater than actual used number of threads d
	double maximumAllacationMemory; //program does not run if calculated allocation memory exceeds this value
	bool saveAllDataBool;
	bool avOpt;
	int avNum;
	bool testBool;
	bool paperBool;
	double a_Corr;
	double b_Corr;

	int corrFuncNr;
	int potNr;
	int noiseNr;
	int sdeSolverNr;
	int Kramers;

	//potential Settings
	double potw;//frequency of the harmonic potential
	double potK; //for potential1
	double potDepth; //for potenital2
	double potWidth; //for potential2
	double potA; //for potential 3
	double potB; //for potential 3
	double potC; //for potential 3
	double pota; //for potential 4
	double potMy; //for potential 4
	double potNy; //for potential 4
	double potMax; //for potential 4
	double potStartDepth; //for potential 5
	double potEndDepth; //for potential 5
	double potStartTime; //for potential 5
	double potEndTime; //for potential 5
	double xc; //for potential 7
	double xb; //for potential 7
	double Ub; //for potential 7

	//energy animation (ea) settings
	bool eaBool; // energy Animation?
	bool eaPotTimeDependent; //if true a potential file is generated for each time step of the animation
	double fps; //unbedingt wieder ändern !!
	//int fps; //files per second
	int eaStride; //number of frames skipped for the energy Animation - if set to 1 a file is generated for all time steps
	int eaPotentialN;
	double eaPotBorderL;
	double eaPotBorderR;
	
	//random number generator
	const gsl_rng_type * T;
	gsl_rng * res_rng;
	long seed;

	//initial position and velocity distribution
	int initCondNr;
	double x0;
	double potB1; // height of the inverse harmonic function to be overcome by the particle
	double potB_eff;
	int B_eff;
	double v0;

	//
	double xCutoff; //if absoulute position of particles
					// becomes greater than this value,
					// they will stay at this position.
					// Set to negative value for not making use of this option
	//Kramer
	bool KramersRate;
	double dx;
	double rxborder;
	double lxborder;

	
	//velocity distribution animation (vda) settings
	bool vdaBool; // velocity distribution animation?
	int vdaStride;

	//theoretical velocity distribution settings
	double vDistRangeBeg;
	double vDistRangeEnd;
	int vDistNBins;
	double vDistStartTime;
	
	//position distribution animation (pda) settings
	bool pdaBool; 
	int pdaStride;
	double pDistRangeBeg;
	double pDistRangeEnd;
	int pDistNBins;
	double pDistStartTime;
	
	//kinEnergy distribution animation (eda) settings
	bool edaBool; 
	int edaStride;
	double eDistRangeBeg;
	double eDistRangeEnd;
	int eDistNBins;
	double eDistStartTime;
	
	//transition options
        bool transition;
	int npoints;
	//dependent variables
	double dt;

	  //----prepare simulations----

	double gamma;

	int nSettling; //must be at least as high as I(t) needs to be at approximately 0
	int nFourier; // must be larger than nSettling to avoid boundary
						//  effects of the FFT and be a power of 2 for the FFT


	SimulationOptions();
	void setInitValues(int argc, char *argv[]);
	void setDependentVariables();
	void readInput(int argc, char *argv[]);

 };
