#pragma once

#include <stdio.h>
#include <iostream>


#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>


#include "simulationOptions.h"

using namespace std;

class Results{

public:
	vector<double> x0Vec;
	vector<double> v0Vec;
	vector< vector<double> >  allX;
	vector< vector<double> >  allV;

	vector<double> tVec;


	Results(const SimulationOptions&);
	void setInitValues(const SimulationOptions&);
	vector<double> fillThermalized(const SimulationOptions&);
 };
 
 class KappaSimulation{
 public:
	vector<double> transition;
	vector<double> kappa;
	int indexNr;
	vector<vector<double> > xAvVec;
	vector<vector<double> > xSquaredAvVec;
	vector<vector<double> > vAvVec;
	vector<vector<double> > transitionAvVec;
	vector<vector<double> > fluxAvVec;
	vector<vector<double> > fluxPositiveAvVec;
	vector<vector<double> > fluxNegativeAvVec;
	vector<vector<double> > fluxPaperAvVec;
	vector<vector<double> > kinEnergyAvVec;
	vector< vector< vector<double> > > pDistAvVec;
	vector< vector< vector<double> > > vDistAvVec;
	int indexAvNum;
	vector<vector<double> > CorrAvVec;
	vector<vector<double> > noiseAvVec;
	vector<vector<double> > noiseSqrAvVec;
	vector<vector<double> > corrXAv;
	vector<vector<double> > corrVAv;
	
	KappaSimulation(const SimulationOptions&);
   
 };
