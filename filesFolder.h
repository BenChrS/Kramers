#pragma once

#include <stdio.h>
#include <iostream>
#include "simulationOptions.h"

using namespace std;

// appends "(0)" to the filename.
//if this file already exists "(n)" is appended with n being the first number for which the file does not exist yet
//this filename is returned
string existingFile(const string&, const string&);


// writes all_t and all_y in data.size() columns in filename, the headerString is appended for file information
int writeToFile(const vector<double>&, const vector< vector<double> >&,
		 const string&, const string &headerString = "", const string &folderName = "", const bool &checkExistence = false, const int &precision = 5, const string &fileExtension = ".txt");

// writes tVec and data in different columns in filename, the headerString is appended for file information
int writeToFile(const vector<double>&, const vector<double>&, const string&, const string &headerString = "",
		const string &folderName = "", const int &precision = 5, const bool &checkExistence = false, const string &fileExtension = ".txt");

// writes tVec and data in different columns in filename, the headerString is appended for file information
int writeToFile( const string&, const string& , const string &folderName = "",  const bool &checkExistence =false, const string &fileExtension = ".txt");

int writeInitValToFile(const SimulationOptions&,const string &fileName, const string &folderName= "",const bool &checkExistence = false, const string &fileExtension = ".txt");

int writeInitValToFile1(const SimulationOptions&,const string &fileName, const string &folderName= "",const bool &checkExistence = false, const string &fileExtension = ".txt");

class Filenames{
public:
	string infoFile;
	string corrFunc;
	string fourier;
	string corW;
	string dispKernel;
	string whiteNoise;
	string noise;
	string corrAv;
	string colNoise;
	string noiseAv;
	string quadraticVeloc;
	string quadraticVelocAv;
	string kinEnergyTheo;
	string potEnergyAv;
	string kinEnergyAv;
	string totalEnergyAv;
	string pDist;
	string pDist1;
	string xDistTheo;
	string allKinEnDist;
	string x;
	string veloc;
	string xAv;
	string xAvTheo;
	string xAvTotal;
	string velocAv;
	string velocAvTotal;
	string vSqrAvTotal;
	string xSqrAv;
	string xSqrAvTotal;
	string xSqrAvTheo;
	string xSqrTheo;
	string vDist;
	string vDist1;
	string vDistTheo;
	string rightDist;
	string rightDistTotal;//averaged over multiple realisations
	string reinitPart;
	string rightDistKramers;
	string kappaTransprob;
	string potential;
	string allX;
	string allVeloc;
	string specXTheo;
	string specX;
	string specVTheo;
	string specV;
	string corrX;
	string corrV;
	string finalRightDist;
	string finalRightDistExact;
	string flux;
	string fluxPositive;
	string fluxNegative;
	string fluxTotal;
	string fluxPositiveTotal;
	string fluxNegativeTotal;
	string fluxPaperTotal;
	string test;
	string test1;
	string test2;
	string SimOpt;
	string SimOpt1;
	string corrAvTotal;
	string noiseAvTotal;
	string corrXTotal;
	string corrVTotal;
	string velocInitDist;
	string xInitDist;
	string velocInitDist1;
	string xInitDist1;
	Filenames();
 };

class Foldernames{
public:

	string main;
	string ea;
	string vda;
	string pda;

	Foldernames(const string&, const SimulationOptions&);
	void initFolder(const SimulationOptions&);
	int createFolder(const string&);
};
