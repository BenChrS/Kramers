#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <fstream>
#include <sstream>

#include "filesFolder.h"
#include "simulationOptions.h"

using namespace std;

string existingFile(const string &fileName, const string &fileExtension)
// appends "(0)" to the filename.
//if this file already exists "(n)" is appended with n being the first number for which the file does not exist yet
//this filename is returned
{
  string fileNameOut = "";
  string fileNameComplete = "";
  fileNameOut.append(fileName);
  fileNameOut.append("(0)");
  fileNameComplete = fileNameOut;
  fileNameComplete.append(fileExtension);
  int i = 1;
  while (boost::filesystem::exists(fileNameComplete))
  {
    fileNameOut = "";
    fileNameOut.append(fileName);
    fileNameOut.append("(");
    fileNameOut.append(boost::lexical_cast<string>( i ));
    fileNameOut.append(")");
    fileNameComplete = fileNameOut;
    fileNameComplete.append(fileExtension);
    i++;
  }
  return fileNameOut;
}

int writeToFile(const vector<double>& tVec, const vector< vector<double> >& data,
		 const string &fileName, const string &headerString, const string &folderName, const bool &checkExistence, const int &precision, const string &fileExtension)
// writes all_t and all_y in data.size() columns in filename, the headerString is appended for file information
{
  string fileNameLong = "";
  if (folderName != "")
  {
      fileNameLong.append(folderName);
      fileNameLong.append("/");
  }
  fileNameLong.append(fileName);
  if (checkExistence)
  {
    fileNameLong = existingFile(fileNameLong, fileExtension);
  }
  fileNameLong.append(fileExtension);
  if (boost::filesystem::exists(fileNameLong))
  {
    printf("%s already exists. Creating new.\n", fileNameLong.c_str());
  }
  ofstream myfile (fileNameLong.c_str());
    if (myfile.is_open())
    {
      myfile << headerString << endl;
      for (int i = 0; i < tVec.size(); i++)
      {
	myfile << setprecision(precision) << tVec.at(i) << ' ';
	for (int j = 0; j < data.size(); j++)
	{
	  myfile <<  setprecision(precision) << data.at(j).at(i) << ' ';
	}
	myfile << endl;
      }
    }
    else printf("Unable to open file: %s\n", fileNameLong.c_str());
    myfile.close();
    return 1;
}

int writeToFile(const vector<double>& tVec, const vector<double>& data, const string &fileName, const string &headerString,
		const string &folderName, const int &precision, const bool &checkExistence, const string &fileExtension)
// writes tVec and data in different columns in filename, the headerString is appended for file information
{
  string fileNameLong = "";
  if (folderName != "")
  {
      fileNameLong.append(folderName);
      fileNameLong.append("/");
  }
  fileNameLong.append(fileName);
  if (checkExistence)
  {
    fileNameLong = existingFile(fileNameLong, fileExtension);
  }
  fileNameLong.append(fileExtension);
  if (boost::filesystem::exists(fileNameLong))
  {
    printf("%s already exists. Creating new.\n", fileNameLong.c_str());
  }
  int nLines = min(tVec.size(), data.size());

  ofstream myfile (fileNameLong.c_str());
    if (myfile.is_open())
    {
      myfile << headerString << endl;
      for (int i = 0; i < nLines; i++)
      {
	myfile << setprecision(precision) << tVec.at(i) << ' ';
	myfile <<  setprecision(precision) << data.at(i) << ' ';
	myfile << endl;
      }
    }
    else printf("Unable to open file: %s\n", fileNameLong.c_str());
    myfile.close();
    return 1;
}

int writeToFile( const string &data, const string &fileName, const string &folderName,
		 const bool &checkExistence, const string &fileExtension)
// writes tVec and data in different columns in filename, the headerString is appended for file information
{
  string fileNameLong = "";
  if (folderName != "")
  {
      fileNameLong.append(folderName);
      fileNameLong.append("/");
  }
  fileNameLong.append(fileName);
  if (checkExistence)
  {
    fileNameLong = existingFile(fileNameLong, fileExtension);
  }
  fileNameLong.append(fileExtension);
  if (boost::filesystem::exists(fileNameLong))
  {
   printf("%s already exists. Creating new.\n", fileNameLong.c_str());
  }

  ofstream myfile (fileNameLong.c_str());
    if (myfile.is_open())
    {
      myfile << data << endl;
    }
    else printf("Unable to open file: %s\n", fileNameLong.c_str());
    myfile.close();
    return 1;
}

int writeInitValToFile(const SimulationOptions& so,const string &fileName, const string &folderName,const bool &checkExistence, const string &fileExtension)
{
string fileNameLong = "";
  if (folderName != "")
  {
      fileNameLong.append(folderName);
      fileNameLong.append("/");
  }
  fileNameLong.append(fileName);
  if (checkExistence)
  {
    fileNameLong = existingFile(fileNameLong, fileExtension);
  }
  fileNameLong.append(fileExtension);
  if (boost::filesystem::exists(fileNameLong))
  {
   printf("%s already exists. Creating new.\n", fileNameLong.c_str());
  }
    ofstream myfile (fileNameLong.c_str());
    if (myfile.is_open())
    {
    myfile      <<"physical settings:"<<endl
                <<"k_b"<<" = "<<so.k_b<<endl
                <<"mass"<<" = "<<so.mass<<endl
                <<"temperature"<<" = "<<so.temperature<<endl
                <<"D"<<" = "<<so.D<<endl
                <<"tau"<<" = "<<so.tau<<endl
                <<"a"<<" = "<<so.a<<endl
                <<"chi"<<" = "<<so.chi<<endl
                <<"alpha"<<" = " << so.alpha << endl
                <<endl
                <<"statistical/program settings:"<<endl
                <<"t0"<<" = "<<so.t0<<endl
                <<"t1"<<" = "<<so.t1<<endl
                <<"tSettling"<<" = "<<so.tSettling<<endl
                <<"timeSettled"<<" = "<<so.timeSettled <<endl
                <<"nStepsFactor"<<" = "<<so.nStepsFactor <<endl
                <<"nStepsTwo"<<" = "<<so.nStepsTwo <<endl
                <<"nSteps"<<" = "<< so.nSteps<<endl
                <<"npTen"<<" = "<< so.npTen<<endl
                <<"npTwo"<<" = "<<so.npTwo <<endl
                <<"np"<<" = "<<so.np <<endl
                <<"optimisationsBool"<<" = "<< so.optimisationBool<<endl
                <<"gPrecisionTen"<<" = "<<so.gPrecisionTen<<endl
                <<"gPrecision"<<" = "<<so.gPrecision<<endl
                <<"dPrecisionTen"<<" = "<<so.dPrecisionTen<<endl
                <<"dPrecision"<<" = "<<so.dPrecision<<endl
                <<"maxThreadNumber"<<" = "<<so.maxThreadNumber<<endl
                <<"maximumAllocationMemory"<<" = "<<so.maximumAllacationMemory<<endl
                <<"saveAllDataBool"<<" = "<<so.saveAllDataBool<<endl
		<<"avOpt"<<" = "<<so.avOpt<<endl
		<<"avNum"<<" = "<<so.avNum<<endl
                <<endl
                <<"corrFuncNr"<<" = "<<so.corrFuncNr<<endl
                <<"potNr"<<" = "<<so.potNr<<endl
                <<"noiseNr"<<" = "<<so.noiseNr<<endl
                <<"sdeSolverNr"<<" = "<<so.sdeSolverNr<<endl
                <<"Kramers"<<" = "<<so.Kramers<<endl
                <<"initCondNr"<<" = "<<so.initCondNr<<endl
                <<endl
                <<"potential settings:"<<endl
                <<"potw"<<" = "<<so.potw<<endl
                <<"potK"<<" = "<<so.potK<<endl
                <<"potDepth"<<" = "<<so.potDepth<<endl
                <<"potWidth"<<" = "<<so.potWidth<<endl
                 <<endl
                <<"potB"<<" = "<<so.potB<<endl
                <<"potC"<<" = "<<so.potC<<endl
                <<"pota"<<" = "<<so.pota<<endl
                <<"potMy"<<" = "<<so.potMy<<endl
                <<"potNy"<<" = "<<so.potNy<<endl
                <<"potStartDepth"<<" = "<<so.potStartDepth<<endl
                <<"potEndDepth"<<" = "<<so.potEndDepth<<endl
                <<"potStartTime"<<" = "<<so.potStartTime<<endl
                <<"potEndTime"<<" = "<<so.potEndTime<<endl
                <<endl
                <<"xc" <<" = " << so.xc << endl
                <<"xb" <<" = " << so.xb << endl
                <<"Ub" <<" = " << so.Ub << endl
                << endl
                <<"energy animation (ea) settings:"<<endl
                <<"eaBool"<<" = "<<so.eaBool<<endl
                <<"eaPotTimeDependant"<<" = "<<so.eaPotTimeDependent<<endl
                <<"fps"<<" = "<<so.fps<<endl
                <<"eaStride"<<" = "<<so.eaStride<<endl
                <<"eaPotentialN"<<" = "<<so.eaPotentialN<<endl
                <<"eaPotBorderL"<<" = "<<so.eaPotBorderL<<endl
                <<"eaPotBorderR"<<" = "<<so.eaPotBorderR<<endl
                <<endl
                <<"random number generator:"<<endl
                <<"T"<<" = "<<so.T<<endl
                <<"rng"<<" = "<<so.res_rng<<endl
                <<"seed"<<" = "<<so.seed<<endl
                <<endl
                <<"Initial conditions:"<<endl
                <<"x0"<<" = "<<so.x0<<endl
                <<"v0"<<" = "<<so.v0<<endl
                <<"potB1"<<" = "<<so.potB1<<endl
                <<"aWhite"<<" = "<<so.aWhite<<endl
                <<"aColoured"<<" = "<<so.aColoured<<endl
                <<"potB_eff"<<" = "<<so.potB_eff<<endl
                <<"rxborder"<<" = "<<so.rxborder<<endl
                <<"lxborder"<<" = "<<so.lxborder<<endl
                <<endl
                <<"transition options:"<<endl
                <<"transition"<<" = "<<so.transition<<endl
                <<"npoints"<<" = "<<so.npoints<<endl
                <<endl
                <<"Position distribution animation(pda):"<<endl
                <<"pdaBool"<<" = "<<so.pdaBool<<endl
                <<"pdaStride"<<" = "<<so.pdaStride<<endl
                <<"pDistRangeBeg"<<" = "<<so.pDistRangeBeg<<endl
                <<"pDistRangeEnd"<<" = "<<so.pDistRangeEnd<<endl
                <<"pDistNBins"<<" = "<<so.pDistNBins<<endl
                <<"pDistStartTime"<<" = "<<so.pDistStartTime<<endl
                <<endl
                <<"Velocity distribution animation(pda):"<<endl
                <<"vdaBool"<<" = "<<so.vdaBool<<endl
                <<"vdaStride"<<" = "<<so.vdaStride<<endl
                <<"vDistRangeBeg"<<" = "<<so.vDistRangeBeg<<endl
                <<"vDistRangeEnd"<<" = "<<so.vDistRangeEnd<<endl
                <<"vDistNBins"<<" = "<<so.vDistNBins<<endl
                <<"vDistStartTime"<<" = "<<so.vDistStartTime<<endl
                <<endl
                <<"dt"<<" = "<<so.dt<<endl
                <<"gamma"<<" = "<<so.gamma<<endl
                <<"nSettling"<<" = "<<so.nSettling<<endl
                <<"nFourier"<<" = "<<so.nFourier<<endl;

    }
}

Filenames::Filenames(){
	this->colNoise = "colNoise";
	this->corW = "corW";
	this->corrAv = "corrAv";
	this->corrFunc =  "corrFunc";
	this->corrV = "corrV";
	this->corrX = "corrX";
	this->dispKernel = "dispKernel";
	this->finalRightDist = "finalRightDist";
	this->finalRightDistExact = "finalRightDistExact";
	this->fourier = "fourier";
	this->infoFile = "infoFile";
	this->kinEnergyAv = "kinEnergyAv";
	this->kinEnergyTheo = "kinEnergyTheo";
	this->noise = "noise";
	this->noiseAv = "noiseAv";
	this->potEnergyAv = "potEnergyAv";
	this->potential = "potential";
	this->quadraticVeloc = "quadraticVeloc";
	this->quadraticVelocAv = "quadraticVelocAv";
	this->rightDist = "rightDist";
	this->rightDistTotal = "rightDistTotal";
    this->reinitPart = "reinitPart";
    this->rightDistKramers = "rightDistKramers";
	this->kappaTransprob = "kappaTransprob";
	this->allVeloc = "allVeloc";
	this->allX = "allX";
	this->specV = "specV";
	this->specVTheo = "specVTheo";
	this->specX = "specX";
	this->specXTheo = "specXTheo";
	this->totalEnergyAv = "totalEnergyAv";
	this->vDist = "vDist";
	this->vDist1 = "vDist1";
	this->vDistTheo = "vDistTheo";
	this->veloc = "veloc";
	this->velocAv = "velocAv";
	this->velocAvTotal = "velocAvTotal";
	this->whiteNoise = "whitNoise";
	this->pDist = "pDist";
	this->pDist1 = "pDist1";
	this->x = "x";
	this->xAv = "xAv";
	this->xAvTheo = "xAvTheo";
	this->xAvTotal = "xAvTotal";
	this->xSqrAv = "xSqrAv";
	this->xSqrAvTotal = "xSqrAvTotal";
	this->xSqrAvTheo = "xSqrAvTheo";
	this->xSqrTheo = "xSqrTheo";
	this->flux = "flux";
	this->fluxPositive = "fluxPositive";
	this->fluxNegative ="fluxNegative";
	this->fluxTotal = "fluxTotal";
	this->fluxPositiveTotal = "fluxPositiveTotal";
	this->fluxNegativeTotal = "fluxNegativeTotal";
	this->fluxPaperTotal = "fluxPaperTotal";
	this->test = "test";
	this->test1 = " test1";
	this->test2 = " test2";
	this->SimOpt= "SimOpt";
	this->corrAvTotal="corrAvTotal";
	this->noiseAvTotal="noiseAvTotal";
	this->corrXTotal="corrXTotal";
	this->corrVTotal="corrVTotal";

}

Foldernames::Foldernames(const string &mainFolderName, const SimulationOptions &so){
	this->main = mainFolderName;
	this->ea = this->main;
	this->vda = this->main;
	this->pda = this->main;
	this->ea.append("/ea");
	this->vda.append("/vda");
	this->pda.append("/pda");

	initFolder(so);
}

void Foldernames::initFolder(const SimulationOptions &so){
	this->createFolder(this->main);
	  if(so.eaBool)
		this->createFolder(this->ea);
	  if (so.vdaBool)
		this->createFolder(this->vda);
	  if (so.pdaBool)
		this->createFolder(this->pda);
	  
}

int Foldernames::createFolder(const string &folderName)
{
  boost::filesystem::path dir(folderName);
  if (boost::filesystem::create_directory(dir))
  {
    return 0;
  }
  else
  {
    printf("could not create directory or already exists: %s\n", folderName.c_str());
    return 1;
  }
}

