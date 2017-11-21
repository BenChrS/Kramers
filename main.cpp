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
#include <fstream>

#include "filesFolder.h"
#include "simulationOptions.h"
#include "results.h"
#include "potential.h"
#include "sdeSolver.h"
#include "vectorCalculations.h"
#include "afterMath.h"
#include "noiseDiss.h"

using namespace std;
using namespace std::tr1;
using namespace std::tr1::placeholders;

// +++++++++++++++ technical functions +++++++++++++++++++++++++++++++++++
//returns current Date and time as String in the format "DD.MM.YYYY-hh:mm:ss"
const string currentDateTime() 
{
  time_t now = time(0);
  struct tm tstruct;
  char buf[80];
  tstruct = *localtime(&now);
  strftime(buf, sizeof(buf), "%d.%m.%Y-%X", &tstruct);
  return buf;
}

void printStatus(const Results &results, const NoiseDiss &noiseDiss, const SimulationOptions &so){
	printf("results.tVec.size()=%d\n", (int) results.tVec.size());
	printf("corrFVec.size()=%d\n", (int) noiseDiss.corrFVec.size());
	printf("dispKernel.size()=%d\n", (int) noiseDiss.dispKernel.size());
	printf("nSteps=%d\n", so.nSteps);
}

int printWarning(SimulationOptions &so){
	if (so.optimisationBool)
	{
	  printf("optimisationBool set to true. dPrecision and gPrecision are manually set.\n"
			  "Make sure that they don't change the outcome when they're are set to smaller values.\n"
			  "Also pay attention to correlation Functions if they rise again back from 0.\n\n");
	}

	printf("program allocates approximately %.4f GBytes!\n", (3 * so.nSteps * (so.np+10) * sizeof(double))/1073741824.0);
	if ( (3 * so.nSteps * (so.np+10) * sizeof(double))/1073741824.0 > so.maximumAllacationMemory)
	{
	  printf("Simulation exceeds maximum allocation memory of %.4f GBytes!\n"
				"This value is manually set so set the variable maximumAllacationMemory to a higher\n"
				"value if you think there is more memory available.\n", so.maximumAllacationMemory);
	  return -1;
	}
	return 0;
}

int main(int argc, char** argv) {  //argc: 1, argv[0]: ./colnoise
	cout << currentDateTime() << endl;
  
	string currentDateString = currentDateTime();
	stringstream headerStringStream;
	headerStringStream << "# " << currentDateString;
	string headerString = headerStringStream.str();

	//----Variable settings----
	const string versionNr = "5.0";
	
// 	if(argc != 2)
// 	{
// 	 cout << "inputFile is missing" << endl;
// 	 return -1;
// 	}
	SimulationOptions so;
//  	so.inputfile=argv[1];
        Filenames filenames;
//	Foldernames foldernames("test", so);
  	Foldernames foldernames(currentDateString, so);

	if (printWarning(so) == -1){
		return -1;
	}
	Results results(so); //x, v, t
	KappaSimulation ksim(so);
	NoiseDiss noiseDiss(so); // noise and dissipation functions
	
	Potential potential(so); //potential and force
        SdeSolver sdeSolver(so, noiseDiss);
	writeInitValToFile(so,filenames.SimOpt,foldernames.main); //initial values
	writeInitValToFile1(so,filenames.SimOpt1,foldernames.main); //initial values (reduced form)

	printStatus(results, noiseDiss, so);

	writeToFile(results.tVec, noiseDiss.corrFVec, filenames.corrFunc, headerString, foldernames.main);
	writeToFile(noiseDiss.tGVec, noiseDiss.G, filenames.fourier, headerString, foldernames.main);
	writeToFile(results.tVec, noiseDiss.dispKernel, filenames.dispKernel, headerString, foldernames.main);
         
// 	  fstream velocDist,posDist;
//           velocDist.open("velocDist.dat", ios::out);
//           for(int k=0; k < results.v0Vec.size(); k++)
//             {
// 	    velocDist << k << " " << results.v0Vec.at(k) << endl; 
// 	    }
// 	  velocDist.close();
// 	  posDist.open("posDist.dat", ios::out);
// 	  for(int k=0; k < results.x0Vec.size(); k++)
//             {
// 	    posDist << k << " " << results.x0Vec.at(k) << endl; 
// 	    }
// 	  posDist.close();
	  
	

	//----actual simulations----

if(so.avOpt==true && so.testBool==false)  // wichtig: BEACHTE so.avOpt==true !!
{
  cout << "first loop" << endl;
  int j;
for(j=0;j<so.avNum;j++)
{
  printf("%.0f %%\n", (double)(j*100)/so.avNum);
ksim.indexAvNum = j;
int i;
	#pragma omp parallel for default(shared) private(i)// Parallelization of the particle loop with OpenMP
	for (i = 0; i < so.np; i++)
	{
	    //printf("%d \n", omp_get_thread_num());
	      //usr information and percentage done
		if ((i == 0) && (!omp_get_thread_num()))
		  printf("working together with %d threads\n", omp_get_num_threads());
		//if (i % (int)ceil((max(double(so.np)/omp_get_num_threads()/100.0, 1.0))) == 0 && !omp_get_thread_num())
		 // printf("%.0f %%\n", (double)(i*omp_get_num_threads()*100)/so.np);
		
		//cout << "omp " << omp_get_thread_num() << endl;
		//calculate noise
		noiseDiss.noiseFunc(so.dt, noiseDiss.randGen.at(omp_get_thread_num()), noiseDiss.allNoise.at(i));
	        //calculate all steps
		sdeSolver.sdeSolverFunc(results.tVec, results.x0Vec.at(i), results.v0Vec.at(i), so.mass,noiseDiss.allNoise.at(i), results.allX.at(i), results.allV.at(i), potential.forceFunc);
	}  
//----saveAllData----
	if (so.saveAllDataBool)
	{
		writeToFile(results.tVec, results.allX, filenames.allX, headerString, foldernames.main);
		writeToFile(results.tVec, results.allX, filenames.allVeloc, headerString, foldernames.main);
	}
doAfterMath(filenames,foldernames,so, potential, results, noiseDiss, headerString,ksim);
}
}
//===============================================================================================
else
{
    if(so.potNr==1 &&  so.initCondNr==0 && so.transition==true)
	{
	cout << "transition loop" << endl;
	double dE_kin= 3.0/so.npoints;
	double E_kin = 0.0;
	double v0 = 0.0;
	int n;
    for(n=0;n<=so.npoints;n++)
	{
	  printf("%.0f %%\n", (double)(n*100)/(so.npoints+1));
		ksim.indexNr = n;
		E_kin += dE_kin;
		v0 = sqrt((2*(E_kin))/so.mass);
	int i;
	#pragma omp parallel for default(shared) private(i) // Parallelization of the particle loop with OpenMP
	for (i = 0; i < so.np; i++)
	{
		//usr information and percentage done
		if ((i == 0) && (!omp_get_thread_num()))
		  printf("working together with %d threads\n", omp_get_num_threads());
		//if (i % (int)ceil((max(double(so.np)/omp_get_num_threads()/100.0, 1.0))) == 0 && !omp_get_thread_num())
		  //printf("%.0f %%\n", (double)(i*omp_get_num_threads()*100)/so.np);

		//calculate noise
		noiseDiss.noiseFunc(so.dt, noiseDiss.randGen.at(omp_get_thread_num()), noiseDiss.allNoise.at(i));

		// calculate all steps
		sdeSolver.sdeSolverFunc(results.tVec, results.x0Vec.at(i), results.v0Vec.at(i), so.mass,noiseDiss.allNoise.at(i), results.allX.at(i), results.allV.at(i), potential.forceFunc);
	}
	//----saveAllData----
	if (so.saveAllDataBool)
	{
		writeToFile(results.tVec, results.allX, filenames.allX, headerString, foldernames.main);
		writeToFile(results.tVec, results.allX, filenames.allVeloc, headerString, foldernames.main);
	}

	doAfterMath(filenames,foldernames,so, potential, results, noiseDiss, headerString, ksim);
	fill(results.v0Vec.begin(), results.v0Vec.end(), v0);
	}
	writeToFile(ksim.kappa,ksim.transition,filenames.kappaTransprob, headerString,foldernames.main);
	
	}


writeToFile(results.tVec, results.v0Vec,filenames.test, headerString,foldernames.main);
}
cout << currentDateTime() << endl;
return 0;
}

