#include "noiseDiss.h"

#include <stdio.h>
#include <iostream>

#include <tr1/functional>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>

#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>

//#include "randGenXML.hxx"
#include "vectorCalculations.h"

using namespace std;
using namespace std::tr1;
using namespace std::tr1::placeholders;

// no noise
void generateNoNoise(vector<double>& output)
{
  fill(output.begin(), output.end(), 0.0);
}

// generate a sequence of white noise
void generateWhiteNoise (const double& dt, gsl_rng* r,  vector<double>& output)
{ //r: noiseDiss.randGen.at(omp_get_thread_num())
  for (int step = 0; step < output.size(); step++)
  {
    output.at(step) = gsl_ran_gaussian(r, 1.0/sqrt(dt));
  }
}

// generate a sequence of colored noise
void generateColoredNoise(const double& dt, gsl_rng* r, const vector<double>& G, vector<double>& colNoise)
//as done in Xu's code with long prehistory and future correlation
{
  vector<double> whiteNoise(colNoise.size() + 2*G.size() -1, 0.0);
  generateWhiteNoise(dt, r, whiteNoise);
  for (int i = 0; i < colNoise.size(); i++)
  {
    colNoise.at(i) = 0;
    for (int j = 0; j < 2*G.size()-1; j++)
    {
      colNoise.at(i) += 1.0*G.at(abs(j-(int)G.size()+1))*whiteNoise.at(j+i)*dt;
    }
  }
}

// generate exponential correlation function
//I(t)=0.5*D/tau*exp(-t/tau)
vector<double> generateCorrFunc0(const double& t0, const double& dt, const int& nSteps, const double& D, const double& tau)

{
  vector<double> output(nSteps, 0.0);
  double t = t0;
  for (int i = 0; i < nSteps; i++)
  {
    output.at(i) = 0.5*D/tau*exp(-t/tau);
    t += dt;
    //cout << nSteps << endl;
  }
  return output;
}

// generate Gaussian correlation function
// I(t)=1/(a*sqrt(Pi))*e^(-(t/a)^2)
vector<double> generateCorrFunc1(const double& t0, const double& dt, const int& nSteps, const double& D, const double& a)
{
  vector<double> output(nSteps, 0.0);
  double t = t0;
  for (int i = 0; i < nSteps; i++)
  {
    output.at(i) = D/(a*sqrt(M_PI))*exp(-pow(t/a,2));
    t += dt;
  }
  return output;
}

// generate fancy correlation function
vector<double> generateCorrFunc2(const double& t0, const double& dt, const int& nSteps, const double& D, const double& chi)
{
  vector<double> output(nSteps, 0.0);
  double t = t0;
  for (int i = 0; i < nSteps; i++)
  {
    if(t==0)
    {
     output.at(i) = (D/chi); 
    }
    else
    {
    output.at(i) = (D/chi)*sin(abs(t)/chi)/(abs(t)/chi)*exp(-abs(t)/chi);
    }
    t += dt;
  }
  return output;
}

// generate correlation function for massless theory
vector<double> generateCorrFunc3(const double& t0, const double& dt, const int& nSteps, const double& D, const double& alpha,const double& temperature, const double& k, const double& mass)
{
  vector<double> output(nSteps, 0.0);
  double t = t0;
  for (int i = 0; i < nSteps; i++)
  {
    output.at(i) = 2.0*k*temperature/4.0*pow(alpha,2.0)*(1-alpha/sqrt(mass)*t)*exp(-alpha/sqrt(mass)*t);
    t += dt;
  }
  return output;
}

// generate dissipation kernel
//Gamma(t)=0.5/(k*T)*I(t) -- dissipation-fluctuation theorem
vector<double> generateDissipationKernel(const int& nMax, const vector<double>& corrFunction,const int& corrFuncNr, const double& minAlpha,
					 const double& temperature, const double& k, const bool& optimisationBool, const double& dPrecision)
//Gamma(t)=0.5/(k*T)*I(t) -- dissipation-fluctuation theorem
{
  int maxN = max(nMax, (int)corrFunction.size()); //corrFunction.size()=nFourier nMax=nSteps
  vector<double> output(maxN, 0.0);
  int i = 0;
  //cout << maxN << " " << nMax << " " << corrFunction.size() << endl;
  if (optimisationBool)
  {
    if(corrFuncNr==0 || corrFuncNr==1 || corrFuncNr==2)
    {
      do
      {
	output.at(i) = 0.5/(k*temperature)*corrFunction.at(i);
	i++;
      }
      while ((i < maxN) && (output.at(i-1) > dPrecision)); // dPrecision=10**(-3)
      int finalN = i;
      output.resize(finalN);
      if (output.at(finalN-1) > dPrecision)
	{
	  printf("Attention! Last entry of dissipationKernel (%.4f)"
    		  "is greater than minimum wanted (%.4f) using optimisation.\n"
    		  "Make nFourier higher for exacter final values of dissipationKernel.\n", output.at(finalN-1), dPrecision);
	}
    }
  
    
    else
    {
      do
      {
	output.at(i) = 0.5/(k*temperature)*corrFunction.at(i);
	i++;
      }
      while (i <maxN && output.at(i-1)>=-0.3); //  muss an Minimum von Gamma angepasst werden
      do
      {
	output.at(i) = 0.5/(k*temperature)*corrFunction.at(i);
	i++;
      }
      while ((i < maxN) && (fabs(output.at(i-1)) > dPrecision));
      
      int finalN = i;
      output.resize(finalN);
      if (output.at(finalN-1) > dPrecision)
	{
	  printf("Attention! Last entry of dissipationKernel (%.4f)"
    		  "is greater than minimum wanted (%.4f) using optimisation.\n"
    		  "Make nFourier higher for exacter final values of dissipationKernel.\n", output.at(finalN-1), dPrecision);
	}
    }
  }
  
  else
  {
    for(int j=0; j<maxN; j++)
    {
      output.at(j) = 0.5/(k*temperature)*corrFunction.at(j); 
    }
  }
  return output;
}

// calculate the function G(t)
void calcG(const double& dt, const vector<double>& corrFunc,const int& corrFuncNr,
		     const bool& timeSaving, const double& gPrecision, vector<double>& GVec)
{
  int steps = corrFunc.size();
//   cout<< "steps" << " " <<steps << endl;
//   cout << corrFunc.size() << endl;
  vector<double> doubleCorrFunc(2*steps, 0.0);
  mirrorData(doubleCorrFunc, corrFunc,false, false);
  double* G_array =(double*) calloc(2*steps, sizeof(double));
  int nG;
  copy(doubleCorrFunc.begin(), doubleCorrFunc.end(), G_array);
   
   double dw;
   cout << "doubleCorrFunc.size()" << " " << doubleCorrFunc.size() << endl;
   dw = M_PI/dt/(doubleCorrFunc.size()/2+1);
   double omegaMax = 12.0;
   vector<double> specOmega(doubleCorrFunc.size()/2+1);
   vector<double> time(steps);
   vector<double> Mei(doubleCorrFunc.size()/2);
   linspace(0.0, dw, specOmega);
   linspace(0.0, dt, time);
       fstream A,B;
      A.open("fourier.dat", ios::out);
//       B.open("Mei.dat", ios::out);
  // in-place halfcomplex fourier transform
  gsl_fft_real_radix2_transform (G_array, 1, 2*steps);
  
  
  // prepare for inverse fourier transform
  for (int i = 0; i < steps+1; i++)
  {
     G_array[i] = G_array[i]*dt;
     A << specOmega.at(i) << " " << G_array[i] << endl;
     G_array[i] = sqrt(fabs(G_array[i]));
  }
  // set imaginary part to 0, since the correlation function is symmetric with respect to the midpoint
  // can be avoided, when the time interval is symmetric with respect to the origin (shift needed)
  for (int i = steps+1; i < 2*steps; i++)
  {
    G_array[i] = 0;
  }
  // in-place inverse halfcomplex fourier transform
  gsl_fft_halfcomplex_radix2_inverse(G_array, 1, 2*steps);
  int i = 0;
  if (timeSaving)
  {
    if(corrFuncNr==0 || corrFuncNr==1 || corrFuncNr==2)
    {
    do
    {
      G_array[i] = G_array[i]/dt;
     // B << time.at(i) << " " << G_array[i] << endl;
      i++;
    }
    while ((i < steps) && (G_array[i-1] > gPrecision));
    nG = i;
    if (G_array[nG-1] > gPrecision)
     printf("Attention! Last entry of G (%.4f) is greater \n"
    		 "than minimum wanted (%.4f) using optimisation.\n"
    		 "Make nFourier higher for exacter final values of G.\n",
			 	 G_array[nG-1], gPrecision);
    }
    else{
     do
    {
      G_array[i] = G_array[i]/dt;
      i++;
    }
    while ((i < steps) && G_array[i-1]>=-0.06); // muss an Minimum von G(t) angepasst werden
    
    do
    {
      G_array[i] = G_array[i]/dt;
      i++;
    }
    while ((i < steps) && fabs(G_array[i-1]) > gPrecision);
    nG = i; 
    }
  }
  else
  {
    for (int i = 0; i < steps; i++)
    {
       G_array[i] = G_array[i]/dt;
       //B << time.at(i) << " " << G_array[i] << endl;
    }
    nG = steps;
  }
  GVec.resize(nG);
  GVec.assign(G_array,G_array + nG);
}

void initRandGen(const long& seed, vector<gsl_rng *>& randGen){
	  const gsl_rng_type * T;
	  gsl_rng_env_setup();
	  T = gsl_rng_default;
	  //cout << "rand "<< time(NULL) << endl;
	  for (int i = 0; i < randGen.size(); i++)
	  {
		randGen.at(i) = gsl_rng_alloc (T);
		gsl_rng_set (randGen.at(i), seed + i);
	  }
}

NoiseDiss::NoiseDiss(SimulationOptions& so){
	  // select and define local correlation function, 0: exponential, 1: Gaussian
	  switch(so.corrFuncNr)
	  {
		case 0: this->corrFunc = bind(generateCorrFunc0, _1, _2, _3, so.D, so.tau); break;
		case 1: this->corrFunc = bind(generateCorrFunc1, _1, _2, _3, so.D, so.a); break;
		case 2: this->corrFunc = bind(generateCorrFunc2, _1, _2, _3, so.D, so.chi); break;
		case 3: this->corrFunc = bind(generateCorrFunc3, _1, _2, _3, so.D, so.alpha,so.temperature,so.k_b,so.mass); break;
		default: printf("corrFuncNr unknown\n");
		//todo ErrorMessage
		break;
	  }

	  
	  
    fstream Gam,GamOm;
    Gam.open("Gam.dat", ios::out);
    GamOm.open("GamOm.dat", ios::out);
    double dom;
    double om=0.0;
    double tkernel=0.0;
    dom = M_PI/so.dt/(2.0*so.nFourier/2+1);
    vector<double> specOm(2.0*so.nFourier/2+1);
    vector<double> kernel(so.nSteps);
    
    for(int i=0; i<kernel.size(); i++)  //Kernel
    {
	
      kernel.at(i)=1.0/4.0*so.alpha*so.alpha*(1-so.alpha/sqrt(so.mass)*tkernel)*exp(-so.alpha/sqrt(so.mass)*tkernel);
      Gam << tkernel << " " << kernel.at(i) << endl ;
      tkernel += so.dt; 
      
    }
    
     for(int i=0; i<specOm.size(); i++)  //FourierTrafo des Kernels
    {
      specOm.at(i)=so.D/(1.0+so.tau*so.tau*om*om);  //FourierTrafo von Corr0
      //specOm.at(i)=so.D*exp(-(so.a*om/2.0)*(so.a*om/2.0));   //FourierTrafo von Corr1
      //specOm.at(i)=2.0*so.k_b*so.temperature*pow(so.alpha,3.0)*pow(om,2.0)/(sqrt(so.mass)*pow(pow(om,2.0)+pow(so.alpha,2.0)/so.mass,2.0));
      GamOm << om << " " << specOm.at(i) << endl ;
      om += dom; 
      
    }
	  // ----prepare simulation functions and dependent variables----
	  this->tFourierVec.resize(so.nFourier, 0.0);  // all times for nFourier steps       

	  this->corrFVec.resize(so.nFourier, 0.0); // I(t) for t=[0, nFourier*dt]	//Correlation function
	  this->dispKernel.resize(so.nSteps, 0.0); //Gamma(t) for t=[0, nSteps*dt]	//dissipationkernel
	  this->noise.resize(so.nSteps, 0.0);
	  this->allNoise.resize(so.np, vector<double>(so.nSteps, 0.0));
	  linspace(0.0, so.dt, this->tFourierVec);
	  this->corrFVec = this->corrFunc(0.0, so.dt, so.nFourier); // erstes Argument ist t0, zweites "nSteps"                        //erzeugt korrelationsfunktion
      this->dispKernel = generateDissipationKernel(so.nSteps, corrFVec,so.corrFuncNr,so.minAlpha, so.temperature, so.k_b, so.optimisationBool, so.dPrecision);	//erzeugt Dissipationskernel
      //generateDissipationKernel(so.nSteps, corrFVec, so.temperature, so.k_b, so.optimisationBool, so.dPrecision);
	  calcG(so.dt, corrFVec,so.corrFuncNr, so.optimisationBool, so.gPrecision, this->G);
	  this->tGVec.resize(G.size(), 0.0);
	  linspace(so.t0, so.dt, tGVec);

	  switch(so.noiseNr)
	  {
		case 0: noiseFunc = bind(generateNoNoise, _3); break;
		case 1: noiseFunc = bind(generateWhiteNoise, _1, _2, _3); break;
		case 2: noiseFunc = bind(generateColoredNoise, _1, _2, G, _3); break;

		default: printf("noiseNr unknown\n");
		//todo ErrorMessage
		break;
	  }
	  this->seed = time(NULL);
	  this->randGen.resize(so.maxThreadNumber);   //maxThreadNumber = 24
	  initRandGen(this->seed, this->randGen);
}


