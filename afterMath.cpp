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
#include "afterMath.h"
//#include "kappasimulation.h"

using namespace std;
using namespace std::tr1;
using namespace std::tr1::placeholders;

#include "vectorCalculations.h"
#include "afterMath.h"

// probability density Distribution of velocity for one specific time (same for position)
// rangeBeg: left boundary of first bin
// rangeEnd: right boundary of last bin
// nBins: number of bins
void calcVDistribution(vector<double>& vec,vector<double>& vec1, const vector<double>& data, const double& rangeBeg, const double& rangeEnd, const int& nBins)
{
//   gsl_histogram * h = gsl_histogram_alloc (nBins);
//   gsl_histogram_set_ranges_uniform (h, rangeBeg, rangeEnd);
//   for (int i = 0; i < data.size(); i++)
//   {
//     gsl_histogram_increment (h, data.at(i));
//   }
//   cout << gsl_histogram_mean (h) << endl;
	vec.resize(nBins, 0.0);
	vec1.resize(nBins, 0.0);
  double range = rangeEnd-rangeBeg;
  int j;
  for (int i = 0; i < data.size(); i++)
  {
    if ((data.at(i) >= rangeBeg) && (data.at(i) < rangeEnd))
    {
      j = floor((data.at(i)-rangeBeg)*nBins/range);
      //cout << j << endl;
      vec.at(j) += 1;   
    }
  }
  vec1=vec;
  int k=0;
  for(int i = 0; i<vec.size();i++)
  {
  k +=vec.at(i);  //number of particles in the range 
  }
  //cout << k << endl;
  /*fstream g;
  g.open("Dist.dat",ios::out);
  for(int j=0; j<vec.size();j++)
  {
   g << vec.at(j) << endl; 
  }*/
  double delta = (rangeEnd-rangeBeg)/nBins;
  transform(vec.begin(), vec.end(), vec.begin(), bind1st(multiplies<double>(), 1.0/(data.size()*delta)));
  //cout << data.size() << delta << data.size()*delta << endl;
  //normalised distribution
  transform(vec1.begin(), vec1.end(), vec1.begin(), bind1st(multiplies<double>(), 1.0/k));
  
}
//averaged probability density Distribution of velocity over range from rowBeg to rowEnd
//averaging over different times?
void calcvDistributionMean(vector<double>& vec,vector<double>& vec1, const vector< vector<double> >& data, const int& rowBeg, const int& rowEnd, const double& rangeBeg,const double& rangeEnd, const int& nBins)
{
  vector< vector<double> > temp(rowEnd-rowBeg, vector<double>(nBins, 0.0));
  vector< vector<double> > temp1(rowEnd-rowBeg, vector<double>(nBins, 0.0));
  for (int i = 0; i < rowEnd-rowBeg; i++)
  {
	vector<double> row;
	getRow(row, data, i + rowBeg);
    calcVDistribution(temp.at(i),temp1.at(i), row, rangeBeg, rangeEnd, nBins);
  }
  buildAverage(vec, temp);
  buildAverage(vec1,temp1);
}

// returns the Maxwell Boltzmann probability density distribution for one dimension v in the range [rangeBeg, rangeEnd]
void calcVDistTheo(vector<double>& vec, const double& temperature, const double& mass, const double& k, const double& rangeBeg, const double& rangeEnd, const int& nBins)
{
  vec.resize(nBins, 0.0);
  double v;
  double delta;
  delta = (rangeEnd-rangeBeg)/nBins;
  v = rangeBeg + 0.5*delta;
  for (int i = 0; i < nBins; i++)
  {
    vec.at(i) = sqrt(mass/(2.0*M_PI*k*temperature))*exp(-mass*v*v/(2.0*k*temperature));
    v += delta;
  }
}

// returns the Maxwell Boltzmann probability density distribution for one dimension x in the range [rangeBeg, rangeEnd]
void calcXDistTheo(vector<double>& vec, const double& temperature, const double& mass,const double& omega, const double& k, const double& rangeBeg, const double& rangeEnd, const int& nBins)
{
  vec.resize(nBins, 0.0);
  double x;
  double delta;
  delta = (rangeEnd-rangeBeg)/nBins;
  x = rangeBeg + 0.5*delta;
  for (int i = 0; i < nBins; i++)
  {
    vec.at(i) = sqrt((mass*omega*omega)/(2.0*M_PI*k*temperature))*exp(-mass*omega*omega*x*x/(2.0*k*temperature));
    x += delta;
  }
}

// computes the average position of the particles
void calcxAv(vector<double>& vec,const vector< vector<double> >& data)
{
 vec.resize(data.at(0).size(),0.0);
 for (int i=0;i <data.at(0).size();i++)
 {
   for(int j=0;j <data.size();j++)
   {
     vec.at(i)+=data.at(j).at(i);
   }
 }
transform(vec.begin(), vec.end(), vec.begin(), bind1st(multiplies<double>(), 1.0/data.size())); 
}

//computes the average position of the particles over avnum simulations

void calcxAv1(vector<double>& vec,const vector< vector<double> >& xAv)
{
vec.resize(xAv.at(0).size(),0.0);
 for (int i=0;i <xAv.at(0).size();i++)
 {
   for(int j=0;j <xAv.size();j++)
   {
     vec.at(i)+=xAv.at(j).at(i);
   }
 }
 transform(vec.begin(), vec.end(), vec.begin(), bind1st(multiplies<double>(), 1.0/xAv.size()));  
}
//Markovian case
void calcxAvTheo(vector<double>& vec,const double dt,const double w, const double gamma, const double x0, const double v0,const int nSteps)
{
 vec.resize(nSteps,0.0);
 double l1,l2,a;
 double t=0.0;
 a=sqrt(pow(gamma,2.0)+4*pow(w,2.0));
 l1=-gamma/2+a/2;
 l2=-gamma/2-a/2;
 
 for(int i=0; i<nSteps;i++)
 {
   vec.at(i)=1/(2*a)*((gamma*(exp(l1*t)-exp(l2*t))+a*(exp(l1*t)+exp(l2*t)))*x0-2*v0*(exp(l2*t)-exp(l1*t)));
   t += dt;
 }
}
void calcxSqrAvTheo(vector<double>& vec,const double dt,const double w, const double gamma,
		    const double mass,const double D,const double x0, const double v0,const int nSteps)
{
 vec.resize(nSteps,0.0);
 double l1,l2,a;
 double t=0.0;
 a=sqrt(pow(gamma,2.0)+4*pow(w,2.0));
 l1=-gamma/2+a/2;
 l2=-gamma/2-a/2;
 for(int i=0; i<nSteps;i++)
 {
   vec.at(i)=pow(1/(2*a)*((gamma*(exp(l1*t)-exp(l2*t))+a*(exp(l1*t)+exp(l2*t)))*x0-2*v0*(exp(l2*t)-exp(l1*t))),2.0)-
             D/(4*pow(a,2.0)*pow(mass,2.0))*((1-exp(2*l2*t))/l2+(1-exp(2*l1*t))/l1-4*(1-exp((l1+l2)*t))/(l1+l2));
   t += dt;
 }
}

//Non-Markovian case,w=0.6,gamma=0.6,x0=-1,v0=0.548286
/*void calcxAvTheo1(vector<double>& vec,const double dt,const int nSteps)
{
 vec.resize(nSteps,0.0);
 double t;
 for(int i=0; i<nSteps;i++)
 {
   vec.at(i)=((0.485252+2*0.09115857)*sin(0.394001*t)-(0.288551+2*0.251535)*cos(0.394001*t))*exp(-0.48214*t)-(0.71145-0.503069)*exp(0.464281*t);
   t += dt;
 }
}*/
void calcxAvTheo1(vector<double>& vec,const double dt,const double x0,const double v0,const int nSteps)
{
 vec.resize(nSteps,0.0);
 double t;
 for(int i=0; i<nSteps;i++)
 {
   vec.at(i)=(2*(-0.242625*x0+0.16704*v0)*sin(0.394002*t)+2*(0.144276*x0-0.458767*v0)*cos(0.394002*t))*exp(-0.48214*t)+(0.71145*x0+0.917535*v0)*exp(0.464281*t);
   t += dt;
 }
}

void calcxAvTheo2(vector<double>& vec,const double dt,const double x0,const double v0,const int nSteps)
{
 vec.resize(nSteps,0.0);
 double t;
 for(int i=0; i<nSteps;i++)
 {
   vec.at(i)=(2*(-1.1602*x0+1.12114*v0)*sin(0.0847125*t)+2*(0.184452*x0-0.458152*v0)*cos(0.0847125*t))*exp(-0.361347*t)+(0.631096*x0+0.916304*v0)*exp(0.522694*t);
   t += dt;
 }
}

// computes the value of kappa (ratio between E_kin and potB1) for given v0 //RAUSHAUEN!!
void comp_kappa(vector<double>& kappa,const double v0,const double& mass,const double& potB1,const int& k)
{
double E_kin = mass*pow(v0,2)/2;
kappa.at(k) = E_kin/potB1;
}


//gnuplot file für Potenzial
// erstmal nur für zeitunabhängiges Potenzial
void GnuPotential(vector <double>& xVec, vector <double>& potVec, function<double(const double&, const double&)> potExtFunc,
	       const int& n, const double& x0, const double& xn)
{
	int t = 0;
	double deltaX= xn - x0;
	double x = x0;
	double dx = deltaX/n;
	
	for(int l=0; l<n ; l++)
	{
	 xVec.at(l) = x;
	 potVec.at(l)= potExtFunc(t,xVec.at(l));
	 x += dx;
	}
}

//total number of particles to the left of the barrier (->n_a) for comparison with analytical value
void countLeft(vector<double>& vec,const vector< vector<double> >& data, const double& border)
{
vec.resize(data.at(0).size(),0.0);
for (int i = 0; i < data.at(0).size(); i++)    //i Zeitpunkt
      {
        for (int j = 0; j < data.size(); j ++) //j-tes Teilchen
        {
          if (data.at(j).at(i) < border)	//befindet sich j-tes Teilchen zum Zeitpunkt i links der Barriere?
          {
            vec.at(i) += 1;
	  }
	}
      }
}

// makes statistic of how many particles are on the right of border over all times
void countRight(vector<double>& vec,vector<double>& vec1, vector<double>& transprob, vector<double>& vec3,const vector< vector<double> >& data, const double& border,const double& Kborder,const bool& Kramers,const int& k,const double& dt)
{
 //countRight(rightDist,reinitPart,ksim.transition,rightDistKramers, results.allX, so.xb,so.rxborder,so.KramersRate,i);
  //vec1 :: total number of particles one the right side, vec:: relative number of particles in the right side
    //if Kramers=true :: vec1 reinitialized particles are taken into account
  //vec3: um backscattering und Neu-Initialiserung bereinigte Teilchenzahl auf der rechten Seite
  vec.resize(data.at(0).size(),0.0);
  vec1.resize(data.at(0).size(),0.0);
  vec3.resize(data.at(0).size(),0.0);
  double t=0.0;
 
//   vector<double> leftParticles;
//   vector<double> fluxLeft;
//   vector<double> fluxRight;
  
//   leftParticles.resize(data.at(0).size(),0.0);  
  
//   fluxLeft.resize(data.at(0).size(),0.0);//fluxLeft behinhaltet schon neu-Initialiserung plus backscattering
//   fluxRight.resize(data.at(0).size(),0.0);
//   fluxLeftRight(fluxRight,fluxLeft,border,data,dt); //fluxLeft:  Anzahl der Teilchen, die von i-1 nach i auf linke Seite gewandert sind
  
//   countLeft(leftParticles,data,border);  //Zähle Teilchen auf linker Seite für jeden Zeitschritt i -> für Normierung ?!?
  
  if(Kramers)
  {
    fstream countRight,kramersClean;
    countRight.open("kramersFlux.dat",ios::out);
    kramersClean.open("clean.dat",ios::out);
    for (int i = 0; i < data.at(0).size(); i++)
    {
        for (int j = 0; j < data.size(); j ++)
        {
          if (data.at(j).at(i) > border)
          {
            vec.at(i) += 1;
            if(data.at(j).at(i) >=Kborder)
            {
                vec1.at(i) += 1;   //number of reinitialized particles per time interval
            }
          }
        }
								      //akkumuliere Neu-Initialiserungen!!
    }


//     for(int i=0;i<vec.size();i++)
//     {
// 	for(int j=0;j<i+1;j++)
// 	{
// 	  vec.at(i) += vec1.at(j);		//Neu-Initialisierung werden für jeden Zeitschritt akkumuliert
// 	}
//     }
    for(int i=1;i<vec.size();i++)
    {
	for(int j=0;j<i;j++)
	{
	  vec.at(i) += vec1.at(j);		//Neu-Initialisierung werden für jeden Zeitschritt akkumuliert
	}					//wichtig: für den Zeitschritt i müssen von vec müssen nur die vorherigen Zeitschritte
    }						//j<i aufaddiert werden 
  }
  else
  {
      for (int i = 0; i < data.at(0).size(); i++)
      {
        for (int j = 0; j < data.size(); j ++)
        {
           if (data.at(j).at(i) > border)
           {
             vec.at(i) += 1;
           }
        }
      }
  transform(vec.begin(), vec.end(), vec.begin(), bind1st(multiplies<double>(), 1.0/data.size()));
  transform(vec1.begin(), vec1.end(), vec1.begin(), bind1st(multiplies<double>(), 1.0/data.size()));
  transform(vec3.begin(), vec3.end(), vec3.begin(), bind1st(multiplies<double>(), 1.0/data.size()));
  }
  
  
//   for(int i=0; i<vec.size() ; i++)				/Normierung auf Teilchenzahl wird in calcFlux vorgenommen
//   {
//       vec.at(i) = vec.at(i)/leftParticles.at(i); 
//   }
//   for(int i=0; i<vec.size() ; i++)
//   {
//       if(i==0)
//       {
// 	vec.at(i)=0.0;
//       }
//       else
//       {
// 	vec.at(i) = vec.at(i)/leftParticles.at(i-1);
//       }
//   }
  /*transform(vec.begin(), vec.end(), vec.begin(), bind1st(multiplies<double>(), 1.0/data.size()));
  transform(vec1.begin(), vec1.end(), vec1.begin(), bind1st(multiplies<double>(), 1.0/data.size()));
  transform(vec3.begin(), vec3.end(), vec3.begin(), bind1st(multiplies<double>(), 1.0/data.size()));*/ //möglicherweise durch leftParticles.at(i) teilen!

//transition value for the kth iteration
int n = vec.size()-1;
transprob.at(k) = vec.at(n);
}

void calcFlux(vector<double>& vec,const vector<double>& vec1,const vector<double>& rightDist,const double& dt,const vector< vector<double> >& data, const double& border,const bool& Kramers) 
{/*
  calcFlux(flux,rightDistKramers,rightDist,so.dt);*/
  //vec1: reinitialized particles
  vector<double> leftParticles;
  leftParticles.resize(rightDist.size(),0.0);
  vec.resize(rightDist.size(),0.0);
  countLeft(leftParticles,data,border);  //Zähle Teilchen auf linker Seite für jeden Zeitschritt i -> für Normierung 

  for(int i=0; i< rightDist.size(); i++)
  { 
    
    if(i==0)
    {
    vec.at(i)=0.0;
    }
    else
    {
      vec.at(i)=rightDist.at(i)-rightDist.at(i-1);
    }
  }
    
  if(Kramers)
  {
    for(int i=0; i<vec.size() ; i++)
    {
	if(i==0)
	{
	  vec.at(i)=0.0;
	}
	else
	{
	  vec.at(i) = vec.at(i)/(leftParticles.at(i-1)*dt);
	}
    }
  }
  else
  {
      transform(vec.begin(), vec.end(), vec.begin(), bind1st(multiplies<double>(), 1.0/dt));
  }
}

//counts number of the particles crossing the potential well for every time step
// mit der Option "Kramers" verliert der linksseitige Fluss seine 
void fluxLeftRight(vector <double>& vec,vector <double>& vec1,const double border,const vector< vector<double> >& data,const double& dt)
{ // sollte doch richtig sein !!!!!!
 vec.resize(data.at(0).size(),0.0); 
 vec1.resize(data.at(0).size(),0.0);
 vector<double> vecleft;
 vecleft.resize(data.at(0).size(),0.0); 
 vector<int> vec_side;
 vec_side.resize(data.size(),2.0);
 int k=0;
 int l=0;
 
  for (int i = 0; i < data.at(0).size(); i++) // zähle Teilchen links von der Schwelle; i Zeitpunkte
      {
        for (int j = 0; j < data.size(); j ++) //j Teilchen
        {
           if (data.at(j).at(i) <= border)     //für alle Teilchen wird geschaut, auf welcher Seite sie sich zum Zeitpunkt i befinden
           {
             vecleft.at(i) += 1;   // Anzahl der Teilchen zum Zeitpunkt i
           }
        }
      }
      
 for (int i=0; i<data.at(0).size(); i++){   //erster Durchlauf weist jedem Teilchen seine augenblickliche Position zu ; 0:links, 1:rechts
   for(int j=0; j<data.size(); j++){        //bei zweitem und nachfolgedem Durchlauf tragen dei Teilchen zum rechtsseitigen Fluss be
     if(vec_side.at(j)==0 && data.at(j).at(i)>border){  // die bei vorherigem Zeitschritt auf der linken Seite waren ( vec.at(i)=0), deren Position
       k +=1;}						// zum aktuellen Zeitpukt aber auf der rechten Seite ist
       //cout << k << endl;
     if(vec_side.at(j)==1 && data.at(j).at(i)<border){
       l +=1;}
       //cout << l << endl;
     if(data.at(j).at(i)<border){       // if-else-Schleife weist j-ten Teilchen seine Position zum Zeitpunkt i zu
     vec_side.at(j)=0;}			// Teilchen links der Grenze
     else{					
     vec_side.at(j)=1;}			//Teilchen rechts der Grenze
   }
   //cout << i << " " << k << " " << l << endl;
   vec.at(i)=k;  // Anzahl der Teilchen, die von i-1 nach i auf rechte Seite gewandert sind
   vec1.at(i)=l; // Anzahl der Teilchen, die von i-1 nach i auf linke Seite gewandert sind
   k=0;
   l=0;
 }

  for(int i=0; i<vec.size(); i++)  // Normierung der Flüsse zum Zeitpunkt i auf Teilchenzahl der linken Seite zum vorhergehenden Zeitpunkt i-1
  {
   if(i==0)
   {
    vec.at(i)=vec.at(i);
   }
   else 
   {
    vec.at(i)=vec.at(i)/vecleft.at(i-1); 
   }
  }
  
 transform(vec.begin(), vec.end(), vec.begin(), bind1st(multiplies<double>(), 1.0/dt)); 
//  transform(vec.begin(), vec.end(), vec.begin(), bind1st(multiplies<double>(), 1.0/(data.size()*dt)));
 for(int a=0;a<vec.size();a++) // negative rechtsseitige Flüsse müssen vermieden werden, da diese lediglich numerischen Ursprungs sind!
	{
	  if(vec.at(a)<0){
	  vec.at(a)=0;}
	}  
//  transform(vec1.begin(), vec1.end(), vec1.begin(), bind1st(multiplies<double>(), 1.0/data.size()));
//  for(int a=0;a<vec1.size();a++)
// 	{
// 	  if(vec1.at(a)<0){
// 	  vec1.at(a)=0;}
// 	} 
 /*fstream f;
 f.open("flux.dat",ios::out);
 for(int n=0;n<vec.size();n++){
 f << n << " " << vec.at(n) << " " << vec1.at(n) << endl;}*/
}

//statistic of FluxPaper oder eigene Methode 
void averageFluxLeftRight(vector<double>& vec, double& averageRate,double& variance, const double& nSteps)
{
  int i;
  int k=0; 
  int l=0;
  averageRate = 0.0;
  variance=0.0;
  l=nSteps*0.6;
  for(i=l; i < vec.size(); i++) // bei t=25 hat der Fluss seine Grenzwert erreicht (nSteps=1000)
  {
    averageRate += vec.at(i); 
    k += 1;
  }
  cout << k << endl;
  averageRate = averageRate/k; // Mittelung über alle verbleibenden 800 Zeitschritte
  
  for(i=l; i < vec.size();i++)
  {
    variance +=  pow(vec.at(i)-averageRate,2.0);
  }
  variance=1.0/sqrt(k)*sqrt(1.0/(k-1)*variance);
}

//KramersRate wie in Paper berechnet
void KramersFluxPaper(vector<double>& rate, const vector < vector<double> >& allX, const double& dt, const int& np, const double& Kborder)
{
  
//   fstream Paper1,Paper2;
//   Paper1.open("ParticlesPaper.dat",ios::out);
//   Paper2.open("fluxPaper.dat",ios::out);
  
  int i;  //Zeiten
  int j;	//Teilchen
  int absParticles=0;  				//Zahl absorbierter Teilchen pro Zeitschritt
  int absParticlesTotal=0; 			//Gesamtzahl absorbierter Teilchen
  int nSteps=0;
  
  vector<int> abs(allX.at(0).size(), 0);			//Zahl absorbierter Teilchen für jedes Zeitinkrement
  vector<int> absTotal(allX.at(0).size(), 0);
  vector<int> alrAbs(allX.size(), 0);				//Teilchen schon absorbiert? 0: nein; 1: ja
  rate.resize(allX.at(0).size(), 0.0);

  //für jeden Zeitschritt werden nacheinander alle Teilchen durchgegangen 
  for(i=0;i<allX.at(0).size();i++) //Zeitschritte
  {
    for(j=0;j<allX.size();j++)		//Teilchen j
    {
      if(allX.at(j).at(i)>=Kborder && alrAbs.at(j)==0)
      {
	absParticles++;
	absParticlesTotal++;
	alrAbs.at(j)=1;		//Teilchen, die schon absorbiert wurden, werden nicht nochmal gezählt
	
      }
    }
//     Paper1 << "absPraticles " << absParticles << " absParticlesTotal " << absParticlesTotal << endl;
    abs.at(i)=absParticles;
    absTotal.at(i)=absParticlesTotal;
    absParticles=0;
  }

  for(i=0; i<rate.size(); i++)
  {
    rate.at(i)=1.0/(np-absTotal.at(i))*abs.at(i)/dt; //per Definition größer oder gleich Null!
    /*Paper2 << "rate.at(i) " << rate.at(i) << endl;*/ 
//     if(np==absTotal.at(i))
//     {
//     cout << "divergent, da n0=na" << endl;
//     }
  }

}



// returns the averaged potential energy over all times
void calcPotEnergyAv(vector<double>& vec, const vector< vector<double> >& allX, const vector<double>& tVec, function<double(double, double)> potentialFunc)
{
  vector< vector<double> > potEnergy(allX.size(), vector<double>(allX.at(0).size(), 0.0));
  for (int i =0; i < allX.size(); i++)
  {
    for (int j = 0; j < allX.at(0).size(); j++)
    {
      potEnergy.at(i).at(j) = potentialFunc(tVec.at(j), allX.at(i).at(j));
    }
  }
  buildAverage(vec, potEnergy);
}

//calculates kinetic Energy out of a velocity Vector
void calcKinEnergy(vector<double>& vec, const vector<double>& velocVec, const double& mass)
{
  vec.resize(velocVec.size(), 0.0);
  vectorSquared(vec, velocVec);
  transform(vec.begin(), vec.end(), vec.begin(), bind1st(multiplies<double>(), 0.5*mass));
}

void calcKinEnergyAv(vector<double>& vec, const vector< vector<double> >& allVeloc, const double& mass)
{
  vector< vector<double> > tempSquared(allVeloc.size(), vector<double>(allVeloc.at(0).size(), 0.0));
  vec.resize(allVeloc.at(0).size(), 0.0);
  for (int i = 0; i < allVeloc.size(); i++){
	  calcKinEnergy(tempSquared.at(i), allVeloc.at(i), mass);
  }
  buildAverage(vec, tempSquared);
//   transform(output.begin(), output.end(), output.begin(), bind1st(multiplies<double>(), 0.5*mass));
}

// calculates potential energy of all particles at one specific time t
// written for prepareEnergyAnimation function
void calcPotEnergy(vector<double>& vec, const vector<double>& xVec, const double& t, function<double(double, double)> potentialFunc)
{
  vec.resize(xVec.size(), 0.0);
  for (int i = 0; i < xVec.size(); i++)
  {
    vec.at(i) = potentialFunc(t, xVec.at(i));
  }
}

//calculates total Energy of all particles at one specific time t
//written for prepareEnergyAnimation function
void calcTotalEnergy(vector<double>& vec, const vector<double>& xVec, const vector<double>& velocVec, const double& t,
			       const double& mass, function<double(double, double)> potentialFunc)
{
  vector<double> kinEnergyVec(xVec.size(), 0.0);
  vector<double> potEnergyVec(xVec.size(), 0.0);
  vec.resize(xVec.size(), 0.0);
  calcKinEnergy(kinEnergyVec, velocVec, mass);
  calcPotEnergy(potEnergyVec, xVec, t, potentialFunc);
  add2Vec(vec, kinEnergyVec, potEnergyVec);
}

//I(0,t') = Xi(0)*Xi(t')
void calcCorr(vector<double>& vec, const vector<double>& data, const int& startN, const int& endN)
{
  vec.resize(data.size(), 0.0);
  transform(data.begin()+startN, data.end()-endN, vec.begin(), bind1st(multiplies<double>(), data.at(startN)));
}

//calculates correlation function I(t)=xi(0)*xi(t) and averages it
void buildCorrAv(vector<double>& vec, const vector< vector<double> >& allNoise, const int& startN, const int& endN)
{
  vector< vector<double> > tempVec(allNoise.size(), vector<double>(allNoise.at(0).size(), 0.0));
  for (int i=0; i < allNoise.size();i++){
	  calcCorr(tempVec.at(i), allNoise.at(i), startN, endN);
  }
  buildAverage(vec, tempVec);
}

//----animation functions----

//writes for each time one file for gnuplot with position as x-axis and total energy at y-axis
int prepareEnergyAnimation(const vector< vector<double> >& allX, const vector< vector<double> >& allVeloc, const vector<double>& tVec,
			       const double& mass, function<double(double, double)> potentialFunc, const int& stride, const string& headerLine, const string& folderName)
{
  vector<double> totalEnergy(allX.size(), 0.0);
  vector<double> yVec(allX.size(), 0.0); //position of all particles at one specific time
  vector<double> velocVec(allX.size(), 0.0); //velocity of all particles at one specific time

  for (int i = 0; i < allX.at(0).size(); i += stride)
  {
    getRow(yVec, allX, i);
    getRow(velocVec, allVeloc, i);
//     totalEnergy = calcPotEnergy(yVec, tVec.at(i), potentialFunc);
    calcTotalEnergy(totalEnergy, yVec, velocVec, tVec.at(i),  mass, potentialFunc);
    writeToFile(yVec, totalEnergy, "energyAnimation", headerLine, folderName);
  }
  return 0;
}

//same as prepareEnergyAnimation && writes potential function in file for each time step
int prepareEnergyAnimationTimePot(const vector< vector<double> >& allX, const vector< vector<double> >& allVeloc, const vector<double>& tVec,
			       const double& mass, function<double(double, double)> potentialFunc, const int& potentialN, const double& potBorderL,
				  const double& potBorderR, const int& skipN, const string& headerLine, const string& folderName)
{
  prepareEnergyAnimation(allX, allVeloc, tVec, mass, potentialFunc, skipN, headerLine, folderName);
  vector<double> potentialRange(potentialN, 0.0);
  linspace(potBorderL, (potBorderR-potBorderL)/potentialN, potentialRange);
  vector<double> potentialVec(potentialN, 0.0);
  double t;
  for (int i = 0; i < tVec.size(); i += 1 + skipN)
  {
    for (int j = 0; j < potentialN; j++)
    {
      t = tVec.at(i);
      potentialVec.at(j) = potentialFunc(t, potentialRange.at(j));
    }
    writeToFile(potentialRange, potentialVec, "potential", headerLine, folderName);
  }
  return 1;
}

//writes one file for gnuplot of the velocity distribution for each time step
int prepareVda(const vector< vector<double> >& allVeloc, const vector<double>& vDistT,vector <vector <double> >& DistTotal,
	  const double& rangeBeg, const double& rangeEnd, const int& nBins, const int& stride,
	  const string& headerLine, const string& folderName)
{
  vector<double> vDist(nBins, 0.0);
  vector<double> vDistNorm(nBins, 0.0);
  vector<double> velocVec(allVeloc.size(), 0.0); //velocity of all particles at one specific time
   
  int j =0; /* Schleife läuft in "stride"-Schritten, um zu gewährleisten, dass Einträge von DistTotal.at(j) 
  (entspricht ksim.vDistAvVec.at(i)) richtig gefüllt wird, definiere zweiten Schleifenindex der von 0 bis allVeloc.at(0).size()/stride
  geht */
  for (int i = 0; i < allVeloc.at(0).size(); i += stride)
  {
    getRow(velocVec, allVeloc, i);
    calcVDistribution(vDist,vDistNorm, velocVec, rangeBeg, rangeEnd, nBins);
    DistTotal.at(j)=vDist;
    //cout << DistTotal.at(0).size() << " " << DistTotal.size() << endl;
     writeToFile(vDistT, vDist, "vDist", headerLine, folderName); 
//     writeToFile(vDistT, vDistNorm, "vDistNorm", headerLine, folderName);  
    j++;
     //cout << j << endl;
  }
  return 0;
}

//writes one file for gnuplot of the position distribution for each time step
int preparePda(const vector< vector<double> >& allx, const vector<double>& pDistT,vector <vector <double> >& DistTotal,
	  const double& rangeBeg, const double& rangeEnd, const int& nBins, const int& stride,
	  const string& headerLine, const string& folderName)
{
  vector<double> pDist(nBins, 0.0);
  vector<double> pDistNorm(nBins, 0.0);
  vector<double> xVec(allx.size(), 0.0); //velocity of all particles at one specific time
  int j=0;
  for (int i = 0; i < allx.at(0).size(); i += stride)
  {
    getRow(xVec, allx, i);
    calcVDistribution(pDist,pDistNorm, xVec, rangeBeg, rangeEnd, nBins);
    DistTotal.at(j)=pDist;
     writeToFile(pDistT, pDist, "pDist", headerLine, folderName);
//     writeToFile(pDistT, pDistNorm, "pDistNorm", headerLine, folderName);
    j++;
  }
  return 0;
}

//discretized bessel Function K0, used for comparison with numerical results
void besselFunction(vector<double>& vec, const vector<double>& t_vec)
{
  vec.resize(t_vec.size(), 0.0);
  transform(t_vec.begin(), t_vec.end(), vec.begin(), gsl_sf_bessel_K0);
  transform(vec.begin(), vec.end(), vec.begin(), bind1st(multiplies<double>(),1.0/M_PI));
}

// discretized fourier transform of the exponential correlation function
void fourier1(vector<double>& vec, const double& samplingRate, const int& sampleN, const int& nSteps)
{
  vec.resize(nSteps, 0.0);
  double dw = M_PI/((double)samplingRate)/((double)sampleN);
  double w = 0;
  for (int i = 0; i < nSteps; i++)
  {
      vec.at(i) = 1.0/(1+pow(w,2));
      w += dw;
  }
}

// analytically expected kinetic energy for exponential correlation function as described
void calcKinEnergyTheo(vector<double>& vec, const vector<double>& tVec, const double& temperature,
		     const double& mass, const double& D, const double& tau, const double& k)
{
  vec.resize(tVec.size(), 0.0);
  double Q = D/mass/k/temperature;
  if (2*Q*tau < 1)
  {
    printf("cannot calculate analytically as 2*Q*tau < 1\n");
  }
  else
  {
    double omega = sqrt(2*Q*tau-1)/(2*tau);
    for (int i = 0; i < tVec.size(); i++)
    {
      vec.at(i) = k*temperature/mass-k*temperature/(mass*(2*Q*tau-1))*(Q*tau+sqrt(2*Q*tau-1)*sin(2*omega*tVec.at(i))
		      + (Q*tau-1)*cos(2*omega*tVec.at(i)))*exp(-tVec.at(i)/tau);
    }
  }
  transform(vec.begin(), vec.end(), vec.begin(), bind1st(multiplies<double>(), 0.5*mass));
}

//calculates theoretical spectral density of velocity for colored noise and no potential
void calcSpecVTheo0(vector<double>& vec, const vector<double>& omegaVec, const double& k_b, const double& temperature, const double& mass, const double& D, const double& tau)
{
  vec.resize(omegaVec.size(), 0.0);
  complex<double> g_ret;
  complex<double> gamma;
  for (int i = 0; i < omegaVec.size(); i++)
  {
    gamma = D/(4.0*k_b*temperature)/complex<double>(1.0, -omegaVec.at(i)*tau);
    g_ret = 1.0/(2.0*gamma-complex<double>(0.0, mass*omegaVec.at(i)));
    vec.at(i) = 2.0*k_b*temperature * real(g_ret);
  }
}

//calculates theoretical spectral density of velocity for colored noise and harmonic potential k*x^2
void calcSpecXTheo1(vector<double>& vec, const vector<double>& omegaVec, const double& k_b, const double& temperature, const double& mass, const double& D, const double& tau, const double& potK)
{
  vec.resize(omegaVec.size(), 0.0);
  complex<double> g_ret;
  complex<double> gamma;
  for (int i = 0; i < omegaVec.size(); i++)
  {
    gamma = D/(4.0*k_b*temperature)/complex<double>(1.0, -omegaVec.at(i)*tau);
    g_ret = 1.0/(complex<double>(potK-mass*omegaVec.at(i)*omegaVec.at(i), 0.0) - complex<double>(0.0,1.0)*2.0*omegaVec.at(i)*gamma);
    vec.at(i) = 2.0*k_b*temperature / omegaVec.at(i) * imag(g_ret);
  }
}

// returns spectral density of given vector corrVec which is the correlation of some function
void calcSpec(vector<double>& vec, const vector<double>& corrVec, const double& dt)
{
  int steps = corrVec.size();
  vector<double> doubleCorrVec(2*steps, 0.0);
  mirrorData(doubleCorrVec, corrVec, false, false);
  double* specArray =(double*) calloc(2*steps, sizeof(double));
  copy(doubleCorrVec.begin(), doubleCorrVec.end(), specArray);
  gsl_fft_real_radix2_transform (specArray, 1, 2*steps);
  for (int i = 0; i < steps+1; i++)
  {
     specArray[i] = specArray[i]*dt;
  }
  vec.assign(specArray,specArray + steps);
}

//calculates <Delta x^2> for white Noise without potential and v0 = 0
void  calcXSqrTheo0(vector<double>& vec, const vector<double>& tVec, const double& mass, const double& gamma, const double& k_b, const double& temperature)
{
  vec.resize(tVec.size(), 0.0);
  for (int i = 0; i < tVec.size(); i++)
  {
    vec.at(i) = 2*mass*k_b*temperature/(gamma*gamma) * (2*exp(-gamma/mass*tVec.at(i)) - 0.5*exp(-2*gamma/mass*tVec.at(i)) + gamma/mass*tVec.at(i) -1.5);
  }
}

//calculates <1/2*mass* v^2> for white Noise without potential
void calcKinEnergyTheo0(vector<double>& vec, const vector<double>& tVec, const double& mass, const double& gamma, const double& k_b, const double& temperature, const double& v0)
{
  vec.resize(tVec.size(), 0.0);
  for (int i = 0; i < tVec.size(); i++)
  {
    vec.at(i) = 0.5*k_b*temperature + (0.5*mass*v0*v0-0.5*k_b*temperature)*exp(-2*gamma/mass*tVec.at(i));
  }
}

//Kramers



//does all calculations after simulating particles, like statistics and so on
void doAfterMath(const Filenames& filenames,const Foldernames& foldernames, const SimulationOptions& so, const Potential& potential, const Results& results, const NoiseDiss& noiseDiss, const string& headerString, KappaSimulation& ksim){
    //percentage of particles being at y > 0
    vector<double> rightDist(so.nSteps);
    vector<double> reinitPart(so.nSteps);
    vector<double> rightDistTotal(so.nSteps);
    vector<double> rightDistKramers(so.nSteps);
    int j= ksim.indexAvNum;
    int i = ksim.indexNr;
    countRight(rightDist,reinitPart,ksim.transition,rightDistKramers, results.allX, so.xb,so.rxborder,so.KramersRate,i,so.dt);
    ksim.transitionAvVec.at(j) = rightDist;
    
    vector<double> LeftParticles(so.nSteps);
    countLeft(LeftParticles,results.allX,so.xb);
    
 /*   
    for(int i=0; i<LeftParticles.size();i++)  //Schleife gibt na für jeden Zeitpunkt aus!
    {
     cout << LeftParticles.at(i) << endl; 
    }*/
//     int n= rightDistTotal.size()-1;
//         double dE=abs(0.50355-rightDist.at(n));
//         if(dE<=pow(10,-2))
// 	{
// 	  cout << "funktioniert Schleife" << endl;
// 	  fstream f;
//     f.open("correlTime.dat", ios::out | ios::app);
// 	 f << "alpha:" << " " << so.alpha << " " << "rightDist:"<< " " << rightDist.at(n)<< " " << "dE" << " " << dE << endl; 
// 	}
//     if(j==so.avNum-1)
// 	{
// 	buildAverage(rightDistTotal,ksim.transitionAvVec);
// 	writeToFile(results.tVec,rightDistTotal,filenames.rightDistTotal,headerString,foldernames.main);
//         double dE=abs(0.50355-rightDistTotal.at(n));
//         if(dE<=pow(10,-3))
// 	{
// 	  fstream F;
//     F.open("correlTimeTotal.dat", ios::out | ios::app);
// 	 F << "alpha:" << " " << so.alpha << " " << "rightDistTotal:"<< " " << rightDistTotal.at(n)<< " " << "dE" << " " << dE << endl; 
// 	}
// 	}
   // writeToFile(results.tVec, rightDist, filenames.rightDist, headerString, foldernames.main);
    //writeToFile(results.tVec, rightDistKramers, filenames.rightDistKramers, headerString, foldernames.main);
    //writeToFile(results.tVec, reinitPart, filenames.reinitPart, headerString, foldernames.main);
	
    
    
    //===================================Corr3============================================

	//Grenzwert der mittleren kinetische Energie  für masselose Theorie (analytisch)
	//
// 	fstream Energy;
//          Energy.open("kinEnergy.dat", ios::out);
         double v2;
         v2=4.0/9.0*so.temperature/so.mass+4.0/9.0*pow(so.v0,2.0);
//          for (int i=0; i<results.tVec.size();i++)
// 	 {
// 	   Energy << results.tVec.at(i) << " " << 1.0/2.0*so.mass*v2 << endl;
// 	 }
	 
   //-----------------------------Corr3----------------------------------
	 
	 

	//======================================noise=======================================================
	//writeToFile(results.tVec, noiseDiss.allNoise.at(0), filenames.noise, headerString, foldernames.main);

	//noise average
	vector<double> noiseAv(so.nSteps, 0.0);
	vector<double> noiseAvSqr(so.nSteps, 0.0); //squared average
	vector<double> noiseSqrAv(so.nSteps, 0.0); //average of squared vector
	vector<double> noiseAvSqrTot(so.nSteps, 0.0);
	vector<double> noiseSqrAvTot(so.nSteps, 0.0);
	vector<double> noiseAvTotal(so.nSteps, 0.0);
	buildAverage(noiseAv,noiseDiss.allNoise);
	buildAverageSquared(noiseSqrAv,noiseDiss.allNoise);
	ksim.noiseAvVec.at(j)=noiseAv;
	ksim.noiseSqrAvVec.at(j)=noiseSqrAv;
	if(j==so.avNum-1)
	{
	buildAverage(noiseAvSqrTot,ksim.noiseAvVec);
	vectorSquared(noiseAvSqr,noiseAvSqrTot);
	buildAverage(noiseSqrAvTot,ksim.noiseSqrAvVec);
 	writeToFile(results.tVec,noiseAvTotal,filenames.noiseAvTotal,headerString,foldernames.main); //noiseAvTotal ist leer
	}
// 	for(int i=0; i<noiseAvSqrTot.size();i++)
// 	{
// 	 noiseAvTotal.at(i)= noiseSqrAvTot.at(i)-noiseAvSqr.at(i);	
// 	  cout << noiseAvTotal.at(i) << " " <<  noiseSqrAvTot.at(i) << " " << noiseAvSqrTot.at(i) << " " <<  noiseAvSqr.at(i) << endl;
// 	}
// 	writeToFile(results.tVec, noiseAv, filenames.noiseAv, headerString, foldernames.main);

	//========================================correlation of colored noise================================================
	vector<double> corrAv(so.nSteps, 0.0);
	vector<double> corrAvTotal(so.nSteps, 0.0);
	//void buildCorrAv(vector<double>& vec, const vector< vector<double> >& allNoise, const int& startN, const int& endN)
	buildCorrAv(corrAv, noiseDiss.allNoise, 0, 0);
	ksim.CorrAvVec.at(j)=corrAv;
	if(j==so.avNum-1)
	{
	buildAverage(corrAvTotal,ksim.CorrAvVec);
	writeToFile(results.tVec,corrAvTotal,filenames.corrAvTotal,headerString,foldernames.main);
	}
// 	writeToFile(results.tVec, corrAv, filenames.corrAv, headerString, foldernames.main);

	//single simulation - position and velocity
        //writeToFile(results.tVec, results.allX.at(0), filenames.x, headerString, foldernames.main);
        //writeToFile(results.tVec, results.allX.at(1), filenames.x, headerString, foldernames.main); // generator check
// 	writeToFile(results.tVec, results.allV.at(0), filenames.veloc, headerString, foldernames.main);
// 	writeToFile(results.tVec, results.allV.at(20), filenames.veloc, headerString, foldernames.main);
	

	//=======================================flux=========================================================
	vector<double> flux(so.nSteps,0.0);
	vector<double> fluxPositive(so.nSteps,0.0);
	vector<double> fluxNegative(so.nSteps,0.0);
	vector<double> fluxPaper(so.nSteps,0.0);
	vector<double> fluxTotal(so.nSteps,0.0);
 	vector<double> fluxPositiveTotal(so.nSteps,0.0);
 	vector<double> fluxNegativeTotal(so.nSteps,0.0);
	vector<double> fluxPaperTotal(so.nSteps,0.0);
	double averageKramers=0.0;
	double variance=0.0;
// 	fstream E;
//          E.open("fluxtotvorher.dat", ios::out | ios::app);
//  	for(int i=0; i < fluxTotal.size(); i++)
//  	{E << fluxTotal.at(i) << endl;}
	calcFlux(flux,rightDistKramers,rightDist,so.dt,results.allX, so.xb,so.KramersRate);
        fluxLeftRight(fluxPositive,fluxNegative,so.xb,results.allX,so.dt);
	KramersFluxPaper(fluxPaper, results.allX, so.dt, so.np, so.rxborder);
	ksim.fluxAvVec.at(j)=flux;
	ksim.fluxPositiveAvVec.at(j)=fluxPositive;
	ksim.fluxNegativeAvVec.at(j)=fluxNegative;
	ksim.fluxPaperAvVec.at(j)=fluxPaper;
// 	if(j==0){
// 	fstream A;
//         A.open("eins.dat", ios::out | ios::app);
// 	for (int i=0; i< ksim.fluxPositiveAvVec.at(0).size(); i++){
// 	A << results.tVec.at(i) <<  " " << i << " " << ksim.fluxPositiveAvVec.at(j).at(i) << endl; }}
// 	
// 	if(j==1){
// 	fstream B;
//         B.open("zwei.dat", ios::out | ios::app);
// 	for (int i=0; i< ksim.fluxPositiveAvVec.at(0).size(); i++){
// 	B << results.tVec.at(i) <<  " "  << i << " " << ksim.fluxPositiveAvVec.at(j).at(i) << endl; }}
// 	
// 	if(j==2){
// 	fstream C;
//         C.open("drei.dat", ios::out | ios::app);
// 	for (int i=0; i< ksim.fluxPositiveAvVec.at(0).size(); i++){
// 	C << results.tVec.at(i) <<  " " << i << " " << ksim.fluxPositiveAvVec.at(j).at(i) << endl; }}
	
	if(j==so.avNum-1)
	{
	    buildAverage(fluxTotal,ksim.fluxAvVec);
	    writeToFile(results.tVec,fluxTotal,filenames.fluxTotal,headerString,foldernames.main);
	    buildAverage(fluxPositiveTotal,ksim.fluxPositiveAvVec); // wird richtig gemittelt, wenn für jeden Fluss eigener Vektor definiert wird
	    for(int a=0;a<fluxPositiveTotal.size();a++)
	    {
	      if(fluxPositiveTotal.at(a)<0){
	      fluxPositiveTotal.at(a)=0;}
	    }  
	    
// 	    if(j==so.avNum-1)
// 	    {
// 	    averageFluxLeftRight(fluxPositiveTotal, averageKramers); //berechnet den Mittelwert des Flusses (beachte, dass der Zeitpunkt, ab dem gemittelt wird, gegebenenfalls angepasst werden muss!!)
// 	    }
	    
	    writeToFile(results.tVec,fluxPositiveTotal,filenames.fluxPositiveTotal,headerString,foldernames.main);
	    buildAverage(fluxNegativeTotal,ksim.fluxNegativeAvVec);
	    for(int a=0;a<fluxNegativeTotal.size();a++)
	    {
	      if(fluxNegativeTotal.at(a)<0){
	      fluxNegativeTotal.at(a)=0;}
	    }  
	    writeToFile(results.tVec,fluxNegativeTotal,filenames.fluxNegativeTotal,headerString,foldernames.main);
	    
	    buildAverage(fluxPaperTotal,ksim.fluxPaperAvVec);
	    ksim.fluxPaperOut=fluxPaperTotal;
	    for(int a=0;a<fluxPaperTotal.size();a++)
	    {
	      if(fluxPaperTotal.at(a)<0){
	      fluxPaperTotal.at(a)=0;}
	    } 
	    
	    if(so.Kramers==0)
	    {
	      averageFluxLeftRight(fluxPaperTotal, averageKramers,variance,so.nSteps);  
	    }
	    else
	    {
	      averageFluxLeftRight(fluxTotal, averageKramers,variance,so.nSteps); 
	    }
	    
	    writeToFile(results.tVec,fluxPaperTotal,filenames.fluxPaperTotal,headerString,foldernames.main);
	}
	//writeToFile(results.tVec,flux,filenames.flux,headerString,foldernames.main);
	//writeToFile(results.tVec,fluxPositive,filenames.fluxPositive,headerString,foldernames.main);
	//writeToFile(results.tVec,fluxNegative,filenames.fluxNegative,headerString,foldernames.main);
	
#ifdef KRAMERS
	//=============================================kramers==========================================================
	
	  if(j==so.avNum-1)  //erst bei letztem Durchlauf ausgeben!
	  {
	    double kramersTheo;
	    double kramersTheo1;
	    double ratio;
	    double kramersLowFric;
	    double inverseKramersTurnover;
	    double kramersTurnover;
	    double a;
	    double gammaFre;
	    double omegaB=sqrt(4*so.Ub/(so.mass*pow((so.xb-so.xc),2.0)));
	    double omegaC=omegaB;
	    gammaFre = so.gamma/so.mass;
	    
	    cout << "beta=gamma/m : " << gammaFre << endl;
	    cout << "xc " << so.xc << endl;
	    cout << "temperature " << so.temperature << endl;
	    cout << "mass " << so.mass << endl;
	    cout << "mass*omegaB² " << so.mass*pow(omegaB,2.0) << endl;
	    if(omegaB==gammaFre)
	      
	    {cout << "omegaB==gamma/m " << gammaFre << endl;}
	    
	    if(so.noiseNr==1)
	    {
	      cout << "WhiteNoise" << endl; 
	      a=sqrt(4.0)/(0.6*2*M_PI)*(sqrt(5.0/4.0)-1.0/2.0);
	      cout << "a " << a << endl;
	      cout << "exp(-Ub/T) " << exp(-so.Ub/(so.k_b*so.temperature)) << endl;
	      cout << "sqrt(Ub/m) " << sqrt(so.Ub/so.mass) << endl;
	      cout << "1/xc " << 1.0/so.xc << endl;
	      cout << "gamma " << so.gamma << " m*wb " << so.mass*omegaB << endl;
	      
	      
	      kramersTheo = omegaC*so.aWhite/(2.0*M_PI*omegaB)*exp(-so.Ub/(so.k_b*so.temperature));
	      kramersLowFric=so.gamma/so.mass*so.Ub/so.temperature*exp(-so.Ub/(so.k_b*so.temperature));
	      inverseKramersTurnover=1.0/kramersTheo+1.0/kramersLowFric;
	      kramersTurnover = 1.0/inverseKramersTurnover;
	      cout << "Ub " << so.Ub << endl;
	      cout << "dt " << so.dt << endl;
	      cout << "beta*Ib " << so.gamma/so.mass*2*M_PI*so.Ub/so.potw << endl;
	      cout << "kramersTheo " << kramersTheo  << endl;
	      cout << "kramersNumAv " << averageKramers << " Varianz " << variance <<  endl;
	      cout << "kramersTheoWeakFric " << kramersLowFric << endl;
	      cout << "kramersTurnover " << kramersTurnover << endl; 
	      if(so.gamma/so.mass*2*M_PI*so.Ub/so.potw<so.temperature)
	      {
		ratio=averageKramers/kramersLowFric;
		cout << "ratio " <<  averageKramers/kramersLowFric << endl;
		cout << "low friction" << endl;
	      }
	      else{
		ratio=averageKramers/kramersTheo;
		cout << "ratio " <<  averageKramers/kramersTheo << endl;
	      }
	    }
	    else if(so.noiseNr==2)
	    {
	      cout << "ColNoise" << endl;
	      cout << "Realteile der Eigenwerte " << endl
		   << "l1 " << so.aColoured << " " << "l2 " << so.bColoured << " " << "l3 " << so.cColoured << endl;
	      kramersTheo =omegaC*so.aColoured/(2.0*M_PI*omegaB)*exp(-so.Ub/(so.k_b*so.temperature));
	      cout << "kramersTheo " << kramersTheo << endl;
	      cout << "kramersNumAv " << averageKramers << " Varianz " << variance <<  endl;
	      ratio=averageKramers/kramersTheo;
	      cout << "ratio " << ratio << endl;
	    }
	    else
	    {
		cout << "Weder farbiges noch weißes Rauschen! " << endl;   
	    }
	  if(so.noiseNr==1)
	  {
	   fstream outputKramersMarkov;
	   outputKramersMarkov.open("outputKramersMarkov.dat", ios::out | ios::app);
	   //outputKramersMarkov << so.Ub  << " " << so.dt <<  " " <<so.gamma/so.mass << " " <<  averageKramers << " " << variance << " " << ratio << " " << kramersTheo << " " << so.rxborder << endl;
	   outputKramersMarkov << so.gamma/so.mass << " " <<  averageKramers << " " << variance << " " << ratio << " " << kramersTheo << " " << so.rxborder << endl;
	   outputKramersMarkov.close();
	  }
	  else
	  {
	    fstream outputKramersColor;
	    outputKramersColor.open("outputKramersColor.dat", ios::out | ios::app);
	    if(so.corrFuncNr==0)
	    {
	      cout << "CorrFunc0" << endl;
	     // outputKramersColor << so.corrFuncNr << " " << so.tau << " " << so.Ub  << " " << so.dt <<  " " << so.gamma/so.mass << " " <<  averageKramers << " " << variance << " " << ratio <<" " << kramersTheo << " " << so.rxborder<<  endl;
	      outputKramersColor << so.corrFuncNr << " " << so.tau << " " << so.Ub  << " " << so.dt <<  " " << so.gamma/so.mass << " " <<  averageKramers << " " << variance << " " << ratio <<" " << kramersTheo << " " << so.rxborder<<  endl;
	      
	    }
	    else if(so.corrFuncNr==1)
	    {
	      cout << "CorrFunc1" << endl;
	      //outputKramersColor << so.corrFuncNr << " " << so.a << " " << so.Ub  << " " << so.dt <<  " " << so.gamma/so.mass << " " <<  averageKramers << " " << variance << " " << ratio <<" " << kramersTheo << " " << so.rxborder << endl;
	      outputKramersColor << so.gamma/so.mass << " " <<  averageKramers << " " << variance << " " << ratio <<" " << kramersTheo << " " << so.rxborder << endl;   
	    }
	    else
	    {
	      cout << "CorrFunc3" << endl;
// 	      outputKramersColor << so.corrFuncNr << " " << 1.0/(so.alpha/sqrt(so.mass)) << " " << so.Ub  << " " << so.dt <<  " " << so.gamma/so.mass << " " <<  averageKramers << " " << variance << " " << ratio <<" " << kramersTheo << " " << so.rxborder <<  endl;
	      outputKramersColor << "D " << so.D <<  averageKramers << " " << variance << " " << ratio <<" " << kramersTheo << " " << so.rxborder <<  endl;
	     	    }
	    outputKramersColor.close();   
	  }
	  }
	
	//-============================================================================
#endif
	
	
	//============================================kinetic energy========================================
        vector<double> kinEnergyAv(so.nSteps, 0.0);
	calcKinEnergyAv(kinEnergyAv, results.allV, so.mass);
	ksim.kinEnergyAvVec.at(j)= kinEnergyAv;
	vector<double> kinEnergyTotal(so.nSteps,0.0);
	
// 	if(j==0){
//  	fstream A;
//          A.open("eins.dat", ios::out | ios::app);
//  	for (int i=0; i< ksim.fluxPositiveAvVec.at(0).size(); i++){
//  	A << results.tVec.at(i) <<  " " << i << " " << ksim.kinEnergyAvVec.at(j).at(i) << endl; }}
//  	
//  	if(j==1){
//  	fstream B;
//          B.open("zwei.dat", ios::out | ios::app);
//  	for (int i=0; i< ksim.fluxPositiveAvVec.at(0).size(); i++){
//  	B << results.tVec.at(i) <<  " "  << i << " " << ksim.kinEnergyAvVec.at(j).at(i) << endl; }}
	
	if(j==so.avNum-1)
	{
	buildAverage(kinEnergyTotal,ksim.kinEnergyAvVec);
	writeToFile(results.tVec, kinEnergyTotal, filenames.kinEnergyAv, headerString, foldernames.main);
	}
	
	double kinEnergyAvAv;
	double kinEnergyStdDeviation;
	kinEnergyAvAv = buildAverage2(kinEnergyAv.begin() + ceil(so.timeSettled/so.dt), kinEnergyAv.end());
	kinEnergyStdDeviation = stdDeviation(kinEnergyAv.begin() + ceil(so.timeSettled/so.dt), kinEnergyAv.end());
	//headerLine << "\n#<kinEnergy> after " << so.timeSettled << "s = " << kinEnergyAvAv << "\n#stdDeviation = " << kinEnergyStdDeviation;
//         writeToFile(results.tVec, kinEnergyAv, filenames.kinEnergyAv, headerString, foldernames.main);

	
	
	//=========================================energy distribution==================================================
	if(so.avNum==1)
	{
	vector<double> eDistT(so.vDistNBins);
	linspace(so.eDistRangeBeg + 0.5*(so.eDistRangeEnd-so.eDistRangeBeg)/so.eDistNBins,
					  (so.eDistRangeEnd-so.eDistRangeBeg)/so.eDistNBins, eDistT);
	vector<double> allKinEnDist(so.vDistNBins);
	vector<double> allKinEnDist1(so.vDistNBins);
	vector< vector<double> >  allKinEn(so.np, vector<double>(so.nSteps, 0.0));
	for(int a=0;a<so.np;a++)
	{
	 calcKinEnergy(allKinEn.at(a),results.allV.at(a),so.mass); 
	}
	
#ifdef KINEN
	fstream EnergyDist,EnergyDist1;
        EnergyDist.open("kinEnergyDist.dat", ios::out);
	EnergyDist1.open("kinEnergyDist1.dat", ios::out);
	for(int a=0; a<so.nSteps;a++)
	{
	 EnergyDist << a << " " << results.allV.at(0).at(a) << " " << 1.0/2.0*so.mass*results.allV.at(0).at(a)*results.allV.at(0).at(a) << " " << allKinEn.at(0).at(a) << endl;
	 EnergyDist1 << a << " " << results.allV.at(3).at(a) << " " << 1.0/2.0*so.mass*results.allV.at(3).at(a)*results.allV.at(3).at(a) << " " << allKinEn.at(3).at(a) << endl;
	}
	EnergyDist.close();
	EnergyDist1.close();
#endif	
	calcvDistributionMean(allKinEnDist,allKinEnDist1, allKinEn, ceil((so.eDistStartTime-so.t0)/so.dt), so.nSteps, so.eDistRangeBeg, so.eDistRangeEnd, so.eDistNBins);
	
	if(so.edaBool==1)
	{
        writeToFile(eDistT, allKinEnDist, filenames.allKinEnDist, headerString, foldernames.main);
	}
	}
	//=========================================velocity distribution================================================
        vector<double> vDistT(so.vDistNBins);
	linspace(so.vDistRangeBeg + 0.5*(so.vDistRangeEnd-so.vDistRangeBeg)/so.vDistNBins,
					  (so.vDistRangeEnd-so.vDistRangeBeg)/so.vDistNBins, vDistT);
	vector<double> vDist(so.vDistNBins);
	vector<double> vDist1(so.vDistNBins);
	vector<double> velocInitDist(so.vDistNBins);
	vector<double> velocInitDist1(so.vDistNBins);
	calcVDistribution(velocInitDist,velocInitDist1, results.v0Vec, so.vDistRangeBeg, so.vDistRangeEnd, so.vDistNBins);
	calcvDistributionMean(vDist,vDist1, results.allV, ceil((so.vDistStartTime-so.t0)/so.dt), so.nSteps, so.vDistRangeBeg, so.vDistRangeEnd, so.vDistNBins);
	//vDist erst dann sinvoll, wenn System bereits equilibriert ist ,ceil((so.vDistStartTime-so.t0)/so.dt) ist Abschätzung für Schritt, ab dem System im Gleichgewicht ist!!
	vector<double> vDistTheo(so.vDistNBins);
	calcVDistTheo(vDistTheo, so.temperature, so.mass, so.k_b, so.vDistRangeBeg, so.vDistRangeEnd, so.vDistNBins);
	if(so.vdaBool==1)
	{
        writeToFile(vDistT, vDist, filenames.vDist, headerString, foldernames.main);
	writeToFile(vDistT, vDist1, filenames.vDist1, headerString, foldernames.main);
	if(ksim.indexAvNum==0)
	{
	writeToFile(vDistT, vDistTheo, filenames.vDistTheo, headerString, foldernames.main);
	writeToFile(vDistT, velocInitDist, filenames.velocInitDist, headerString, foldernames.main);
	writeToFile(vDistT, velocInitDist1, filenames.velocInitDist1, headerString, foldernames.main);
	}
	}
	
	//======================================position distribution=======================================================
        vector<double> pDistT(so.pDistNBins);
	linspace(so.pDistRangeBeg + 0.5*(so.pDistRangeEnd-so.pDistRangeBeg)/so.pDistNBins,
					  (so.pDistRangeEnd-so.pDistRangeBeg)/so.pDistNBins, pDistT);
	vector<double> pDist(so.pDistNBins);
	vector<double> pDist1(so.pDistNBins);
	vector<double> xInitDist(so.vDistNBins);
	vector<double> xInitDist1(so.vDistNBins);
	calcVDistribution(xInitDist,xInitDist1, results.x0Vec, so.pDistRangeBeg, so.pDistRangeEnd, so.pDistNBins);
	calcvDistributionMean(pDist,pDist1, results.allX, ceil((so.pDistStartTime-so.t0)/so.dt), so.nSteps, so.pDistRangeBeg, so.pDistRangeEnd, so.pDistNBins);
	vector<double> xDistTheo(so.vDistNBins);
	calcXDistTheo(xDistTheo, so.temperature, so.mass, so.potw, so.k_b, so.vDistRangeBeg, so.vDistRangeEnd, so.pDistNBins);
	if(so.pdaBool==1)
	{
	writeToFile(pDistT, pDist, filenames.pDist, headerString, foldernames.main);
	writeToFile(pDistT, pDist1, filenames.pDist1, headerString, foldernames.main);
	if(ksim.indexAvNum==0)
	{
	writeToFile(pDistT, xDistTheo, filenames.xDistTheo, headerString, foldernames.main);
	writeToFile(pDistT, xInitDist, filenames.xInitDist, headerString, foldernames.main);
	writeToFile(pDistT, xInitDist1, filenames.xInitDist1, headerString, foldernames.main);
	}
	}
	
	//===================================averaged position============================================================
	vector<double> xAvTheo(so.nSteps,0.0);
	vector<double> xAv(so.nSteps, 0.0);
	vector<double> xAvTotal(so.nSteps,0.0);
	vector<double> xAvSqrTotal(so.nSteps,0.0);
	vector<double> xSqrAvTheo(so.nSteps, 0.0);
	vector<double> xSqrAv(so.nSteps, 0.0);
	vector<double> xSqrAvTotal(so.nSteps, 0.0);
	vector<double> var(so.nSteps,0.0);
	//vector<double> xSqr1(so.nSteps, 0.0);
	//vector<double> xSqr2(so.nSteps, 0.0);
	//calcxAv(xAv,results.allX);
	if(so.noiseNr==1){
	calcxAvTheo(xAvTheo,so.dt,so.potw, so.gamma, so.x0, so.v0,so.nSteps);
	calcxSqrAvTheo(xSqrAvTheo,so.dt,so.potw, so.gamma,so.mass,so.D, so.x0, so.v0,so.nSteps);
	}
	else if(so.noiseNr==2 && so.corrFuncNr==0){
	  if(so.tau==2){
	calcxAvTheo1(xAvTheo,so.dt,so.x0,so.v0,so.nSteps);}
	  else{
	    calcxAvTheo2(xAvTheo,so.dt,so.x0,so.v0,so.nSteps);
	  }
	}
	//writeToFile(results.tVec,xAvTheo,filenames.xAvTheo,headerString,foldernames.main); //analytical solution for the inverse harmonical potential
        //writeToFile(results.tVec,xSqrAvTheo,filenames.xSqrAvTheo,headerString,foldernames.main);
	//vectorSquared(xSqr1,results.allX.at(0));
	 //vectorSquared(xSqr2,results.allX.at(1));
        buildAverage(xAv, results.allX);
	buildAverageSquared(xSqrAv, results.allX);
	ksim.xAvVec.at(j)= xAv;
	ksim.xSquaredAvVec.at(j)=xSqrAv;
	if(j==so.avNum-1)
	{
	buildAverage(xAvTotal,ksim.xAvVec);
	writeToFile(results.tVec,xAvTotal,filenames.xAvTotal,headerString,foldernames.main);
	buildAverage(xSqrAvTotal,ksim.xSquaredAvVec);
	writeToFile(results.tVec,xSqrAvTotal,filenames.xSqrAvTotal,headerString,foldernames.main);
	
	for(int n=0;n<xAvTotal.size();n++)
	{
	var.at(n)=xSqrAvTotal.at(n)-xAvTotal.at(n)*xAvTotal.at(n) ;
	}
	//writeToFile(results.tVec,var,filenames.xSqrAvTotal,headerString,foldernames.main);
	}
	//buildAverage(xAv, results.allX);
	writeToFile(results.tVec, results.allX.at(0), filenames.x, headerString, foldernames.main);
	writeToFile(results.tVec, results.allX.at(1), filenames.x, headerString, foldernames.main);
	writeToFile(results.tVec, results.allX.at(2), filenames.x, headerString, foldernames.main);
	//writeToFile(results.tVec, xAv, filenames.xAv, headerString, foldernames.main);
	//writeToFile(results.tVec, velocAv, filenames.velocAv, headerString, foldernames.main);
	//writeToFile(results.tVec, xSqrAv, filenames.xSqrAv, headerString, foldernames.main);
	//writeToFile(results.tVec, xSqr1, filenames.test1, headerString, foldernames.main);
	//writeToFile(results.tVec, xSqr2, filenames.test2, headerString, foldernames.main);

    //=====================================averaged velocity !!============================================================== 
	vector<double> velocAv(so.nSteps, 0.0);
	vector<double> vSqrAv(so.nSteps, 0.0);
	vector<double> velocAvTotal(so.nSteps,0.0);
	vector<double> vSqrAvTotal(so.nSteps, 0.0);
	buildAverage(velocAv, results.allV);
	buildAverageSquared(vSqrAv, results.allV);
	ksim.vAvVec.at(j)= velocAv;
	ksim.vSquaredAvVec.at(j)=vSqrAv;
	if(j==so.avNum-1)
	{
	buildAverage(velocAvTotal,ksim.vAvVec);
	writeToFile(results.tVec,velocAvTotal,filenames.velocAvTotal,headerString,foldernames.main);
	buildAverage(vSqrAvTotal,ksim.vSquaredAvVec);
	writeToFile(results.tVec,vSqrAvTotal,filenames.vSqrAvTotal,headerString,foldernames.main);
	}
	writeToFile(results.tVec,results.allV.at(0),filenames.veloc,headerString,foldernames.main);
	writeToFile(results.tVec,results.allV.at(1),filenames.veloc,headerString,foldernames.main);
	writeToFile(results.tVec,results.allV.at(2),filenames.veloc,headerString,foldernames.main);
	
	//====================================potential related=======================================
	
	//gnuplot-file für potenzial
	//============================
	vector <double> xVec(10000,0.0);
	vector <double> potVec(10000,0.0);
	
	//void GnuPotential(vector <double>& xVec, vector <double>& potVec, function<double(const double&, const double&)> potExtFunc,
	//const int& n, const double& x0, const double& xn)
	GnuPotential(xVec,potVec,potential.potentialFunc, 10000, 0.0, 20.0);
	
	if(j==0)
	{
	writeToFile(xVec,potVec,filenames.potential,headerString,foldernames.main);
	}
	
	
	
	if (so.potNr != 0)
	{
	  //potential energy
	  //================
	  vector<double> potEnergyAv(so.nSteps, 0.0);
	  calcPotEnergyAv(potEnergyAv, results.allX, results.tVec, potential.potentialFunc);
	  writeToFile(results.tVec, potEnergyAv, filenames.potEnergyAv, headerString,foldernames.main);

	  //total energy
	  //============
	  vector<double> totalEnergyAv(so.nSteps, 0.0);
	  add2Vec(totalEnergyAv, potEnergyAv, kinEnergyAv);
	  writeToFile(results.tVec, totalEnergyAv, filenames.totalEnergyAv, headerString, foldernames.main, 7);
	}

	//given potential
	//===============
	vector<double> potentialRange(so.eaPotentialN, 0.0);
	linspace(so.eaPotBorderL, (so.eaPotBorderR-so.eaPotBorderL)/so.eaPotentialN, potentialRange);
	vector<double> potentialVec(so.eaPotentialN, 0.0);
	for (int i = 0; i < so.eaPotentialN; i++)
	{
	  potentialVec.at(i) = potential.potentialFunc(0.0, potentialRange.at(i));
	}
	//writeToFile(potentialRange, potentialVec, filenames.potential, headerString, foldernames.main);

	//===================Theoretical curve of kinetic energy without potential=====================================
	if (2*so.D/so.mass/so.k_b/so.temperature*so.tau >= 1 && so.potNr == 0 && so.noiseNr == 2)
	{
	  vector<double> kinEnergyTheo(so.nSteps, 0.0);
	  calcKinEnergyTheo(kinEnergyTheo, results.tVec, so.temperature, so.mass, so.D, so.tau, so.k_b);
	  writeToFile(results.tVec, kinEnergyTheo, filenames.kinEnergyTheo, headerString, foldernames.main);
	}
	else if (so.potNr == 0 && so.noiseNr == 2)
	{
	  printf("did not calculate kinEnergyTheo as 2*Q*tau >= 1!\n");
	}

	//theoretical spectral density for no potential and colored noise and corrFuncNr 0
	if (so.potNr==0 && so.noiseNr==2 && so.corrFuncNr==0)
	{
	  int spectralN = 500;
	  double omegaMax = 12.0;
	  vector<double> specOmegaTheo(spectralN);
	  linspace(0.0, omegaMax/spectralN, specOmegaTheo);
	  vector<double> specVVecTheo(spectralN, 0.0);
	  calcSpecVTheo0(specVVecTheo, specOmegaTheo, so.k_b, so.temperature, so.mass, so.D, so.tau);
	  //writeToFile(specOmegaTheo, specVVecTheo, filenames.specVTheo, headerString, foldernames.main);
	}
	else if(so.potNr==1 && so.noiseNr == 2 && so.corrFuncNr == 0)
	{
	  int spectralN = 500;
	  double omegaMax = 12.0;
	  vector<double> specOmegaTheo(spectralN);
	  linspace(0.0, omegaMax/spectralN, specOmegaTheo);
	  vector<double> specXVecTheo(spectralN, 0.0);
	  calcSpecXTheo1(specXVecTheo, specOmegaTheo, so.k_b, so.temperature, so.mass, so.D, so.tau, so.potK);
	  //writeToFile(specOmegaTheo, specXVecTheo, filenames.specXTheo, headerString, foldernames.main);
	}

	//correlation of x when in equilibrium
	vector<double> settledTVec((int) round((so.tEnd-so.t0-so.timeSettled)/so.dt),0.0);
	linspace(0.0, so.dt, settledTVec);
	vector<double> corrX(settledTVec.size(), 0.0);
	vector<double> corrXTotal(settledTVec.size(), 0.0);

	buildCorrAv(corrX, results.allX, (int) round(so.timeSettled/so.dt), 0);
	ksim.corrXAv.at(j)=corrX;
	if(j==so.avNum-1)
	{
	buildAverage(corrXTotal,ksim.corrXAv);
	writeToFile(settledTVec,corrXTotal,filenames.corrXTotal,headerString,foldernames.main);
	}
	
// 	writeToFile(settledTVec, corrX, filenames.corrX, headerString, foldernames.main);

	//spectral density of x
	double dw;
	corrX.resize(pow(2,floor(log2(corrX.size()))));
	dw = M_PI/so.dt/corrX.size();
	double omegaMax = 12.0;
	vector<double> specXVec(corrX.size(), 0.0);
	vector<double> specOmega(ceil(omegaMax/dw));
	linspace(0.0, dw, specOmega);
	calcSpec(specXVec, corrX, so.dt);
	//writeToFile(specOmega, specXVec, filenames.specX, headerString, foldernames.main);

	//correlation of v when in equilibrium
	vector<double> corrV((int) round((so.tEnd-so.t0-so.timeSettled)/so.dt), 0.0);
	vector<double> corrVTotal(corrV.size(), 0.0);
	buildCorrAv(corrV, results.allV, (int) round(so.timeSettled/so.dt), 0);
	
	ksim.corrVAv.at(j)=corrV;
	if(j==so.avNum-1)
	{
	buildAverage(corrVTotal,ksim.corrVAv);
	writeToFile(settledTVec,corrVTotal,filenames.corrVTotal,headerString,foldernames.main);
	}
//  	writeToFile(settledTVec, corrV, filenames.corrV, headerString, foldernames.main);

	//spectral density of v
	corrV.resize(pow(2,floor(log2(corrV.size()))));
	dw = M_PI/so.dt/corrV.size();
	vector<double> specVVec(corrV.size(), 0.0);
	calcSpec(specVVec, corrV, so.dt);
// 	writeToFile(specOmega, specVVec, filenames.specV, headerString, foldernames.main);

	if (so.noiseNr == 1 && so.potNr == 0)
	{
	  //<Delta x^2>
	  vector<double> xSqrTheoVec(so.nSteps, 0.0);
	  calcXSqrTheo0(xSqrTheoVec, results.tVec, so.mass, so.gamma, so.k_b, so.temperature);
	  //writeToFile(results.tVec, xSqrTheoVec, filenames.xSqrTheo, headerString, foldernames.main);

	  //<1/2mv^2>
	  vector<double> kinEnergyTheo(so.nSteps, 0.0);
	  calcKinEnergyTheo0(kinEnergyTheo, results.tVec, so.mass, so.gamma, so.k_b, so.temperature, so.v0);
	  writeToFile(results.tVec, kinEnergyTheo, filenames.kinEnergyTheo, headerString, foldernames.main);
	}

	//=====================================energy animation==========================================
	if (so.eaBool)
	{
	  printf("Now creating files for energy animation...\n");
	  if(so.eaPotTimeDependent)
	  {
		prepareEnergyAnimationTimePot(results.allX, results.allV, results.tVec, so.mass, potential.potentialFunc, so.eaPotentialN,
			  so.eaPotBorderL, so.eaPotBorderR, so.eaStride, headerString, foldernames.ea);
	  }
	  else
	  {
		prepareEnergyAnimation(results.allX, results.allV, results.tVec, so.mass, potential.potentialFunc, so.eaStride, headerString, foldernames.ea);
	  }
	}
	
	// ===========================velocity distribution animation====================================
	if (so.vdaBool)
	{
	 printf("Now creating files for velocity distribution animation...\n");
	 vector <vector<double> > vDistTotal(ceil(double(so.nSteps)/double(so.vdaStride)),vector<double>(so.vDistNBins,0.0));
	 cout << vDistTotal.size () << " " << ksim.vDistAvVec.at(0).size()<< endl; // size entspricht gerade der Anzahl verschiedener Zeitpunkte
	 prepareVda(results.allV, vDistT,ksim.vDistAvVec.at(j), so.vDistRangeBeg, so.vDistRangeEnd, so.vDistNBins, so.vdaStride, headerString, foldernames.vda);
	if(j==so.avNum-1)
	  {
	  buildAverage3(vDistTotal,ksim.vDistAvVec); //hier wird über alle Verteilungen zu verschiedenen Zeitpunkten gemittelt,
	  for(int k=0; k<vDistTotal.size();k++)	     //nur sinnvoll für Zeiten, in denen das System bereits equilibriert ist
	  {					     //aktuell hat prepareVda() keine Aussagekraft
	   writeToFile(vDistT, vDistTotal.at(k), "vDistTotalNorm", headerString, foldernames.vda); 
	  }
	  }
	}
	if (so.pdaBool)
	{
	 printf("Now creating files for position distribution animation...\n");
	  vector <vector<double> > pDistTotal(ceil(double(so.nSteps)/double(so.pdaStride)),vector<double>(so.pDistNBins,0.0));
	  preparePda(results.allX, pDistT,ksim.pDistAvVec.at(j), so.pDistRangeBeg, so.pDistRangeEnd, so.pDistNBins, so.pdaStride, headerString, foldernames.pda);
	  //cout << "pDistTotal.size()" << " " <<pDistTotal.size() << " " << "pDistTotal.at(0).size()" << " " << pDistTotal.at(0).size() << endl;
	  //cout << "pDistAvVec.at(0).size()" << " " << ksim.pDistAvVec.at(0).size() << " " << "pDistAvVec.at(0).at(0).size()" << " " << ksim.pDistAvVec.at(0).at(0).size() << endl;
	  if(j==so.avNum-1)
	  {
	  buildAverage3(pDistTotal,ksim.pDistAvVec);
	  for(int k=0; k<pDistTotal.size();k++) 
	  {
	   writeToFile(pDistT, pDistTotal.at(k), "pDistTotalNorm", headerString, foldernames.pda); 
	  }
	  }
	}
}
