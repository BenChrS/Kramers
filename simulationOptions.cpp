#include <stdio.h>
#include <iostream>
#include <math.h>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "simulationOptions.h"

using namespace std;

 SimulationOptions::SimulationOptions(){
 	
//  	setInitValues();
//  	setDependentVariables();
 }
void SimulationOptions::readInput(int argc, char *argv[])
	{
	  
	  
	  
	  
	//input parameters	
	string temp_string;
	string inputFile  = argv[2];
	//string inputpath ="/home/schueller/Desktop/Kramers/build/"+inputFile ;
	// read from file:
	ifstream input;
 	input.open( inputFile.c_str());
// 	input.open("input.txt");
		getline(input,temp_string);
		stringstream(temp_string) >> this->tEnd;
		getline(input,temp_string);
		stringstream(temp_string) >> this->nStepsFactor;
		getline(input,temp_string);
		stringstream(temp_string) >> this->nStepsTwo;
// 		getline(input,temp_string);
// 		stringstream(temp_string) >> this->tau;
// 		getline(input,temp_string);
// 		stringstream(temp_string) >> this->a;
// 		getline(input,temp_string);
// 		stringstream(temp_string) >> this->chi;
 		getline(input,temp_string);
 		stringstream(temp_string) >> this->potw;
//   		getline(input,temp_string);
//  		stringstream(temp_string) >> this->D;
 		getline(input,temp_string);
 		stringstream(temp_string) >> this->temperature;
// 		getline(input,temp_string);
// 		stringstream(temp_string) >> this->B_eff;
		getline(input,temp_string);
		stringstream(temp_string) >> this->npTen;
		getline(input,temp_string);
		stringstream(temp_string) >> this->npTwo;
		getline(input,temp_string);
		stringstream(temp_string) >> this->avOpt;
		getline(input,temp_string);
		stringstream(temp_string) >> this->avNum;
// 		getline(input,temp_string);
// 		stringstream(temp_string) >> this->transition;
// 		getline(input,temp_string);
// 		stringstream(temp_string) >> this->npoints;
		getline(input,temp_string);
		stringstream(temp_string) >> this->edaBool;
		getline(input,temp_string);
		stringstream(temp_string) >> this->pdaBool;
		getline(input,temp_string);
		stringstream(temp_string) >> this->vdaBool;
		getline(input,temp_string);
		stringstream(temp_string) >> this->corrFuncNr;
		getline(input,temp_string);
		stringstream(temp_string) >> this->potNr;
		getline(input,temp_string);
		stringstream(temp_string) >> this->noiseNr;
		getline(input,temp_string);
		stringstream(temp_string) >> this->sdeSolverNr;
		getline(input,temp_string);
		stringstream(temp_string) >> this->Kramers;
		getline(input,temp_string);
		stringstream(temp_string) >> this->initCondNr;
		getline(input,temp_string);
		stringstream(temp_string) >> this->xc;
		getline(input,temp_string);
		stringstream(temp_string) >> this->Ub;
		getline(input,temp_string);
		stringstream(temp_string) >> this->rxborder;

// 		getline(input,temp_string);
// 		stringstream(temp_string) >> this->testBool;
// 		getline(input,temp_string);
// 		stringstream(temp_string) >> this->a_Corr;
// 		getline(input,temp_string);
// 		stringstream(temp_string) >> this->b_Corr;
// 		getline(input,temp_string);
		//stringstream(temp_string) >> this->pota;
		//getline(input,temp_string);
		//stringstream(temp_string) >> this->potMy;
		//getline(input,temp_string);
		//stringstream(temp_string) >> this->potNy;
	input.close();
	
	setInitValues(argc,argv);
	setDependentVariables();  
	}
void SimulationOptions::setInitValues(int argc, char *argv[]){
	this->commentLine = "Standard_5.0";
        
	//physical settings
	this->k_b = 1; // k_boltzmann constant
	//this->temperature = 1.0; //2.0
	//this->mass =40.0/(25.0*6*6);//1.0;//0.1;//0.1*this->temperature;//1.0/16.0*this->temperature; //1.0; //mass of particle
	string D_str =argv[3];
	this->D =atof( D_str.c_str() );
	string corrTime=argv[4];
        this->tau = atof( corrTime.c_str() ); // see paper /only important for first correlation function
	this->a =atof( corrTime.c_str() ); // only important for second correlation function
	this->chi=atof( corrTime.c_str() );//correlation time for third correlation function
//  	this->alpha=30.0*sqrt(this->mass); // correlation time for massless theory // Einheit Masse/Zeit alpha=sqrt(mass)*alpha' 
			//alpha=sqrt(mass) für alpha'=1 , alpha' ist inverse Korrelationszeit
	this->minAlpha=-this->alpha*this->alpha/4.0*exp(-2.0);  //Minimum des Kernels
	//cout << "minAlpha " << minAlpha << endl;

	//statistical/program settings
 	this->t0 = 0.0; //time interval [t0, t1]
        //this->t1 = 10.0;     //1.0: für kb*T=0.5 Limes
	this->tSettling = 10.0; // time needed for I(t) to be approximately 0
	//this->timeSettled = (this->t1-this->t0)/3.0; //approximate time particles need to be in equilibrium - only important for kinetic Energy Average - not yet in external call
        //this->nStepsFactor = 1;//round(this->t1-this->t0);
        //this->nStepsTwo =15;   //15: für kb*T=0.5 Limes
	this->nSteps =nStepsFactor*pow(2, this->nStepsTwo); ///2000; // number of final datapos (stored), must be devidable by 2
        //this->npTen=4;
        //this->npTwo=0;
	this->np = pow(10,npTen)*pow(2,this->npTwo); //number of averaged simulations (number of particles)
	this->optimisationBool = true;  //must be set to false for third correlation function
	this->gPrecisionTen = -3;
	this->gPrecision = pow(10,this->gPrecisionTen); //if opitmisationBool is set to true gPrecision will be the smallest value in G(t)
	this->dPrecisionTen = -3;
	this->dPrecision = pow(10,this->dPrecisionTen); //if optimisationBool is set to true dPrecision will be the smallest value in dispKernel(t)
	this->maxThreadNumber = 24; //must be greater than actual used number of threads d
	this->maximumAllacationMemory = 60; //program does not run if calculated allocation memory exceeds this value
	this->saveAllDataBool = false;
	this->testBool = false;
	this-> paperBool = 0;
	this->a_Corr = 6; 
	this->b_Corr = 7;
	//this->avOpt = false;
	//this->avNum = 1;

	//this->corrFuncNr = 1;
        //this->potNr = 1;
        //this->noiseNr = 2;
        //this->sdeSolverNr = 3;
        //this->Kramers = 0;

	//potential Settings
	//this->potw= 5.0; //frequency of the harmonic potential
	this->potK = pow(this->potw,2.0)*this->mass; //for potential1
   	this->potDepth = 2; //for potenital2
	this->potWidth = 4.0; //for potential2
	this->potA = 0.5; //for potential 3
	this->potB = 2.0; //for potential 3
	this->potC = -2.0; //for potential 3
	this->pota = 1; //for potential 4 distance between minimum and maximum
	this->potMy = 2.5; //for potential 4 potential height
	this->potNy = -5; //for potential 4 height of right potential well    
	this->potMax=(3*this->pota*this->potNy)/(8*(2*this->potMy-this->potNy)); // for potential 4 : x-wert des Maximums 
	
#ifdef KRAMERS
 	cout << "max: " << this->potMax << endl;
	cout << "A" << " " << -(-2*potMy+potNy)/(2*pow(pota,4.0)) << endl;
        cout << "B" << " " << -(potNy)/(4*pow(pota,3.0)) << endl;
        cout << "C" << " " << -(2*potMy-potNy)/(pow(pota,2.0)) << endl;
        cout << "D" << " " <<  (3*potNy)/(4*pota) << endl;
        cout << "E" << " " <<  potMy << endl;
#endif
	this->potStartDepth = 5; //for potential 5
	this->potEndDepth = 0.5; //for potential 5
	this->potStartTime = 1.0; //for potential 5
	this->potEndTime = 5.0; //for potential 5
	
	//this->xc=1;
	//this->t1 =1.0; //1.0*this->xc; 
	this->xb=1.6*this->xc;
	//this->Ub=2.5;   // units of wc   //2.5*this->k_b*this->temperature;
	if(this->potNr==7)
	{
	this->mass =4.0*this->Ub/(this->potw*this->potw*pow(0.6*this->xc,2.0));
	}
	else
	{
	 this->mass=0.1;
	}
	this->potK = pow(this->potw,2.0)*this->mass; //for potential1
	this->alpha=30*sqrt(this->mass); // correlation time for massless theory // Einheit Masse/Zeit alpha=sqrt(mass)*alpha'  //a=30 für EkinAv
			//alpha=sqrt(mass) für alpha'=1 , alpha' ist inverse Korrelationszeit
	this->timeSettled = (this->tEnd-this->t0)/3.0; //approximate time particles need to be in equilibrium - only impo		
	this->minAlpha=-this->alpha*this->alpha/4.0*exp(-2.0);  //Minimum des Kernels
// 	cout << "minAlpha " << minAlpha << endl;
	
	
// 	this->mass = this->Ub/2.0; // Test für Skalierungsverhalten Ub/m=2
	
// 	cout << "omegaB: " << sqrt(4*this->Ub/(this->mass*pow((this->xb-this->xc),2.0))) << endl;
// 	cout << "m*omegaB: " << this->mass*pow(sqrt(4*this->Ub/(this->mass*pow((this->xb-this->xc),2.0))),2.0) << endl;
	
// 	this->D = 2*this->k_b*this->temperature*this->mass*sqrt(4*this->Ub/(this->mass*pow((this->xb-this->xc),2.0)));   //nur zum check des Skalierungsverhalten von Kramersrate, Wieder auskommentieren später
// 	cout << "D: " << 2*this->k_b*this->temperature*this->mass*sqrt(4*this->Ub/(this->mass*pow((this->xb-this->xc),2.0))) << endl;

	//energy animation (ea) settings
	this->eaBool = false; // energy Animation?
	this->eaPotTimeDependent = false; //if true a potential file is generated for each time step of the animation
	this->fps = 1; //files per second
	this->eaStride = max(1,(int) floor(this->nSteps/(this->tEnd-this->t0)/this->fps)); //number of frames skipped for the energy Animation - if set to 1 a file is generated for all time steps
	this->eaPotentialN = 100;
	this->eaPotBorderL = -8.0;
	this->eaPotBorderR = 8.0;

 //transition options
        this->transition=false;
	this->npoints=2;	
	
    //Random number generator for initial velocity
    gsl_rng_env_setup();
    this->T= gsl_rng_default;
    this->res_rng=gsl_rng_alloc(T);
    this->seed=time(NULL);
    string seed1_str=argv[5];
    long seed1_lon=atol(seed1_str.c_str());
    this->seed1=time(NULL)+seed1_lon;
	

	//initial conditions (x0, v0)
   
    
    // wähle Anfangsbedingungen für Vergleich mit paper oder beliebig
    if(this->paperBool==0)
    {
      if(this->potNr==7)
      {
      x0=1.0*this->xc;//0.0;
      }
      else
      {
      x0=0.0;
      }
      v0=0.0;
      
      if(this->noiseNr==2)
      {
	//this->potw= 5.0;//sqrt(4*this->Ub/(this->mass*pow((this->xb-this->xc),2.0)));
	this->gamma = this->D/(2.0*this->k_b*this->temperature);
	double time;
	double beta,om;
	double p,q,diskr;
	beta=this->gamma/this->mass;
	time=this->tau;
	om = this->potw;
	 
	p = beta/time-om*om-1.0/(3.0*time*time);
	q = 2.0/(27.0*pow(time,3.0))-1.0/(3.0*time)*(beta/time-om*om)-om*om/time;
	
	//cout << "p: " << p << " " << "q: " << q << endl;
	 
	//cout << "tau: " << time << " " << "wc: " << this->potw << " " << "gamma: " << this->gamma << " " << "beta: " 
	//     << beta <<  endl;
	 
	 //Eigenwerte für inverses harmonisches Potenzial
	 double a,ny,my,d,u,v;
	 a=1.0/(3.0*time);
	 //cout << "a" << " " <<a << endl;
	 ny=-1.0+3.0*beta*time-3.0*pow((this->potw*time),2.0);
// 	 cout << "ny" << " " <<ny << endl;
	 my=-1.0+9.0/2.0*beta*time+9.0*pow((this->potw*time),2.0);
// 	 cout << "my" << " " <<my << endl;
	 
	 diskr =  my*my+pow(ny,3.0);
// 	 cout << "diskr: " << diskr << endl;
	 d=sqrt(pow(my,2.0)+pow(ny,3.0));
// 	 cout << "d" << " " <<d << endl;
	 u=cbrt(my+d);
// 	 cout << u << " " << v << endl;
	 v=cbrt(my-d);
	 
	 if(diskr>0)
	 {
	  this->aColoured=-a+a*u+a*v; //eingevalue for coloured noise 
	  this->bColoured=-1.0/2.0*(1/time+this->aColoured);
	  this->cColoured=-1.0/2.0*(1/time+this->aColoured);
// 	  cout << "l1 " << this->aColoured << endl;
// 	  cout << "l2 " << this->bColoured << endl;
// 	  cout << "l3 " << this->cColoured << endl;
	 }
	 else if(diskr==0)
	 {
	  this->aColoured = 3.0*q/p - a;
	  this->bColoured = -3.0*q/(2.0*p)-a;
	  this->cColoured = -3.0*q/(2.0*p)-a;
// 	  cout << "l1 " << this->aColoured << endl;
// 	  cout << "l2 " << this->bColoured << endl;
// 	  cout << "l3 " << this->cColoured << endl;
	 }
	 else
	 {
	  this->aColoured = sqrt(-4.0/3.0*p)*cos(1.0/3.0*acos(-q/2.0*sqrt(-27.0/pow(p,3.0))))-a;
	  this->bColoured = -sqrt(-4.0/3.0*p)*cos(1.0/3.0*acos(-q/2.0*sqrt(-27.0/pow(p,3.0)))+M_PI/3.0)-a;
	  this->cColoured = -sqrt(-4.0/3.0*p)*cos(1.0/3.0*acos(-q/2.0*sqrt(-27.0/pow(p,3.0)))-M_PI/3.0)-a;
// 	  cout << "l1 " << this->aColoured << endl;
// 	  cout << "l2 " << this->bColoured << endl;
// 	  cout << "l3 " << this->cColoured << endl;
	   
	}
      }
      else
      {
	//this->potw= 5.0;//sqrt(4*this->Ub/(this->mass*pow((this->xb-this->xc),2.0)));
	this->gamma = this->D/(2.0*this->k_b*this->temperature);
	double beta;
	beta=this->gamma/this->mass;
	this->aWhite=0.5*(sqrt(pow((beta),2)+4*pow(this->potw,2))-beta); //eigenvalue for white noise
      }
    }
    else
    { //Berechnung der Anfangsgeschwindigkeit über die entsprechenden Eigenwerte 
	if(this->potNr==1 & this->potK<0 & this->transition==false)
	{
// 	  cout << potK  << " " << "..." << endl;
	  this->x0 = -1;	
	  this->potB1 = 0.5*this->mass*pow(this->potw,2)*pow(this->x0,2); // height of the inverse harmonic function to be overcome by the particle
	  this->temperature = this->potB1/2; // setze Temperatur mit Höhe des Potenzials in Verbindung um zu gewährleisten, dass T auf jeden Fall kleiner ist als Schwelle !
	  this->gamma = this->D/(2.0*this->k_b*this->temperature);
	  if(this->noiseNr==1)
	  {
	    this->aWhite=0.5*(sqrt(pow((this->gamma),2)+4*pow(this->potw,2))-this->gamma); //eigenvalue for white noise
	    this->potB_eff = pow((this->potw/this->aWhite),2)*this->potB1;
	  }
	  
	  else if(this->noiseNr==2)
	  {
// // 	    cout << "noiseNr==2" << endl;
	    double time;
	  
	    time=this->tau;
	      
	      
	      //cout << "time" << " " <<time << endl;
	    double a,ny,my,d,u,v;
	    a=1.0/(3.0*time);
	    //cout << "a" << " " <<a << endl;
	    ny=-1.0+3.0*this->gamma*time-3.0*pow((this->potw*time),2);
	    //cout << "b" << " " <<b << endl;
	    my=-1.0+9.0/2.0*this->gamma*time+9.0*pow((this->potw*time),2);
	    //cout << "c" << " " <<c << endl;
	    d=sqrt(pow((my),2)+pow((ny),3));
	    //cout << "d" << " " <<d << endl;
	    u=pow((my+d),1.0/3.0);
	    v=pow((my-d),1.0/3.0);
	    this->aColoured=-a+a*u+a*v; //eingevalue for coloured noise 
	    this->bColoured=-1/2*(1/time+this->aColoured);
	    this->cColoured=-1/2*(1/time+this->aColoured);
	    //cout << "aColoured" << " " << aColoured<< endl;
	    this->potB_eff = pow((this->potw/this->aColoured),2)*this->potB1;
	  }
	  switch(this->B_eff)
	  {
			       case 0:{
				  this->v0=0; }break;
			       case 1:{
				  this->v0 = sqrt(this->potB_eff/this->mass);}break;
			       case 2:{
				  this->v0 = sqrt(2*this->potB_eff/this->mass);}break;
			       case 3:{
			          this->v0 = sqrt(4*this->potB_eff/this->mass);}break;
	  }

	}
    }
// 	cout << " v0:" << " " << this->v0 << endl;
	
	//position distribution animation (pda) settings
	//this->pdaBool = false; 
	this->pdaStride = max(1,(int) floor(this->nSteps/(this->tEnd-t0)/this->fps)); /* entweder für jeden einzelnen Zeitschritt
	 (d.h. pdaStride=1) oder nur einen Bruchteil davon abhängig von der Wahl von fps */
	this->pDistRangeBeg = -4;//-10*sqrt(this->k_b*this->temperature/(this->mass*this->potw*this->potw));//-10.0;
	this->pDistRangeEnd = -this->pDistRangeBeg;
	this->pDistNBins = max(2,(int)round(this->np/40.0));
	this->pDistStartTime = (tEnd-t0)/3.0;
	
	//velocity distribution animation (vda) settings
	//this->vdaBool = false; 
	this->vdaStride = max(1,(int) floor(this->nSteps/(this->tEnd-t0)/this->fps)); // bestimmt zusammen mit nSteps die Anzahl verschiedener Zeitpunkte, für die Verteilung erzeugt wird
	//theoretical velocity distribution settings
	this->vDistRangeBeg = -15;//-10*sqrt(this->k_b*this->temperature/this->mass);
	this->vDistRangeEnd = -this->vDistRangeBeg;
	this->vDistNBins = max(2,(int)round(this->np/40.0)); // es sollen mindestens 40 Teilchen in einem Bin enthalten sein 
	this->vDistStartTime = (tEnd-t0)/3.0;
	
	//kinEnergy distribution animation (eda) settings
	//this->edaBool = true; 
	this->edaStride = max(1,(int) floor(this->nSteps/(this->tEnd-t0)/this->fps)); // bestimmt zusammen mit nSteps die Anzahl verschiedener Zeitpunkte, für die Verteilung erzeugt wird
	//theoretical velocity distribution settings
	this->eDistRangeBeg = 0;//-10*sqrt(this->k_b*this->temperature/this->mass);
	this->eDistRangeEnd = 30;
	this->eDistNBins = max(2,(int)round(this->np/40.0)); // es sollen mindestens 40 Teilchen in einem Bin enthalten sein 
	this->eDistStartTime = (tEnd-t0)/3.0;

	this->xCutoff = -1.0; //if absolute position of particles
					// becomes greater than this value,
					// they will stay at this position.
					// Set to negative value for not making use of this option
	//Kramers
	switch(this->Kramers)
	{
	  case 0:
	 {
	   this->KramersRate = false;
	 }break;
	  case 1:
	 {
	   this->KramersRate = true;
       //this->dx=0.0;
  //         this->rxborder =2.6*this->xc; //1.7 tail-Paper;//sqrt(2*this->potB_eff/abs(this->potK))+this->dx;
	   this->lxborder = -10000.0;           
	   if(this->potNr==1 & this->potK==-1)
	   {
              this->lxborder = -2.0;//-(this->rxborder+2*this->rxborder-this->dx);
	   }
	 }break;
	}
//	this->rxborder=2.6*this->xc;

}



void SimulationOptions::setDependentVariables(){
this->dt = (tEnd-t0)/((double)nSteps);
 //this->dt = (10*xc)/((double)nSteps); //mit passendem Skalierungsverhalten für Pot 7
// cout << "dt " << dt << " mass " << mass << endl;

  //----prepare simulations----

this->gamma = this->D/(2.0*this->k_b*this->temperature);
// cout << "gamma/m: " << this->gamma/this->mass << endl;
this->nSettling = pow(2,ceil(log2(this->tSettling/this->dt))); //must be at least as high as I(t) needs to be at approximately 0
// cout <<" nSettling " << this->nSettling << endl; 
this->nFourier =this->nSettling*pow(2, 1); // must be larger than nSettling to avoid boundary
					//  effects of the FFT and be a power of 2 for the FFT
// cout << "nFourier " << this->nFourier << endl;
}



