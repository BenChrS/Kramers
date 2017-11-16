#include <stdio.h>
#include <iostream>
#include <stdio.h>

#include <tr1/functional>

#include <iterator>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>

#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>

#include <omp.h>

#include "simulationOptions.h"
#include "sdeSolver.h"
#include "noiseDiss.h"
//-----------------------------
#include <fstream>
//--------------------------------
using namespace std;
using namespace std::tr1;
using namespace std::tr1::placeholders;
	
// Adams-Bashforth scheme for white noise
//solves the Langevin equation with white noise and friction coefficient gamma
int adamsBWhiteNoise(const vector<double>& tVec, const double& dt, const double& x0, const double& v0, const double& mass,const double& potw,
		     const double& k_b, const double& temperature,const vector<double>& whiteNoise, const double gamma, 
		     const double D,vector<double>& xVec, vector<double>& velocVec,function<double(const double&, const double&)> forceExtFunc,
		     const bool& KramersRate,const double rxborder,const double lxborder,gsl_rng* res_rng,gsl_rng* noise_rng,const int& initCond,
		     const NoiseDiss& noiseDiss)
{
  int nSteps = xVec.size();
  double force1 = 0.0;
  double force2 = 0.0;
  double force3 = 0.0;
  int stepZero=0;
  int end = 0;
  vector <double> Noisecp;
  Noisecp.resize(whiteNoise.size(),0.0);
  Noisecp = whiteNoise;

//   fstream Noise, force;
//   force.open("force.dat",ios::out);
//   Noise.open("Noise.dat", ios::app);
//   
//    for(int i = 0; i < Noisecp.size(); i++)
// 	      {
// 		Noise << tVec.at(i) << "     "  << whiteNoise.at(i)  << endl;
// 	      }
  
  // --steps for t = dt, Euler-Cromer Method--
  xVec.at(0) = x0;
  velocVec.at(0) = v0;

  force1 = forceExtFunc(tVec.at(0), xVec.at(0)) - gamma*velocVec.at(0);
//   force1 = mass*forceExtFunc(tVec.at(0), xVec.at(0)) - gamma*velocVec.at(0);

  // Euler-Cromer is more stable for harmonic oscillations, the velocity is evoluated first and is used for x in the same step
  velocVec.at(1) = velocVec.at(0) + dt/mass*(sqrt(D) * Noisecp.at(0) + force1);
  xVec.at(1) = xVec.at(0) + dt*velocVec.at(1);

  force2 = force1;

  // --steps for t = 2dt, Adams-Bashforth n=2--

  force1 = forceExtFunc(tVec.at(1), xVec.at(1)) - gamma*velocVec.at(1);
//   force1 = mass*forceExtFunc(tVec.at(1), xVec.at(1)) - gamma*velocVec.at(1);


  xVec.at(2) = xVec.at(1) + dt*(1.5*velocVec.at(1)-0.5*velocVec.at(0));
  velocVec.at(2) = velocVec.at(1) + dt/mass*(1.5*(sqrt(D) * Noisecp.at(1) + force1)
						-0.5*(sqrt(D) * Noisecp.at(0) + force2));

  force3 = force2;
  force2 = force1;

  // --steps for t > 2dt, Adams-Bashforth n=3--
  for (int step = 3; step < nSteps; step++)
  {
    if(KramersRate == true)
	  {
	    while(xVec.at(step-1)>=rxborder || xVec.at(step-1)<=lxborder)
	    {
	     stepZero = step;
             switch(initCond)
               {
                case 0:
                       {
                        xVec.at(step) = x0;
                        velocVec.at(step) = v0;
                       }break;

                case 1:
                       {
                        xVec.at(step) = x0;
                        velocVec.at(step) = gsl_ran_gaussian(res_rng, sqrt(k_b*temperature/mass));//random initial velocity
			// könnte problematisch sein ?!
                         //cout << velocVec.at(step) << endl ;
                       }break;

                case 3:
                       {
                        xVec.at(step) = gsl_ran_gaussian(res_rng, sqrt(k_b*temperature/(mass*potw*potw)));//random initial velocity;
                        velocVec.at(step) = gsl_ran_gaussian(res_rng, sqrt(k_b*temperature/mass));//random initial velocity;
                       }break;


                }

              //fstream Noise1;
	      //Noise1.open("Noise1.dat", ios::app);
             
             //noise
             noiseDiss.noiseFunc(dt, noise_rng, Noisecp); // könnte problematisch sein ?!
	     
// 	     for(int i = 0; i < Noisecp.size(); i++)
// 	      {
// 		Noise1 <<tVec.at(i) << "          " << Noisecp.at(i) << endl;
// 	      } 
	     if(step==nSteps-1)
	     {
	       end = 1;
	       break;
	     }
	   
	      force1 = 0.0;
	      force2 = 0.0;
	      force3 = 0.0;
	    
	      force1 = forceExtFunc(tVec.at(step), xVec.at(step))-gamma*velocVec.at(step);
// 	      force1 = mass*forceExtFunc(tVec.at(step), xVec.at(step))-gamma*velocVec.at(step);
	    
	      velocVec.at(step+1) = velocVec.at(step) + dt/mass*(sqrt(D) * Noisecp.at(step-stepZero) + force1);
	      xVec.at(step+1) = xVec.at(step) + dt*velocVec.at(step+1);
	    
	      if(step==nSteps-2)
	      {
	       end = 1;
	       break; 
	      }
	    
	      force2 = force1;

	      // --steps for t = 2dt, Adams-Bashforth n=2--
	      force1 = forceExtFunc(tVec.at(step+1), xVec.at(step+1))- gamma*velocVec.at(step+1);
// 	      force1 = mass*forceExtFunc(tVec.at(step+1), xVec.at(step+1))- gamma*velocVec.at(step+1);

	      xVec.at(step+2) = xVec.at(step+1) + dt*(1.5*velocVec.at(step+1)-0.5*velocVec.at(step));
	      velocVec.at(step+2) = velocVec.at(step+1) + dt/mass*(1.5*(sqrt(D) * Noisecp.at(step+1-stepZero) + force1)
						 -0.5*(sqrt(D) * Noisecp.at(step-stepZero) + force2));

	      force3 = force2;
	      force2 = force1;
	    
	      if(step==nSteps-3)
	       {
	         end = 1;
	         break; 
	       }
	       step += 3;
	    }
	  }
	  
	  if(end == 1)
	  {
	   break; 
	  }

          force1 = forceExtFunc(tVec.at(step-1), xVec.at(step-1)) - gamma*velocVec.at(step-1);
// 	  force1 = mass*forceExtFunc(tVec.at(step-1), xVec.at(step-1)) - gamma*velocVec.at(step-1);

          xVec.at(step) = xVec.at(step-1) + dt*(23.0/12.0*velocVec.at(step-1)
							  -16.0/12.0*velocVec.at(step-2)
							  +5.0/12.0*velocVec.at(step-3));
          velocVec.at(step) = velocVec.at(step-1)
				    + dt/mass*(23.0/12.0*(sqrt(D) * Noisecp.at(step-1-stepZero) + force1)
					      -16.0/12.0*(sqrt(D) * Noisecp.at(step-2-stepZero) + force2)
					      +5.0/12.0*(sqrt(D) * Noisecp.at(step-3-stepZero) + force3));
         force3 = force2;
         force2 = force1;
  }
//   for(int i = 0; i < Noisecp.size(); i++)
// 	      {
// 		force << tVec.at(i) << "     "  << forceExtFunc(tVec.at(i), xVec.at(i))  << endl;
// 	      }
	      
// 	    cout << velocVec.at(0) + dt/mass*(sqrt(D) * whiteNoise.at(0)+forceExtFunc(tVec.at(0),xVec.at(0))-gamma*velocVec.at(0) ) << endl ;
// 	    cout << xVec.at(0)+dt*velocVec.at(1) << endl;
// 	    cout << dt  << " " << xVec.at(0) << " " << velocVec.at(0) << " " <<  D << endl;
  return 1;
}

// Euler scheme for white noise
//solves the Langevin equation with white noise and friction coefficient gamma with the euler method
int eulerWhiteNoise(const vector<double>& tVec, const double& dt, const double& x0, const double& v0, const double& mass,
		    const double& k_b, const double& temperature,const vector<double>& whiteNoise, const double& gamma, 
		    const double& D,vector<double>& xVec, vector<double>& velocVec,function<double(const double&, const double&)> forceExtFunc,
		    const bool& KramersRate,const double rxborder,const double lxborder,gsl_rng* res_rng,gsl_rng* noise_rng,const int& initCond,
		    const NoiseDiss& noiseDiss)
{
  int nSteps = xVec.size();
  double force;
  int stepZero=0;
  int end = 0;
  vector <double> Noisecp;
  Noisecp.resize(whiteNoise.size(),0.0);
  Noisecp = whiteNoise;

  xVec.at(0) = x0;
  velocVec.at(0) = v0;

  for (int step = 1; step < nSteps; step++)
  {
      if(KramersRate == true)
      {	  
	  while(xVec.at(step-1)>=rxborder || xVec.at(step-1)<=lxborder)
	  {
          stepZero = step;
          double force = 0.0;
          switch(initCond)
          {
            case 0:
            {
             xVec.at(step) = x0;
             velocVec.at(step) = v0;
            }break;

            case 1:
            {
             xVec.at(step) = x0;
             velocVec.at(step) = gsl_ran_gaussian(res_rng, sqrt(k_b*temperature/mass));      //random initial velocity
            }break;
          }

          //noise
          noiseDiss.noiseFunc(dt, noise_rng, Noisecp);
	    if(step==nSteps-1)
	    { 
	      end = 1;
	      break;
	    }
	   step++;
	  }
	 }
	 if(end == 1)
	 {
	   break;
	 }
      force = forceExtFunc(tVec.at(step-1), xVec.at(step-1)) - gamma*velocVec.at(step-1);
      xVec.at(step) = xVec.at(step-1) + dt*velocVec.at(step-1);
      velocVec.at(step) = velocVec.at(step-1) + dt/mass*(sqrt(D) * Noisecp.at(step-1-stepZero) + force);
  }
  return 1;
}

// Euler scheme for colored noise
// solves the Langevin-differential equation with colored Noise and dissipation kernel with the euler method
// testing
int eulerColNoise(const vector<double>& tVec, const double& dt, const double& x0, const double& v0, const double& mass,
		  const double& k_b, const double& temperature, const vector<double>& colNoise, 
		  const vector<double>& dissipationKernel, vector<double>& xVec, vector<double>& velocVec,
		  function<double(const double&, const double&)> forceExtFunc,const bool& KramersRate,const double rxborder,
                  const double lxborder,gsl_rng* res_rng,gsl_rng* noise_rng,const int& initCond,const NoiseDiss& noiseDiss)
{
    int nSteps = xVec.size();
    int maxDispN = dissipationKernel.size();
    int lastDispN = 0;
    double dispIntegral = 0.0;
    double force = 0.0;
    int stepZero = 0;
    int end = 0;
    vector <double> colNoisecp;
    colNoisecp.resize(colNoise.size(),0.0);
    colNoisecp = colNoise;
    /*fstream Noise,Noise1;
    Noise.open("Noise.dat", ios::out);
    Noise1.open("Noise1.dat", ios::out);
    for(int i = 0; i < colNoisecp.size(); i++)
    {
    Noise << colNoise.at(i) << "          " << colNoisecp.at(i) << endl;
    }*/

    // --steps for t = dt, Euler-Cromer Method--
    xVec.at(0) = x0;
    velocVec.at(0) = v0;

    for(int step = 1; step < nSteps; step++)
    {
       if(KramersRate == true)
      {	  
	  while(xVec.at(step-1)>=rxborder || xVec.at(step-1)<=lxborder)
	  {
          stepZero = step;
          double dispIntegral = 0.0;
          double force = 0.0;
          //long seed = noiseDiss.seed;  //fixed seed!!
          switch(initCond)
          {
            case 0:
            {
             xVec.at(step) = x0;
             velocVec.at(step) = v0;
            }break;

            case 1:
            {
             xVec.at(step) = x0;
             velocVec.at(step) = gsl_ran_gaussian(res_rng, sqrt(k_b*temperature/mass));      //random initial velocity
             //cout << velocVec.at(step) << endl ;
            }break;
          }

          //noise
          noiseDiss.noiseFunc(dt, noise_rng, colNoisecp);
        /*for(int i = 0; i < colNoisecp.size(); i++)
            {
            Noise1 << colNoise.at(i) << "          " << colNoisecp.at(i) << endl;
        }*/
	    if(step==nSteps-1)
	    { 
	      end = 1;
	      break;
	    }
	   step++;
	  }
	 }
	 if(end == 1)
	 {
	   break;
	 }
      lastDispN = min(maxDispN, step-stepZero);
      //dispIntegral using trapezoid integral method
      dispIntegral = 0.5*(dissipationKernel.at(lastDispN-1)*velocVec.at(step-lastDispN)+   //Ränder
				dissipationKernel.at(0)*velocVec.at(step-1));
      for (int i = 1; i < lastDispN-1 ; i++)
      {
	dispIntegral += dissipationKernel.at(i)*velocVec.at(step-1-i);
      }

      force = forceExtFunc(tVec.at(step-1), xVec.at(step-1));

      velocVec.at(step) = velocVec.at(step-1) + dt/mass*(-2.0*dispIntegral*dt + colNoisecp.at(step-1-stepZero) + force);
      xVec.at(step) = xVec.at(step-1) + dt*velocVec.at(step-1);
    }
    return 1;
}
    
// Euler-Cromer scheme for colored noise, more stable for harmonic oscillations
// solves the Langevin-differential equation with colored Noise and dissipation kernel with the Euler-Cromer method
int eulerCromerColNoise(const vector<double>& tVec, const double& dt, const double& x0, const double& v0, const double& mass,
			const double& k_b, const double& temperature, const vector<double>& colNoise, const vector<double>& dissipationKernel,
		        vector<double>& xVec, vector<double>& velocVec,function<double(const double&, const double&)> forceExtFunc,
                        const bool& KramersRate,const double rxborder,const double lxborder,gsl_rng* res_rng,gsl_rng* noise_rng,
			const int& initCond,const NoiseDiss& noiseDiss)
{
    int nSteps = xVec.size();
    int stepZero = 0;
    int end = 0;
    int maxDispN = dissipationKernel.size();
    int lastDispN = 0;
    double dispIntegral = 0.0;
    double force = 0.0;
    vector <double> colNoisecp;
    colNoisecp.resize(colNoise.size(),0.0);
    colNoisecp = colNoise;

    // --steps for t = dt, Euler-Cromer Method--
    xVec.at(0) = x0;
    velocVec.at(0) = v0;
    //------------------------
    /*fstream Noise,Noise1;
    Noise.open("Noise.dat", ios::out);
    Noise1.open("Noise1.dat", ios::out);
    for(int i = 0; i < colNoisecp.size(); i++)
    {
    Noise << colNoise.at(i) << "          " << colNoisecp.at(i) << endl;
    }*/
    /*fstream diss,x,v,vO,xO,f,colour,colour1;
    diss.open("dissInt.dat", ios::out);
    x.open("xVec.dat", ios::out);
    v.open("vVec.dat", ios::out);
    vO.open("v0.dat" , ios::out);
    xO.open("x0.dat" , ios::out); 
    f.open("force.dat", ios::out);
    colour.open("Noise.dat", ios::out);
    colour1.open("Noise1.dat", ios::out);
    vO << v0 << " " << velocVec.at(0) << " " << vVector.at(0) << endl;
    for(int i = 0; i < colNoise.size(); i++)
    {
     colour << colNoise.at(i)  << endl;
    }*/
        //-------------------------------

    for(int step = 1; step < nSteps; step++)
    {
      if(KramersRate == true)
      {	  
	  while(xVec.at(step-1)>=rxborder || xVec.at(step-1)<=lxborder)
	  {
        stepZero = step;
	    double dispIntegral = 0.0;
	    double force = 0.0;
        //long seed = noiseDiss.seed;  //fixed seed!!
	    switch(initCond)
	    {
	      case 0:
	      {
           xVec.at(step) = x0;
           velocVec.at(step) = v0;
	      }break;
          case 1:
	      {
           xVec.at(step) = x0;
           //const gsl_rng_type * T;
           //gsl_rng_env_setup();
           //gsl_rng * r;
           //T = gsl_rng_default;
           //r = gsl_rng_alloc(T);
           //gsl_rng_set (r, seed+step);
           velocVec.at(step) = gsl_ran_gaussian(res_rng, sqrt(k_b*temperature/mass)); //random initial velocity
           //cout << velocVec.at(step) << endl ;
	      }break;
	    }
	    	    
	    //noise
        //gsl_rng * r1;
	    //gsl_rng_set (r1, seed+step);
	    noiseDiss.noiseFunc(dt, noise_rng, colNoisecp);
        /*for(int i = 0; i < colNoisecp.size(); i++)
            {
            Noise1 << colNoise.at(i) << "          " << colNoisecp.at(i) << endl;
        }*/
	    if(step==nSteps-1)
	    { 
	      end = 1;
	      break;
	    }
	   step++;
	  }
	 }
	 if(end == 1)
	 {
	  break; 
	 }
      lastDispN = min(maxDispN, step-stepZero);
      //dispIntegral using trapezoid integral method
      dispIntegral = 0.5*(dissipationKernel.at(lastDispN-1)*velocVec.at(step-lastDispN)+
				  dissipationKernel.at(0)*velocVec.at(step-1));
	  for (int i = 1; i < lastDispN-1 ; i++)
	  {
	    dispIntegral += dissipationKernel.at(i)*velocVec.at(step-1-i);
	  }
	
      force = forceExtFunc(tVec.at(step-1), xVec.at(step-1));

      velocVec.at(step) = velocVec.at(step-1) + dt/mass*(-2.0*dispIntegral*dt + colNoisecp.at(step-1-stepZero) + force);
      xVec.at(step) = xVec.at(step-1) + dt*velocVec.at(step);
     //----------------------------------------------------- 
     /* if(step==5 | step== 20 | step== 30)
	 {
	  stepZero=step;
	  xVec.at(step)=x0;
	  velocVec.at(step)=v0;
	  const gsl_rng_type * T;
	    gsl_rng_env_setup();
	    T = gsl_rng_default;
	    long seed = time (NULL);
	    gsl_rng * r;
	    r = gsl_rng_alloc(T);
	    gsl_rng_set (r, seed+step);
	    noiseDiss.noiseFunc(so.dt, r, colNoisecp);
	    for(int i = 0; i < colNoisecp.size(); i++)
            {
            Noise1 << colNoise.at(i) << "          " << colNoisecp.at(i) << endl;
            }
	 }*/
    //-----------------------------
      /*dissInt.push_back(dispIntegral);
      xVector.push_back(xVec.at(step));
      vVector.push_back(velocVec.at(step));
      ForceVec.push_back(force);
      cout << step << " " << tVec.at(step) << endl;
      //cout << dissInt.at(step) << " " << xVector.at(step) << " " << vVector.at(step) << endl;
      /*dissInt.resize(step,dispIntegral);
      xVector.resize(step,xVec.at(step));
      vVector.resize(step,velocVec.at(step));*/
      
    }
  
    /*for(int i = 0; i< dissInt.size();i++)
    {
    diss << tVec.at(i) << "    " << dissInt.at(i) << endl;
    }
    for(int i = 0; i< xVector.size();i++)
    {
    x << tVec.at(i) << "    " << xVector.at(i) << endl;
    }
    for(int i = 0; i< vVector.size();i++)
    {
    v << tVec.at(i) << "    " << vVector.at(i) << endl;
    }
    for(int i = 0; i< ForceVec.size();i++)
    {
    f << tVec.at(i) << "    " << ForceVec.at(i) << endl;
    }*/
    //----------------------------------------
    return 1;
}

// Adams-Bashforth scheme for colored noise
// solves the Langevin-differential equation with colored Noise and dissipation kernel with the Adams-Bashforth-Method
// without implicit integral solving
// energy is not perfectly conserved
// if absolute position of particle becomes larger than xCutoff, all following positions will be set to the same value
// set xCutoff to negative value for turning that feature off
int adamsBColNoise(const vector<double>& tVec, const double& dt, const double& x0, const double& v0, const double& mass,const double& potw,
		     const double& k_b, const double& temperature,const vector<double>& colNoise, const vector<double>& dissipationKernel,
		     vector<double>& xVec, vector<double>& velocVec,function<double(const double&, const double&)> forceExtFunc, 
		     const double& xCutoff,const bool& KramersRate, const double& rxborder, const double lxborder,gsl_rng* res_rng,gsl_rng* noise_rng,
		     const int& initCond,const NoiseDiss& noiseDiss)
{
    int nSteps = xVec.size();
    int maxDispN = dissipationKernel.size();
    int lastDispN = 0;
    double dispIntegral1 = 0.0;
    double dispIntegral2 = 0.0;
    double dispIntegral3 = 0.0;
    double force1 = 0.0;
    double force2 = 0.0;
    double force3 = 0.0;
    int stepZero = 0;
    int end = 0;
    vector <double> colNoisecp;
    colNoisecp.resize(colNoise.size(),0.0);
    colNoisecp = colNoise;
    
    /*fstream Noise,Noise1;
    Noise.open("Noise.dat", ios::out);
    Noise1.open("Noise1.dat", ios::out);
    for(int i = 0; i < colNoisecp.size(); i++)
    {
    Noise << colNoise.at(i) << "          " << colNoisecp.at(i) << endl;
    }*/
     
    // --steps for t = dt, Euler-Cromer Method--
    xVec.at(0) = x0;
    velocVec.at(0) = v0;
    //cout << x0 <<" " <<  v0 << endl;

    dispIntegral1 = dissipationKernel.at(0)*velocVec.at(0);
    force1 = forceExtFunc(tVec.at(0), xVec.at(0));

    velocVec.at(1) = velocVec.at(0) + dt/mass*(-2.0*dispIntegral1*dt + colNoisecp.at(0) + force1);
    xVec.at(1) = xVec.at(0) + dt*velocVec.at(1);

    dispIntegral2 = dispIntegral1;
    force2 = force1;

    // --steps for t = 2dt, Adams-Bashforth n=2--
    dispIntegral1 = 0.5*(dissipationKernel.at(1)*velocVec.at(0)
			  + dissipationKernel.at(0)*velocVec.at(1));
    force1 = forceExtFunc(tVec.at(1), xVec.at(1));

    xVec.at(2) = xVec.at(1) + dt*(1.5*velocVec.at(1)-0.5*velocVec.at(0));
    velocVec.at(2) = velocVec.at(1) + dt/mass*(1.5*(-2.0*dispIntegral1*dt + colNoisecp.at(1) + force1)
						 -0.5*(-2.0*dispIntegral2*dt + colNoisecp.at(0) + force2));
    dispIntegral3 = dispIntegral2;
    dispIntegral2 = dispIntegral1;

    force3 = force2;
    force2 = force1;

    // --steps for t > 2dt, Adams-Bashforth n=3--
    for (int step = 3; step < nSteps; step++)
    {
	if ((abs(xVec.at(step-1)) > xCutoff) && (xCutoff > 0)){
	  xVec.at(step) = xVec.at(step-1);
	}
	else{
	  if(KramersRate == true)
	  {
	    while(xVec.at(step-1)>=rxborder || xVec.at(step-1)<=lxborder)
	    {
	     stepZero = step;
         switch(initCond)
         {
           case 0:
           {
            xVec.at(step) = x0;
            velocVec.at(step) = v0;
           }break;

           case 1:
           {
            xVec.at(step) = x0;
            velocVec.at(step) = gsl_ran_gaussian(res_rng, sqrt(k_b*temperature/mass));//random initial velocity
            //cout << velocVec.at(step) << endl ;
           }break;
					case 3:
           {
            xVec.at(step) = gsl_ran_gaussian(res_rng, sqrt(k_b*temperature/(mass*potw*potw)));//random initial velocity;
            velocVec.at(step) = gsl_ran_gaussian(res_rng, sqrt(k_b*temperature/mass));//random initial velocity;
           }break;


         }

         //noise
         noiseDiss.noiseFunc(dt, noise_rng, colNoisecp);
	    /*for(int i = 0; i < colNoisecp.size(); i++)
            {
            Noise1 << colNoise.at(i) << "          " << colNoisecp.at(i) << endl;
	    }*/
		    
	    if(step==nSteps-1)
	    {
	      end = 1;
	      break;
	    }
	    
	    dispIntegral1 = 0.0;
	    dispIntegral2 = 0.0;
	    dispIntegral3 = 0.0;
	    force1 = 0.0;
	    force2 = 0.0;
	    force3 = 0.0;
	    
	    dispIntegral1 = dissipationKernel.at(step-stepZero)*velocVec.at(step);
	    force1 = forceExtFunc(tVec.at(step), xVec.at(step));
	    
	    velocVec.at(step+1) = velocVec.at(step) + dt/mass*(-2.0*dispIntegral1*dt + colNoisecp.at(step-stepZero) + force1);
	    xVec.at(step+1) = xVec.at(step) + dt*velocVec.at(step+1);
	    
	    if(step==nSteps-2)
	    {
	     end = 1;
	     break; 
	    }
	    
	    dispIntegral2 = dispIntegral1;
	    force2 = force1;

	    // --steps for t = 2dt, Adams-Bashforth n=2--
	    dispIntegral1 = 0.5*(dissipationKernel.at(step+1-stepZero)*velocVec.at(step)
			  + dissipationKernel.at(step-stepZero)*velocVec.at(step+1));
	    force1 = forceExtFunc(tVec.at(step+1), xVec.at(step+1));

	    xVec.at(step+2) = xVec.at(step+1) + dt*(1.5*velocVec.at(step+1)-0.5*velocVec.at(step));
	    velocVec.at(step+2) = velocVec.at(step+1) + dt/mass*(1.5*(-2.0*dispIntegral1*dt + colNoisecp.at(step+1-stepZero) + force1)
						 -0.5*(-2.0*dispIntegral2*dt + colNoisecp.at(step-stepZero) + force2));
	    dispIntegral3 = dispIntegral2;
	    dispIntegral2 = dispIntegral1;

	    force3 = force2;
	    force2 = force1;
	    
	    if(step==nSteps-3)
	    {
	     end = 1;
	     break; 
	    }
	     step += 3;
	    }
	  }
	  
	  if(end == 1)
	  {
	   break; 
	  }

	  lastDispN = min(maxDispN, step-stepZero);
	  //dispIntegral using trapezoid integral method ("1/2-weighting" of boundary values and "1-weighting" inside the integration domain)
	  dispIntegral1 = 0.5*(dissipationKernel.at(lastDispN-1)*velocVec.at(step-lastDispN)+
				  dissipationKernel.at(0)*velocVec.at(step-1));
	  for (int i = 1; i < lastDispN-1 ; i++)
	  {
	    dispIntegral1 += dissipationKernel.at(i)*velocVec.at(step-1-i);
	  }

	  force1 = forceExtFunc(tVec.at(step-1), xVec.at(step-1));

	  xVec.at(step) = xVec.at(step-1) + dt*(23.0/12.0*velocVec.at(step-1)
							      -16.0/12.0*velocVec.at(step-2)
							      +5.0/12.0*velocVec.at(step-3));
	  velocVec.at(step) = velocVec.at(step-1)
				      + dt/mass*(23.0/12.0*(-2.0*dispIntegral1*dt +colNoisecp.at(step-1-stepZero) + force1)
						  -16.0/12.0*(-2.0*dispIntegral2*dt +colNoisecp.at(step-2-stepZero) + force2)
						  +5.0/12.0*(-2.0*dispIntegral3*dt +colNoisecp.at(step-3-stepZero) + force3));
	  dispIntegral3 = dispIntegral2;
	  dispIntegral2 = dispIntegral1;
	  force3 = force2;
	  force2 = force1;
	}
    }
    return 1;
}

SdeSolver::SdeSolver(const SimulationOptions& so, const NoiseDiss& noiseDiss){
	  switch(so.sdeSolverNr)
	  {
		case 0: this->sdeSolverFunc = bind(adamsBWhiteNoise, _1, so.dt, _2, _3, _4,so.potw,so.k_b,so.temperature, _5, so.gamma, so.D, _6, _7, _8,so.KramersRate,so.rxborder,so.lxborder,so.res_rng,noiseDiss.randGen.at(omp_get_thread_num()),so.initCondNr,noiseDiss); break;
        case 1: this->sdeSolverFunc = bind(eulerColNoise, _1, so.dt, _2, _3, _4, so.k_b,so.temperature,_5, noiseDiss.dispKernel, _6, _7, _8,so.KramersRate,so.rxborder,so.lxborder,so.res_rng,noiseDiss.randGen.at(omp_get_thread_num()),so.initCondNr,noiseDiss); break;
        case 2: this->sdeSolverFunc = bind(eulerCromerColNoise, _1, so.dt, _2, _3, _4, so.k_b,so.temperature, _5, noiseDiss.dispKernel, _6, _7, _8,so.KramersRate,so.rxborder,so.lxborder,so.res_rng,noiseDiss.randGen.at(omp_get_thread_num()),so.initCondNr,noiseDiss); break;
        case 3: this->sdeSolverFunc = bind(adamsBColNoise, _1, so.dt, _2, _3, _4,so.potw, so.k_b,so.temperature, _5, noiseDiss.dispKernel, _6, _7, _8, so.xCutoff,so.KramersRate,so.rxborder,so.lxborder,so.res_rng,noiseDiss.randGen.at(omp_get_thread_num()),so.initCondNr,noiseDiss); break;
		case 4: this->sdeSolverFunc = bind(eulerWhiteNoise, _1, so.dt, _2, _3, _4,so.k_b,so.temperature, _5, so.gamma, so.D, _6, _7, _8,so.KramersRate,so.rxborder,so.lxborder,so.res_rng,noiseDiss.randGen.at(omp_get_thread_num()),so.initCondNr,noiseDiss); break;
		default: printf("sdeSolverNr unknown\n");
		//todo ErrorMessage
		break;
	  }
}
