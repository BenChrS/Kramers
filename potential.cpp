
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

#include "simulationOptions.h"
#include "potential.h"

using namespace std;
using namespace std::tr1;
using namespace std::tr1::placeholders;

// ----potential-functions----

// no force
double forceExt0(const double& t, const double& x)
{
  return 0.0;
}

//no potential
double potential0(const double& t, const double& x)
{
  return 0.0;
}

//force of harmonic oscillator
double forceExt1(const double& t, const double& x, const double& k)
{
  return -k*x;
}

//potential of harmonic oscillator
double potential1(const double& t, const double& x, const double& k)
{
  return 1.0/2.0*k*pow(x,2);
}

//force of double well potential, potDepth is the well height and potWidth is the distance between minimum and maximum of the potential
double forceExt2(const double& t, const double& x, const double& potDepth, const double& potWidth)
{
  return 4.0*potDepth/pow(potWidth,2)*x - 4.0*potDepth/pow(potWidth,4)*pow(x,3);
}

//potential of double well potential, potDepth is the well height and potWidth is the distance between minimum and maximum of the potential
double potential2(const double& t, const double& x, const double& potDepth, const double& potWidth)
{
  return -2.0*potDepth/pow(potWidth,2)*pow(x,2) + potDepth/pow(potWidth,4)*pow(x,4) + potDepth;
}

//force of assymetric mexican hat
double forceExt3(const double& t, const double& x, const double& a, const double& b, const double& c)
{
  return -4.0*a*pow(x,3) - 3.0*b*pow(x,2) - 2.0*c*x;
}

//force of assymetric mexican hat
double potential3(const double& t, const double& x, const double& a, const double& b, const double& c)
{
  return a*pow(x,4) + b*pow(x,3) + c*pow(x,2);
}

//Mexican hat ; my y-Achsenaschnittspunkt, a: Abstand Maximum/Minimum, ny: Höhe der rechten Potenzialmulde
double forceExt4(const double& t, const double& x,const double& pota, const double& potMy, const double& potNy)
{
  double A = -(-2*potMy+potNy)/(2*pow(pota,4.0));
  double B = -(potNy)/(4*pow(pota,3.0));
  double C = -(2*potMy-potNy)/(pow(pota,2.0));
  double D = (3*potNy)/(4*pota);
  return -4*A*pow(x,3.0)-3*B*pow(x,2.0)-2*C*x-D;
}

//x**3-potential, potDepth is the well height and potWidth is the distance between minimum and maximum of the potential
double potential4(const double& t, const double& x,const double& pota, const double& potMy, const double& potNy)
{
  double A = -(-2*potMy+potNy)/(2*pow(pota,4.0));
  double B = -(potNy)/(4*pow(pota,3.0));
  double C = -(2*potMy-potNy)/(pow(pota,2.0));
  double D = (3*potNy)/(4*pota);
  double E = potMy;
  return A*pow(x,4.0)+B*pow(x,3.0)+C*pow(x,2.0)+D*x+E;
}



 //kombinierte Potenzialschwelle aus 2 harmonischen Funktionen;
 double forceExt5(const double& t, const double& x,const double& pota, const double& potMy, const double& potNy)
 {
   double A = -(-2*potMy+potNy)/(2*pow(pota,4.0));
   double C = -(2*potMy-potNy)/(pow(pota,2.0));
   double waSqr = 12*A*pow(pota,2.0)+2*C;
   
   if(x<0)
   {
     return -waSqr*(x+pota);
   }
   else if(x==0)
   {
    return 0; 
   }
   else
   {
    return -waSqr*(x-pota); 
   }
 }
 
 //kombinierte Potenzialschwelle aus 2 harmonischen Funktionen;
 double potential5(const double& t, const double& x,const double& pota, const double& potMy, const double& potNy)
 {
   double A = -(-2*potMy+potNy)/(2*pow(pota,4.0));
   double C = -(2*potMy-potNy)/(pow(pota,2.0));
   double waSqr = 12*A*pow(pota,2.0)+2*C;
   
   if(x<=0)
   {
    return 1.0/2*waSqr*pow(x+pota,2.0); 
   }
   else
   {
    return 1.0/2*waSqr*pow(x-pota,2.0);
   }
 }

//kombinierte Potenzialschwelle aus 2 harmonischen Funktionen mit modifizierter Frequenz; Test 
double forceExt6(const double& t, const double& x,const double& pota)
{
  
  if(x<0)
  {
    return -50.0*(x+pota);
  }
  else if(x==0)
  {
   return 0; 
  }
  else
  {
   return -50.0*(x-pota); 
  }
}

//kombinierte Potenzialschwelle aus 2 harmonischen Funktionen mit modifizierter Frequenz, derzeitige Frequenz nur für spezielle 
//Parameterwahl;
double potential6(const double& t, const double& x,const double& pota)
{
 
  if(x<=0)
  {
   return 25.0*pow(x+pota,2.0); 
  }
  else
  {
   return 25.0*pow(x-pota,2.0);
  }
}

//two smoothly joint parabolas (Paper: thermal decay of a metastable state)
double forceExt7(const double& t, const double& x,const double xc, const double xb, const double Ub,const double mass)
{
  
double xm, C0, wc;
xm = (xb+xc)/2.0;
wc= sqrt(4.0*Ub/(mass*pow((xb-xc),2.0)));
C0 = mass*pow(wc,2.0);

if(x<xm)
{
  return -C0*(x-xc);
}
else
{
  return +C0*(x-xb);
}

}

//two smoothly joint parabolas (Paper: thermal decay of a metastable state)
//for arbitrary mass
double potential7(const double& t, const double& x,const double xc, const double xb, const double Ub,const double mass)
{
double xm, C0, wc;
xm = (xb+xc)/2.0;
wc= sqrt(4.0*Ub/(mass*pow((xb-xc),2.0)));
C0 = mass*pow(wc,2.0);

if(x<xm)
{
  return C0*pow((x-xc),2.0)/2.0;
}
else
{
  return Ub-C0*pow((x-xb),2.0)/2.0;
}
}

// //two smoothly joint parabolas (Paper: thermal decay of a metastable state)
// //only if mass=1*temperature
// double potential7(const double& t, const double& x,const double xc, const double xb, const double Ub)
// {
// double xm, C0;
// xm = (xb+xc)/2.0;
// C0 = 4*Ub/pow((xb-xc),2.0);
// 
// if(x<xm)
// {
//   return C0*pow((x-xc),2.0)/2.0;
// }
// else
// {
//   return Ub-C0*pow((x-xb),2.0)/2.0;
// }
// }



Potential::Potential(SimulationOptions& so){
	// select and define local potential, 0: no force, 1: harmonic, 2: double well, 3: asymmetric mexican hat
	  switch(so.potNr)
	  {
		case 0: {forceFunc = forceExt0; potentialFunc = potential0;} break;
		case 1: {forceFunc = bind(forceExt1, _1, _2, so.potK);
			  potentialFunc = bind(potential1, _1, _2, so.potK);} break;
		case 2: {forceFunc = bind(forceExt2, _1, _2, so.potDepth, so.potWidth);
			  potentialFunc = bind(potential2, _1, _2, so.potDepth, so.potWidth);} break;
		case 3: {forceFunc = bind(forceExt3, _1, _2, so.potA, so.potB, so.potC);
			  potentialFunc = bind(potential3, _1, _2, so.potA, so.potB, so.potC);} break;
		case 4: {forceFunc = bind(forceExt4, _1, _2, so.pota,so.potMy,so.potNy);
			  potentialFunc = bind(potential4, _1, _2,so.pota,so.potMy,so.potNy);} break;
		case 5: {forceFunc = bind(forceExt5, _1, _2, so.pota,so.potMy,so.potNy);
			  potentialFunc = bind(potential5, _1, _2,so.pota,so.potMy,so.potNy);} break;
	        case 6: {forceFunc = bind(forceExt6, _1, _2, so.pota);
			  potentialFunc = bind(potential6, _1, _2,so.pota);} break;
		case 7: {forceFunc = bind(forceExt7, _1, _2, so.xc,so.xb,so.Ub,so.mass);
			  potentialFunc = bind(potential7, _1, _2,so.xc,so.xb,so.Ub,so.mass);} break;
		default: {printf("potNr unknown\n");
				//todo error message
				} break;
		}
}
