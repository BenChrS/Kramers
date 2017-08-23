#pragma once

#include <stdio.h>
#include <iostream>
#include <stdio.h>

#include <tr1/functional>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>

#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>

#include "simulationOptions.h"

using namespace std;
using namespace std::tr1;

// ----potential-functions----
// no force
double forceExt0(const double& t, const double& x);

//no potential
double potential0(const double& t, const double& x);

//force of harmonic oscillator
double forceExt1(const double& t, const double& x, const double& k);


//potential of harmonic oscillator
double potential1(const double& t, const double& x, const double& k);


//force of double well potential, potDepth is the well height and potWidth is the distance between minimum and maximum of the potential
double forceExt2(const double& t, const double& x, const double& potDepth, const double& potWidth);

//potential of double well potential, potDepth is the well height and potWidth is the distance between minimum and maximum of the potential
double potential2(const double& t, const double& x, const double& potDepth, const double& potWidth);


//force of asymmetric Mexican hat
double forceExt3(const double& t, const double& x, const double& a, const double& b, const double& c);


//potential of asymmetric Mexican hat
double potential3(const double& t, const double& x, const double& a, const double& b, const double& c);

//force of asymmetric Mexican hat, pota: distance between 0 an minima, potMy: f(-pota), potNy: y-distance between minimum and y-value at x=0
double forceExt4(const double& t, const double& x,const double& pota, const double& potMy, const double& potNy);

//potential of asymmetric Mexican hat, pota: distance between 0 an minima, potMy: f(-pota), potNy: y-distance between minimum and y-value at x=0
//not working
double potential4(const double& t, const double& x,const double& pota, const double& potMy, const double& potNy);

//kombinierte Potenzialschwelle aus 2 harmonischen Funktionen;
double forceExt5(const double& t, const double& x,const double& pota, const double& potMy, const double& potNy);

//kombinierte Potenzialschwelle aus 2 harmonischen Funktionen;
double potential5(const double& t, const double& x,const double& pota, const double& potMy, const double& potNy);

//kombinierte Potenzialschwelle aus 2 harmonischen Funktionen mit modifizierter Frequenz;
double forceExt6(const double& t, const double& x,const double& pota);

//kombinierte Potenzialschwelle aus 2 harmonischen Funktionen mit modifizierter Frequenz;
double potential6(const double& t, const double& x,const double& pota);

//two smoothly joint parabolas (Paper: thermal decay of a metastable state)
double forceExt7(const double& t, const double& x,const double xc, const double xb, const double Ub,const double mass);

//two smoothly joint parabolas (Paper: thermal decay of a metastable state)
double potential7(const double& t, const double& x,const double xc, const double xb, const double Ub, const double mass);

class Potential{
public:
	function<double(const double&, const double&)> forceFunc;
	function<double(const double&, const double&)> potentialFunc;

	Potential(SimulationOptions&);
};
