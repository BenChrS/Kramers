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

#include "vectorCalculations.h"

using namespace std;
using namespace std::tr1;
using namespace std::tr1::placeholders;

// returns vector with length nSteps with equally spaced values
// starting from t0 with distance dt to the next value
void linspace(const double &t0, const double &dt, vector<double> &vec)
{
  double t = t0;
  for (int i = 0; i < vec.size(); i++)
  {
    vec.at(i) = t;
    t += dt;
  }
}

// generate a symmetric or asymmetric vector of the input data (double length)
// in the first half of the returned vector the original content of data is stored, in the second half the content is stored in reversed order
// if negative is set to true, the reversed half gets a minus sign
void mirrorData(vector<double> &vec, const vector<double>& data, const bool& negative, const bool& allZero)
{
	vec.resize(2*data.size(), 0.0);
  int sign = 1;
  int allZeroSign = 1;
  if (negative)
  {
    sign = -1;
  }
  if (allZero)
  {
    allZeroSign = 0;
  }
  copy(data.begin(), data.end(), vec.begin());
  copy(data.begin(), data.end(), vec.begin() +data.size());
  reverse(vec.begin() +data.size(), vec.end());
  transform(vec.begin() +data.size(), vec.end(), vec.begin() +data.size(), bind1st(multiplies<double>(),sign*allZeroSign));
}

// returns the row of a two dimensional vector
void getRow(vector<double> &vec, const vector< vector<double> >& data, const int& row)
{
  vec.resize(data.size());
  for (int i = 0; i<data.size(); i++)
  {
    vec.at(i) = data.at(i).at(row);
  }
}

//calculates the square of a vector
void vectorSquared(vector<double> &vec, const vector<double>& data)
{
  vec.resize(data.size());
  transform(data.begin(), data.end(), data.begin(), vec.begin(), multiplies<double>());
}

//adds 2 vectors
void add2Vec(vector<double> &vec, const vector<double>& vec1, const vector<double>& vec2)
{
  vec.resize(vec1.size(), 0.0);
  for (int i = 0; i<vec1.size(); i++){
		  vec.at(i) = vec1.at(i) + vec2.at(i);
  }
}

// builds average over vectors of vectors, all inner vectors must have the same size
void buildAverage(vector<double> &vec, const vector< vector<double> >& data)
{
  vec.resize(data.at(0).size(), 0.0); // vec.resize dient der Anpassung der Dimensionen der Vektoren
  for (int i = 0; i < data.size(); i++){ // ist die Dimension von data größer, so werden die fehlenden Einträge ergänzt und mit 0.0 initialisiert, die übrigen Werte bleiben erhalten
	  add2Vec(vec, vec, data.at(i)); //does this work???
  }
  transform(vec.begin(), vec.end(), vec.begin(), bind1st(multiplies<double>(), 1.0/data.size()));
}

// builds average of one vector
double buildAverage2(vector<double>::iterator first, vector<double>::iterator last)
{
  double output;
  output = accumulate(first, last, 0.0, plus<double>());
  return (output/distance(first,last));
}

//builds average over vectors of vectors of vectors
void buildAverage3(vector<vector<double> > &vec, const vector< vector< vector<double> > >& data)
{ /* data entspricht z.B. ksim.vDistAvVec , vec entspricht einem Eintrag des "äußeren" Vektors mit avNum Einträgen */
  vec.resize(data.at(0).size(), vector<double>(data.at(0).at(0).size(),0.0)); /*Vektor mit nSteps/stride Einträgen, die 
  jeweils nBins Einträge haben*/
  for (int i = 0; i < data.at(0).size(); i++) // nSteps/stride-Ebene --> i entspricht Zeitpunkt
  {
    for(int j=0; j < data.at(0).at(0).size(); j++) //nBins-Ebene --> j entspricht Bin
    {
	  for(int k=0; k< data.size();k++) //avNum-Ebene --> k entspricht Durchgang der Simulation
	  {
	   vec.at(i).at(j) = vec.at(i).at(j) + data.at(k).at(i).at(j); 
	  }
    }
  }
  for (int i=0; i < vec.size();i++)
  {
  transform(vec.at(i).begin(), vec.at(i).end(), vec.at(i).begin(), bind1st(multiplies<double>(), 1.0/data.size()));
  /* jede Verteilung zum Zeitpunkt i wird genau data.size()=k mal erzeugt, Mittelung erfolgt, indem durch diese Anzahl 
   geteilt wird*/ 
  }
    
  }

// builds squared average over vector of vectors, all inner vectors must have the same size
void buildAverageSquared(vector<double> &vec, const vector< vector<double> >& data)
{
  vector< vector<double> > squaredTemp(data.size(), vector<double>(data.at(0).size(), 0.0));
  for (int i = 0; i < data.size(); i++)
  {
     vectorSquared(squaredTemp.at(i), data.at(i));
  }
  buildAverage(vec, squaredTemp);
}

// builds average over squared vector, all inner vectors must have the same size
double buildAverageSquared2(const vector<double>& data)
{
   vector<double> squaredTemp(data.size(), 0.0);
   vectorSquared(squaredTemp, data);
   return buildAverage2(squaredTemp.begin(), squaredTemp.end());
}

// returns the standard deviation of the given data
double stdDeviation( vector<double>::iterator first, vector<double>::iterator last)
{
  double mean = buildAverage2(first, last);
  vector<double> diff(distance(first,last));
  transform(first, last, diff.begin(),bind2nd(minus<double>(), mean));
  double sq_sum = inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
  double stdev = sqrt(sq_sum / distance(first,last));
  return stdev;
}
