#pragma once

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

using namespace std;
using namespace std::tr1;

// returns vector with length nSteps with equally spaced values
// starting from t0 with distance dt to the next value
void linspace(const double &t0, const double &dt, vector<double> &vec);

// in the first half of the returned vector the original content of data is stored, in the second half the content is stored in reversed order
// if negative is set to true, the reversed half gets a minus sign
void mirrorData(vector<double> &vec, const vector<double>& data, const bool& negative = false, const bool& allZero = false);

// returns the row of a two dimensional vector
void getRow(vector<double> &vec, const vector< vector<double> >& data, const int& row);

//calculates the square of a vector
void vectorSquared(vector<double> &vec, const vector<double>& data);

//adds 2 vectors
void add2Vec(vector<double>& vec, const vector<double>& vec1, const vector<double>& vec2);

// builds average over vectors of vectors, all inner vectors must have the same size
void buildAverage(vector<double> &vec, const vector< vector<double> >& data);

// builds average of one vector
double buildAverage2(vector<double>::iterator first, vector<double>::iterator last);

void buildAverage3(vector<vector<double> > &vec, const vector< vector< vector<double> > >& data);

// builds squared average over vector of vectors, all inner vectors must have the same size
void buildAverageSquared(vector<double> &vec, const vector< vector<double> >& data);

// builds average over squared vector, all inner vectors must have the same size
double buildAverageSquared2(const vector<double>& data);

// returns the standard deviation of the given data
double stdDeviation( vector<double>::iterator first, vector<double>::iterator last);
