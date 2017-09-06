#include <stdio.h>
#include <iostream>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>

#include "simulationOptions.h"
#include "vectorCalculations.h"
#include "results.h"

using namespace std;

Results::Results(const SimulationOptions &so){
	setInitValues(so);
	this->tVec.resize(so.nSteps, 0.0);
	linspace(so.t0, so.dt, this->tVec);
	this->allX.resize(so.np, vector<double>(so.nSteps, 0.0));
	this->allV.resize(so.np, vector<double>(so.nSteps, 0.0));

}

void Results::setInitValues(const SimulationOptions &so){
	this->x0Vec.resize(so.np, 0.0);
	this->v0Vec.resize(so.np, 0.0);

	switch(so.initCondNr)
		  {
			// fixed space and velocity values
			case 0:
			{
				fill(this->x0Vec.begin(), this->x0Vec.end(), so.x0);
				fill(this->v0Vec.begin(), this->v0Vec.end(), so.v0);
			} break;
			// fixed space value and thermal velocity
			case 1:
			{
				fill(this->x0Vec.begin(), this->x0Vec.end(), so.x0);
                this->v0Vec = fillThermalized(so);
			} break;
			// Gaussian (thermal) distribution in space and fixed velocity
			case 2:
			{
                this->x0Vec = fillThermalized(so);
				fill(this->v0Vec.begin(), this->v0Vec.end(), so.v0);
			} break;
			// Gaussian (thermal) distribution in space and velocity
			case 3:
			{
                this->x0Vec = fillThermalized(so);
                this->v0Vec = fillThermalized(so);
			} break;
			default: {
				printf("initCondNr unknown\n");
				//todo: add error message;
			} break;

		  }
}

// the thermalised velocity of the "Brownian particle" is distributed according to a Gaussian
vector<double> Results::fillThermalized(const SimulationOptions &so)
{
  vector<double> output(so.np);
  //const gsl_rng_type * T;
  //gsl_rng_env_setup();
  //T = gsl_rng_default;
  //long seed = time (NULL);
  //gsl_rng * r;
  //so.res_rng = gsl_rng_alloc (T);
  cout << time(NULL) << endl;       //bei wenigen Teilchen ist time(NULL) an dieser Stelle identisch mit time(NULL) in noiseDiss
  gsl_rng_set (so.res_rng, so.seed);
  for (int i = 0; i < so.np; i++)
  {
    output.at(i) = gsl_ran_gaussian(so.res_rng, sqrt(so.k_b*so.temperature/so.mass));
  }
  //gsl_rng_free(so.res_rng);  //würde bei nochmaligem Aufruf von so.res zu einem Fehler führen, weil dieses nicht mehr existent
  return output;
}

KappaSimulation::KappaSimulation(const SimulationOptions &so){
	this->transition.resize(so.nSteps,0.0);
	this->kappa.resize(so.nSteps,0.0);
	this->indexNr=0;
	this->xAvVec.resize(so.avNum,vector<double>(so.nSteps, 0.0));
	this->xSquaredAvVec.resize(so.avNum,vector<double>(so.nSteps, 0.0));
	this->vAvVec.resize(so.avNum,vector<double>(so.nSteps, 0.0));
	this->vSquaredAvVec.resize(so.avNum,vector<double>(so.nSteps, 0.0));
	this->transitionAvVec.resize(so.avNum,vector<double>(so.nSteps, 0.0));
	this->fluxAvVec.resize(so.avNum,vector<double>(so.nSteps, 0.0));
	this->fluxPositiveAvVec.resize(so.avNum,vector<double>(so.nSteps, 0.0));
	this->fluxNegativeAvVec.resize(so.avNum,vector<double>(so.nSteps, 0.0));
	this->fluxPaperAvVec.resize(so.avNum,vector<double>(so.nSteps, 0.0));
	this->kinEnergyAvVec.resize(so.avNum,vector<double>(so.nSteps, 0.0));
	this->pDistAvVec.resize(so.avNum,vector< vector<double> >(ceil(double(so.nSteps)/double(so.pdaStride)),vector<double>(so.pDistNBins,0.0)));
	this->vDistAvVec.resize(so.avNum,vector< vector<double> >(ceil(double(so.nSteps)/double(so.vdaStride)),vector<double>(so.vDistNBins,0.0)));
	/* Vektoren p/vDistAvVec haben avNum Einträge, von denen wiederum jeder Eintrag nSteps/Stride Einträge hat (entspricht Zahl der betrachteten Zeitpunkte)
	 jeder dieser nSteps/stride Einträge hat wiederum nBins Einträge  */ 
	this->indexAvNum=0;
	this->CorrAvVec.resize(so.avNum,vector<double>(so.nSteps, 0.0));
	this->noiseAvVec.resize(so.avNum,vector<double>(so.nSteps, 0.0));
	this->noiseSqrAvVec.resize(so.avNum,vector<double>(so.nSteps, 0.0));
	this->corrXAv.resize(so.avNum,vector<double>((int) round((so.t1-so.t0-so.timeSettled)/so.dt), 0.0));
	this->corrVAv.resize(so.avNum,vector<double>((int) round((so.t1-so.t0-so.timeSettled)/so.dt), 0.0));



}
