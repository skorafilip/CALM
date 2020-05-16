#ifndef _CALM_H_
#define _CALM_H_

#include <iostream>
#include <vector>
#include <list>
#include <string>
#include "TRandom2.h"
#include "TGenPhaseSpace.h"
#include "ParticleType.h"
#include "ParticleDB.h"
#include "Particle.h"
#include "THGlobal.h"
#include "TF1.h"
//--------reggae
#include "reggae.h"
#include "specrel.h"

using namespace std;

enum eEventType {GLOBAL, MINIJETS_GLOBAL, MINIJETS_LOCAL,GLOBAL_REGGAE, MINIJETS_GLOBAL_REGGAE, MINIJETS_LOCAL_REGGAE};

class CALM {
	public:
		CALM();
		~CALM();
		int   GenerateParticles(ParticleDB* aPartDB,int aMultBinMin, int aMultBinMax, double aEnergy, std::list<Particle>* aParticles, eEventType aEventType = GLOBAL);
		void   SetMultiplicities(ParticleDB* aDB,int aMultBinMin, int aMultBinMax);
		void   Randomize();
	private:
		TRandom2* mRandom;
      // values constant (what can be generated etc)
		int		 mNpart; // number of kinds of particles (pions, kaons, nuclides etc.)
		double* mNmean; // mean values for the particle yields dN/dy of each particle count -> taken from ALICE
		double mRapidityInterval; //interval of rapidity
		double* mXYZ; // to generate distance from primary vertex 
		int* mNpartkinds; // number of particles for each kind
		string** mNames; // names of particles to be generated
      TF1 *singleEnergyDistr; // energy distribution of single particle -> taken from Pythia
      TF1 *pionsMultDistr;
      TF1 *kaonsMultDistr;
      TF1 *nucleonsMultDistr;
      TF1 *lambdasMultDistr;
		double Etot; //Total momentum to be distributed among particles
      // values for this event
      vector<string> mParticlesThisEvent;

	private:
      int* GetMultiplicitiesOfPartciles(int pythiaMult,int aMultBinMin, int aMultBinMax, int& Nsum);
	  void CheckConservAtionLaws(int* Nrand, vector<vector<int>>& Npart, ParticleDB* aPartDB);
	  //double** GetXYZ(int Nsum);
	  vector<vector<double>> GetXYZ(int Nsum);
	  double GetTotalEnergy(int Nsum);
	  double* GetMasses(int Nsum, ParticleDB* aPartDB);
	  bool TrySetEventDecay(int Nsum, double* masses, TGenPhaseSpace& event, double& TotEnergy);
	  bool FilterUnlikelyEvents(TGenPhaseSpace& event, double& weight);
	  void SaveAllParticles_GLOBAL(int Nsum, double weight, vector<vector<double>> XYZrand, TGenPhaseSpace event, ParticleDB* aPartDB, list<Particle>* aParticles);
	  void SaveAllParticles_GLOBAL_REGGAE(int Nsum, vector4 *avec, vector<vector<double>> XYZrand, ParticleDB* aPartDB, list<Particle>* aParticles);
	  bool DOREGGAE(int Nsum, double* masses, vector4 en, vector4* avec);

	  void SeparateJets(int Nsum, vector<double>* masses, vector<string>* names, ParticleDB* aPartDB);
	  bool TrySetEventDecay_MINIJETS(int Nsum, double* masses0, double* masses1, TGenPhaseSpace& event0, TGenPhaseSpace& event1, double& TotEnergy, double* divideEn, vector<double>* masses);
	  bool FilterUnlikelyEvents_MINIJETS(TGenPhaseSpace& event0, TGenPhaseSpace& event1, double& weight0, double& weight1);
	  void SaveAllParticles_MINIJETS(vector<double>* masses, vector<string>* names, double weight0, double weight1, double TotEnergy, double* divideEn, vector<vector<double>> XYZrand, TGenPhaseSpace event0, TGenPhaseSpace event1, ParticleDB* aPartDB, list<Particle>* aParticles, eEventType aEventType);
	 
	  bool SeparateJets_LOCAL(int Nsum, vector<double>* masses, vector<string>* names, ParticleDB* aPartDB);


};

#endif
