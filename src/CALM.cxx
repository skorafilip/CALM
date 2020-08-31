#include "CALM.h"

//extern Configurator *sMainConfig;

void AddParticleSums(int &Qsum, int &Ssum, int &Bsum, string particleName, ParticleType *tParticleType)
{
   int tmpInt;
   //charge check
   if (particleName.find("plu") != std::string::npos)
      Qsum++;
   else if (particleName.find("min") != std::string::npos || particleName.find("plb") != std::string::npos)
      Qsum--;
   else if (particleName.find("zer") != std::string::npos || particleName.find("zrb") != std::string::npos)
      ;
   //barion number check
   tmpInt = tParticleType->GetNumberQ() - tParticleType->GetNumberAQ() + tParticleType->GetNumberS() - tParticleType->GetNumberAS();
   if (tmpInt == 3)
      Bsum++;
   else if (tmpInt == -3)
      Bsum--;
   //strangeness check
   tmpInt = tParticleType->GetNumberS() - tParticleType->GetNumberAS();
   if (tmpInt == 1)
      Ssum--; //  for quark s: S=-1
   else if (tmpInt == -1)
      Ssum++;
}


int *CALM::PythiaMult(int aMultBinMin, int aMultBinMax, int &Nsum)
{
   // number of particles generated (for each kind) - from Pytia distribution
   int *Nrand = new int[mNpart];
   int pythiaMult = eventConfig->pythiaMult;

   pionsMultDistr = new TF1("pionsMultDistr", eventConfig->pionsMultDistr.c_str(), eventConfig->pionsMultDistr_xMin, eventConfig->pionsMultDistr_xMax);
   kaonsMultDistr = new TF1("kaonsMultDistr", eventConfig->kaonsMultDistr.c_str(), eventConfig->kaonsMultDistr_xMin, eventConfig->kaonsMultDistr_xMax);
   nucleonsMultDistr = new TF1("nucleonsMultDistr", eventConfig->nucleonsMultDistr.c_str(), eventConfig->nucleonsMultDistr_xMin, eventConfig->nucleonsMultDistr_xMax);
   lambdasMultDistr = new TF1("lambdasMultDistr", eventConfig->lambdasMultDistr.c_str(), eventConfig->lambdasMultDistr_xMin, eventConfig->lambdasMultDistr_xMax);
   // generating the number of particles in each kind
   do
   {
      Nsum = 0;
      Nrand[0] = pionsMultDistr->GetRandom();
      Nrand[1] = kaonsMultDistr->GetRandom();
      Nrand[2] = nucleonsMultDistr->GetRandom();
      Nrand[3] = lambdasMultDistr->GetRandom();
      Nsum = Nrand[0] + Nrand[1] + Nrand[2] + Nrand[3];
   } while (Nsum < aMultBinMin || Nsum > aMultBinMax || (Nrand[1] + Nrand[3]) % 2 != 0 || (Nrand[2] + Nrand[3]) % 2 != 0);
   delete pionsMultDistr;
   delete kaonsMultDistr;
   delete nucleonsMultDistr;
   delete lambdasMultDistr;

   return Nrand;
}

int *CALM::AlicePoissonMult(int aMultBinMin, int aMultBinMax, int &Nsum)
{
   // number of particles generated (for each kind) - from Poisson distribution
   int *Nrand = new int[mNpart];
   int pythiaMult = eventConfig->pythiaMult;

   do
   {
      Nsum = 0;
      for (int i = 0; i < mNpart; ++i)
      {
         Nrand[i] = mRandom->Poisson(mNmean[i] * mRapidityInterval * mNpartkinds[i]);
         Nsum += Nrand[i];
      }
   } while (Nsum < aMultBinMin || Nsum > aMultBinMax || (Nrand[1] + Nrand[3]) % 2 != 0 || (Nrand[2] + Nrand[3]) % 2 != 0);

   return Nrand;
}

int **CALM::GetTypesForParticles(int *Nrand, ParticleDB *aPartDB)
{
   int Qsum, Bsum, Ssum;

   ParticleType *tParticleType;

   int **Npart = new int*[mNpart];
   for(int i = 0; i < mNpart; i++)
     Npart[i] = new int[Nrand[i]];

   do
   {
      Qsum = 0;
      Ssum = 0;
      Bsum = 0;
      // generating the number of specific particles within each kind
      // check of the charge, strangeness and baryon number
      for (int i = 0; i < mNpart; ++i)
      {
         for (int j = 0; j < Nrand[i]; ++j)
         {
            Npart[i][j] = (int)mRandom->Uniform(mNpartkinds[i]);
            tParticleType = aPartDB->GetParticleType(mNames[i][Npart[i][j]].c_str());

            AddParticleSums(Qsum, Bsum, Ssum, mNames[i][Npart[i][j]], tParticleType);
         }
      }
   } while (Qsum != 0 || Ssum != 0 || Bsum != 0);

   return Npart;
}

double ** CALM::GetVertexXYZ(int Nsum)
{
   double ** XYZrand = new double*[Nsum];
   for (int j = 0; j < Nsum; ++j)
   {
      XYZrand[j] = new double[3];
      for (int i = 0; i < 3; ++i)
      {
         XYZrand[j][i] = mRandom->Gaus(0, mXYZ[i]);
      }
   }
   return XYZrand;
}

void CALM::SetTotalEnergy(int Nsum, double aEnergy)
{
   //double Etot;
   singleEnergyDistr = new TF1("singleEnergyDistr", eventConfig->singleEnergyDistr.c_str(), eventConfig->singleEnergyDistr_xMin, eventConfig->singleEnergyDistr_xMax);
   do
   {
      Etot = 0.;
      for (int i = 0; i < Nsum; ++i)
         Etot += singleEnergyDistr->GetRandom(eventConfig->singleEnergyDistr_xMin, eventConfig->singleEnergyDistr_xMax);
   } while (Etot > aEnergy);
   delete singleEnergyDistr;
   //return Etot;
}

double *CALM::GetMasses(int Nsum, ParticleDB *aPartDB)
{
   double *masses = new double[Nsum];
   for (int i = 0; i < Nsum; i++)
   {
      masses[i] = aPartDB->GetParticleType(mParticlesThisEvent[i].c_str())->GetMass();
   }
   return masses;
}

void CALM::SeparateJets(int Nsum, vector<double> *masses, vector<string> *names, ParticleDB *aPartDB)
{
   do
   {
      if (masses[0].size() > 0 || masses[1].size() > 0)
      {
         masses[0].clear();
         masses[1].clear();
         names[0].clear();
         names[1].clear();
      }
      for (int i = 0; i < Nsum; ++i)
      {

         if (mRandom->Integer(2))
         {
            masses[1].push_back(aPartDB->GetParticleType(mParticlesThisEvent[i].c_str())->GetMass());
            names[1].push_back(mParticlesThisEvent[i].c_str());
         }
         else
         {
            masses[0].push_back(aPartDB->GetParticleType(mParticlesThisEvent[i].c_str())->GetMass());
            names[0].push_back(mParticlesThisEvent[i].c_str());
         }
      }
   } while (masses[0].size() < 4 || masses[1].size() < 4);
}

bool CALM::SeparateJets_LOCAL(int Nsum, vector<double> *masses, vector<string> *names, ParticleDB *aPartDB)
{
   ParticleType *tParticleType;
   bool isSuccess = false;
   int control = 0;
   int Qjet[2], Bjet[2], Sjet[2];
   //int tmpInt;

   do
   {
      if (masses[0].size() > 0 || masses[1].size() > 0)
      {
         masses[0].clear();
         masses[1].clear();
         names[0].clear();
         names[1].clear();
      }
      for (int it_clean = 0; it_clean < 2; it_clean++)
      {
         Qjet[it_clean] = 0;
         Sjet[it_clean] = 0;
         Bjet[it_clean] = 0;
      }
      for (int i = 0; i < Nsum; ++i)
      {
         if (mRandom->Integer(2))
         {
            tParticleType = aPartDB->GetParticleType(mParticlesThisEvent[i].c_str());
            masses[1].push_back(tParticleType->GetMass());
            names[1].push_back(mParticlesThisEvent[i].c_str());

            AddParticleSums(Qjet[0], Bjet[0], Sjet[0], mParticlesThisEvent[i], tParticleType);
         }
         else
         {
            tParticleType = aPartDB->GetParticleType(mParticlesThisEvent[i].c_str());
            masses[0].push_back(tParticleType->GetMass());
            names[0].push_back(mParticlesThisEvent[i].c_str());

            AddParticleSums(Qjet[1], Bjet[1], Sjet[1], mParticlesThisEvent[i], tParticleType);
         }
      }
      isSuccess = Qjet[0] == 0 && Sjet[0] == 0 && Bjet[0] == 0 && Qjet[1] == 0 && Sjet[1] == 0 && Bjet[1] == 0 && masses[0].size() >= 4 && masses[1].size() >= 4;
      control++;
   } while (!(isSuccess || control > 100)); //Qjet[0] != 0 || Sjet[0] != 0 || Bjet[0] != 0 || Qjet[1] != 0 || Sjet[1] != 0 || Bjet[1] != 0 || masses[0].size() < 4 || masses[1].size() < 4);
   return isSuccess;
}
bool CALM::TrySetEventDecay(int Nsum, double *masses, TGenPhaseSpace &event, double &TotEnergy)
{
   bool isSuccess = false;
   int control = 0;
   TLorentzVector en;
   do
   {
      // generate total energy
      TotEnergy = Etot;
      en.SetE(TotEnergy);

      isSuccess = event.SetDecay(en, Nsum, masses);
      control++;
   } while (!(isSuccess || control > 10));
   return isSuccess;
}

bool CALM::TrySetEventDecay_MINIJETS(int Nsum, vector<double> *masses, TGenPhaseSpace &event0, TGenPhaseSpace &event1, double &TotEnergy, double *divideEn)
{
   double masses0[masses[0].size()];
   double masses1[masses[1].size()];
   for (int j = 0; j < masses[0].size(); ++j)
   {
      masses0[j] = masses[0][j];
   }
   for (int j = 0; j < masses[1].size(); ++j)
   {
      masses1[j] = masses[1][j];
   }

   bool isSuccess = false;
   int control = 0;
   TLorentzVector en;
   do
   {
      TotEnergy = Etot;
      en.SetE(TotEnergy * (divideEn[0] / (2. * (divideEn[0] + divideEn[1]))));

      isSuccess = (event0.SetDecay(en, masses[0].size(), masses0) && event1.SetDecay(en, masses[1].size(), masses1));
      control++;
   } while (!(isSuccess || control > 10));
   return isSuccess;
}

bool CALM::FilterUnlikelyEvents(TGenPhaseSpace &event, double &weight)
{
   bool isSuccess = false;
   double tmpweight;
   int control = 0;

   do
   {
      weight = event.Generate();
      if (weight != weight)
         weight = 0;
      tmpweight = mRandom->Uniform(1.e-13);

      isSuccess = (weight != 0 && weight >= tmpweight);
      control++;
   } while (!(isSuccess || control > 1e6));
   return isSuccess;
}

bool CALM::FilterUnlikelyEvents_MINIJETS(TGenPhaseSpace &event0, TGenPhaseSpace &event1, double &weight0, double &weight1)
{
   bool isSuccess = false;
   double tmpweight;
   int control = 0;
   do
   {
      weight0 = event0.Generate();
      weight1 = event1.Generate();
      if ((weight0 != weight0) || (weight1 != weight1))
      {
         weight1 = 0;
         weight0 = 0;
      }
      tmpweight = mRandom->Uniform(1.e-13);

      isSuccess = (weight0 != 0 && weight1 != 0 && weight0 * weight1 >= tmpweight);
      control++;
   } while (!(isSuccess || control > 1e6));
   return isSuccess;
}

bool CALM::ReggaeNegativeEnergyCheck(int Nsum, double *masses, vector4 en, vector4 *avec)
{
   long int seed = time(NULL);

   bool checkE = true; //isSuccess; check negative energy
   int control = 0;
   do
   {
      Mconserv(en, Nsum, masses, avec, &seed); //genbod algoritmus
      collision(Nsum, avec, &seed);            //collision algoritmus
      checkE = true;

      //*****************************************
      //check particles for negative energy: in such case re-generate the event

      for (int i = 0; i < Nsum; i++)
      {
         if (avec[i][0] <= 0 || avec[i][1] != avec[i][1] || avec[i][2] != avec[i][2] || avec[i][3] != avec[i][3])
         {
            checkE = false; //cout<<"Negative energy!"<<endl;
         }
      }

      control++;
   } while (!(checkE || control >= 10));

   return checkE;
}

bool CALM::ReggaeNegativeEnergyCheck_MINIJETS(int Nsum, vector<double> *masses, vector4 en, vector4 *avec0, vector4 *avec1)
{
   double masses0[masses[0].size()];
   double masses1[masses[1].size()];
   for (int j = 0; j < masses[0].size(); ++j)
   {
      masses0[j] = masses[0][j];
   }
   for (int j = 0; j < masses[1].size(); ++j)
   {
      masses1[j] = masses[1][j];
   }

   long int seed = time(NULL);

   bool checkE = 1; //isSuccess; check negative energy
   int control = 0;
   do
   {
      checkE = 1;
      //first jet
      Mconserv(en, masses[0].size(), masses0, avec0, &seed); //genbod algoritmus
      collision(masses[0].size(), avec0, &seed);             //collision algoritmus
      //second jet
      Mconserv(en, masses[1].size(), masses1, avec1, &seed); //genbod algoritmus
      collision(masses[1].size(), avec1, &seed);             //collision algoritmus

      //*****************************************
      //check particles for negative energy: in such case re-generate the event

      for (int i = 0; i < masses[0].size(); i++)
      {
         if (avec0[i][0] <= 0 || avec0[i][1] != avec0[i][1] || avec0[i][2] != avec0[i][2] || avec0[i][3] != avec0[i][3])
         {

            checkE = 0; //cout<<"Negative energy!"<<endl;
         }
      }
      for (int i = 0; i < masses[1].size(); i++)
      {
         if (avec1[i][0] <= 0 || avec1[i][1] != avec1[i][1] || avec1[i][2] != avec1[i][2] || avec1[i][3] != avec1[i][3])
         {

            checkE = 0; //cout<<"Negative energy!"<<endl;
         }
      }

      control++;
   } while (!(checkE || control >= 10));

   return checkE;
}

void CALM::SaveAllParticles_GLOBAL(int Nsum, double weight, double **XYZrand, TGenPhaseSpace event, ParticleDB *aPartDB, list<Particle> *aParticles)
{
   Particle *tParticle;
   TLorentzVector *tmp;
   for (int i = 0; i < Nsum; i++)
   {
      tmp = event.GetDecay(i);
      tParticle = new Particle(aPartDB->GetParticleType(mParticlesThisEvent[i].c_str()));
      tParticle->SetParticlePX(tmp->E(), tmp->Px(), tmp->Py(), tmp->Pz(),
                               0, XYZrand[i][0], XYZrand[i][1], XYZrand[i][2],
                               weight, 0);
      aParticles->push_back(*tParticle);
      PRINT_DEBUG_2(mParticlesThisEvent[i] << " , " << endl);

      delete tParticle;
   }
}

void CALM::SaveAllParticles_MINIJETS(vector<double> *masses, vector<string> *names, double weight0, double weight1, double TotEnergy, double *divideEn, double **XYZrand, TGenPhaseSpace event0, TGenPhaseSpace event1, ParticleDB *aPartDB, list<Particle> *aParticles, eEventType aEventType)
{
   double phi, eta, theta, p1[3], p2[3], Ejet1, Ejet2;
   phi = mRandom->Uniform(0, 2 * TMath::Pi());
   eta = mRandom->Uniform(-2., 2.);
   theta = 2 * TMath::ATan(TMath::Exp(-eta));

   Ejet1 = TotEnergy * (divideEn[1] / (2. * (divideEn[0] + divideEn[1]))) / masses[0].size();
   Ejet2 = TotEnergy * (divideEn[1] / (2. * (divideEn[0] + divideEn[1]))) / masses[1].size();
   if (aEventType == MINIJETS_GLOBAL)
   {
      p1[0] = Ejet1 * TMath::Sin(theta) * TMath::Sin(phi);
      p1[1] = Ejet1 * TMath::Sin(theta) * TMath::Cos(phi);
      p1[2] = Ejet1 * TMath::Cos(theta);
      p2[0] = Ejet2 * TMath::Sin(theta) * TMath::Sin(phi);
      p2[1] = Ejet2 * TMath::Sin(theta) * TMath::Cos(phi);
      p2[2] = Ejet2 * TMath::Cos(theta);
   }
   else if (aEventType == MINIJETS_LOCAL)
   {
      p1[0] = TotEnergy / 4. / masses[0].size() * TMath::Sin(theta) * TMath::Sin(phi);
      p1[1] = TotEnergy / 4. / masses[0].size() * TMath::Sin(theta) * TMath::Cos(phi);
      p1[2] = TotEnergy / 4. / masses[0].size() * TMath::Cos(theta);
      p2[0] = TotEnergy / 4. / masses[1].size() * TMath::Sin(theta) * TMath::Sin(phi);
      p2[1] = TotEnergy / 4. / masses[1].size() * TMath::Sin(theta) * TMath::Cos(phi);
      p2[2] = TotEnergy / 4. / masses[1].size() * TMath::Cos(theta);
   }

   Particle *tParticle;
   TLorentzVector *tmp;
   for (int i = 0; i < masses[0].size(); i++)
   {
      tmp = event0.GetDecay(i);
      tParticle = new Particle(aPartDB->GetParticleType(names[0][i]));
      tParticle->SetParticlePX(tmp->E() + Ejet1, tmp->Px() + p1[0], tmp->Py() + p1[1], tmp->Pz() + p1[2],
                               0, XYZrand[i][0], XYZrand[i][1], XYZrand[i][2],
                               weight0 * weight1, 0);
      aParticles->push_back(*tParticle);
      delete tParticle;
   }
   for (int i = 0; i < masses[1].size(); i++)
   {
      tmp = event1.GetDecay(i);
      tParticle = new Particle(aPartDB->GetParticleType(names[1][i]));
      tParticle->SetParticlePX(tmp->E() + Ejet2, tmp->Px() - p2[0], tmp->Py() - p2[1], tmp->Pz() - p2[2],
                               0, XYZrand[masses[0].size() + i][0], XYZrand[masses[0].size() + i][1], XYZrand[masses[0].size() + i][2],
                               weight0 * weight1, 0);
      aParticles->push_back(*tParticle);
      delete tParticle;
   }
}

void CALM::SaveAllParticles_GLOBAL_REGGAE(int Nsum, vector4 *avec, double **XYZrand, ParticleDB *aPartDB, list<Particle> *aParticles)
{
   Particle *tParticle;
   TLorentzVector *tmp = new TLorentzVector();
   double weight = 1;

   for (int i = 0; i < Nsum; i++)
   {
      tmp->SetPxPyPzE(avec[i][1], avec[i][2], avec[i][3], avec[i][0]); //sat values from Reggae

      tParticle = new Particle(aPartDB->GetParticleType(mParticlesThisEvent[i].c_str()));
      tParticle->SetParticlePX(tmp->E(), tmp->Px(), tmp->Py(), tmp->Pz(),
                               0, XYZrand[i][0], XYZrand[i][1], XYZrand[i][2],
                               weight, 0);
      aParticles->push_back(*tParticle);
      PRINT_DEBUG_2(mParticlesThisEvent[i] << " , " << endl);

      delete tParticle;
   }
}

void CALM::SaveAllParticles_MINIJETS_REGGAE(vector<double> *masses, vector<string> *names, vector4 *avec0, vector4 *avec1, double TotEnergy, double **XYZrand, ParticleDB *aPartDB, list<Particle> *aParticles)
{

   double phi, eta, theta, p1[3], p2[3], Ejet1, Ejet2;
   phi = mRandom->Uniform(0, 2 * TMath::Pi());
   eta = mRandom->Uniform(-2., 2.);
   theta = 2 * TMath::ATan(TMath::Exp(-eta));

   p1[0] = TotEnergy / 4. / masses[0].size() * TMath::Sin(theta) * TMath::Sin(phi);
   p1[1] = TotEnergy / 4. / masses[0].size() * TMath::Sin(theta) * TMath::Cos(phi);
   p1[2] = TotEnergy / 4. / masses[0].size() * TMath::Cos(theta);
   Ejet1 = TotEnergy / 4. / masses[0].size();
   p2[0] = TotEnergy / 4. / masses[1].size() * TMath::Sin(theta) * TMath::Sin(phi);
   p2[1] = TotEnergy / 4. / masses[1].size() * TMath::Sin(theta) * TMath::Cos(phi);
   p2[2] = TotEnergy / 4. / masses[1].size() * TMath::Cos(theta);
   Ejet2 = TotEnergy / 4. / masses[1].size();

   Particle *tParticle;
   TLorentzVector *tmp = new TLorentzVector();

   for (int i = 0; i < masses[0].size(); i++)
   {
      tmp->SetPxPyPzE(avec0[i][1], avec0[i][2], avec0[i][3], avec0[i][0]); //sat values from Reggae
      tParticle = new Particle(aPartDB->GetParticleType(names[0][i]));
      tParticle->SetParticlePX(tmp->E() + Ejet1, tmp->Px() + p1[0], tmp->Py() + p1[1], tmp->Pz() + p1[2],
                               0, XYZrand[i][0], XYZrand[i][1], XYZrand[i][2],
                               1, 0);
      aParticles->push_back(*tParticle);
      delete tParticle;
   }
   for (int i = 0; i < masses[1].size(); i++)
   {
      tmp->SetPxPyPzE(avec1[i][1], avec1[i][2], avec1[i][3], avec1[i][0]); //sat values from Reggae
      tParticle = new Particle(aPartDB->GetParticleType(names[1][i]));
      tParticle->SetParticlePX(tmp->E() + Ejet2, tmp->Px() - p2[0], tmp->Py() - p2[1], tmp->Pz() - p2[2],
                               0, XYZrand[masses[0].size() + i][0], XYZrand[masses[0].size() + i][1], XYZrand[masses[0].size() + i][2],
                               1, 0);
      aParticles->push_back(*tParticle);
      delete tParticle;
   }

   delete tmp;
}

CALM::CALM() : mRandom(0), mNames(0), mNmean(0)
{
   Configurator *newConfig = new Configurator("./config.ini");
   newConfig->ReadParameters();
   eventConfig = new ConfigurationHolder(newConfig);

   mRandom = new TRandom2(0);

   //constants
   mNpart = 4; //particle types (pions, kaons, protons, lambdas)

   int Npartkinds[] = {3, 4, 4, 2}; //Look at Names below
   string Names[] = {
       "pi0139plu", "pi0139min", "pi0135zer",
       "Ka0492plu", "Ka0492min", "Ka0492zer", "Ka0492zrb",
       "pr0938plu", "pr0938plb", "ne0939zer", "ne0939zrb",
       "Lm1115zer", "Lm1115zrb"};


   int it = 0;
   mNpartkinds = new int[mNpart];
   mNames = new string *[mNpart];
   for (int i = 0; i < mNpart; i++)
   {
      mNpartkinds[i] = Npartkinds[i];
      mNames[i] = new string[Npartkinds[i]];
      for (int j = 0; j < Npartkinds[i]; j++)
      {
         mNames[i][j] = Names[it];
         PRINT_DEBUG_2("name[" << it << ":" << j << "," << i << "] = " << Names[it]);
         it++;
      }
   }

   //reading configuration
   mNmean = eventConfig->Nmean; //charged particle yields per rapidity unit
   mRapidityInterval = eventConfig->RapidityInterval;
   mXYZ = eventConfig->XYZ;
   if(eventConfig->pythiaMult==1){
      GetMultiplicitiesOfPartciles = &CALM::PythiaMult;
   }
   else{
      GetMultiplicitiesOfPartciles = &CALM::AlicePoissonMult;
   }
}
CALM::~CALM()
{
   delete mRandom;
}

int CALM::GenerateParticles(ParticleDB *aPartDB, int aMultBinMin, int aMultBinMax, double aEnergy, list<Particle> *aParticles, eEventType aEventType)
{
   int *Nrand; // multiplicities for each kind

   int** Npart; // particle to be generated
   
   int Nsum; // number of particles generated

   double **XYZrand; // vertex locations

   //_______distributing the total number of particles for each kind and for the specific particles
   Nrand = (this->*GetMultiplicitiesOfPartciles)(aMultBinMin, aMultBinMax, Nsum);

   //_______randomization type for particle with GLOBAL conservation
   Npart = GetTypesForParticles(Nrand, aPartDB);

   //________rewriting the particles into one list
   for (int i = 0; i < mNpart; ++i)
   {
      for (int j = 0; j < Nrand[i]; ++j)
      {
         mParticlesThisEvent.push_back(mNames[i][Npart[i][j]]);
      }
   }

   //________XYZ generating
   XYZrand = GetVertexXYZ(Nsum);

   //________Energy generating
   SetTotalEnergy(Nsum, aEnergy);

   //________Genbod part
   // generate total momentum for given energy
   double TotEnergy;
   switch (aEventType)
   {
   case GLOBAL:
   default:
   {
      TGenPhaseSpace event;
      double *masses = GetMasses(Nsum, aPartDB);

      if (!TrySetEventDecay(Nsum, masses, event, TotEnergy))
      {
         mParticlesThisEvent.clear();
         return 99;
      }

      double weight = 0;
      if (!FilterUnlikelyEvents(event, weight))
      {
         return 99;
      }

      SaveAllParticles_GLOBAL(Nsum, weight, XYZrand, event, aPartDB, aParticles);

      delete[] masses;
      break;
   }
   case MINIJETS_GLOBAL:
   {
      TGenPhaseSpace event1, event0;
      int amountOfJets = 2;
      vector<double> *masses = new vector<double>[amountOfJets];
      vector<string> *names = new vector<string>[amountOfJets];
      SeparateJets(Nsum, masses, names, aPartDB);

      double *divideEn = eventConfig->divideEn; // 0: energy of particles, 1: boostenergy

      if (!TrySetEventDecay_MINIJETS(Nsum, masses, event0, event1, TotEnergy, divideEn))
      {
         mParticlesThisEvent.clear();
         return 99;
      }

      double weight0 = 0;
      double weight1 = 1;
      if (!FilterUnlikelyEvents_MINIJETS(event0, event1, weight0, weight1))
      {
         return 99;
      }
      SaveAllParticles_MINIJETS(masses, names, weight0, weight1, TotEnergy, divideEn, XYZrand, event0, event1, aPartDB, aParticles, aEventType);
      delete[] masses;
      delete[] names;
      //delete[] divideEn;
      break;
   }
   case MINIJETS_LOCAL:
   {
      TGenPhaseSpace event1, event0;
      int amountOfJets = 2;
      vector<double> *masses = new vector<double>[amountOfJets];
      vector<string> *names = new vector<string>[amountOfJets];
      if (!SeparateJets_LOCAL(Nsum, masses, names, aPartDB))
      {
         mParticlesThisEvent.clear();
         return 99;
      }

      TLorentzVector *tmp;
      TLorentzVector en;
      double *divideEn = eventConfig->divideEn; //new double[2]{1, 1}; // 0: energy of particles, 1: boostenergy
      if (!TrySetEventDecay_MINIJETS(Nsum, masses, event0, event1, TotEnergy, divideEn))
      {
         mParticlesThisEvent.clear();
         return 99;
      }

      double weight0 = 0;
      double weight1 = 1;
      if (!FilterUnlikelyEvents_MINIJETS(event0, event1, weight0, weight1))
      {
         return 99;
      }

      SaveAllParticles_MINIJETS(masses, names, weight0, weight1, TotEnergy, divideEn, XYZrand, event0, event1, aPartDB, aParticles, aEventType);
      delete[] masses;
      delete[] names;
      //delete[] divideEn;
      break;
   }
   case GLOBAL_REGGAE:
   {
      //*********** REGGAE part ****************
      //*** good, but check negative energy

      vector4 en;
      vector4 *avec = new vector4[Nsum];

      double *masses = GetMasses(Nsum, aPartDB);

      // get total momentum
      TotEnergy = Etot; //include mass of the particles in the range

      //set starting values to distribute
      en[0] = TotEnergy;
      en[1] = 0.0;
      en[2] = 0.0;
      en[3] = 0.0; //0 - energy, 1- px, 2-py, 3-px

      if (!ReggaeNegativeEnergyCheck(Nsum, masses, en, avec))
      {
         return 99;
      }

      SaveAllParticles_GLOBAL_REGGAE(Nsum, avec, XYZrand, aPartDB, aParticles);

      delete[] masses;
      delete[] avec;
      break;
   }
   case MINIJETS_GLOBAL_REGGAE:
   {

      int amountOfJets = 2;
      vector<double> *masses = new vector<double>[amountOfJets];
      vector<string> *names = new vector<string>[amountOfJets];

      SeparateJets(Nsum, masses, names, aPartDB);

      vector4 en;
      vector4 *avec0 = new vector4[masses[0].size()];
      vector4 *avec1 = new vector4[masses[1].size()];
      // get total momentum
      TotEnergy = Etot; //include mass of the particles in the range

      //set starting values to distribute
      en[0] = TotEnergy / 4.;
      en[1] = 0.0;
      en[2] = 0.0;
      en[3] = 0.0; //0 - energy, 1- px, 2-py, 3-px
      if (!ReggaeNegativeEnergyCheck_MINIJETS(Nsum, masses, en, avec0, avec1))
      {
         return 99;
      }

      SaveAllParticles_MINIJETS_REGGAE(masses, names, avec0, avec1, TotEnergy, XYZrand, aPartDB, aParticles);

      delete[] masses;
      delete[] names;
      delete[] avec0;
      delete[] avec1;
      break;
   }
   case MINIJETS_LOCAL_REGGAE:
   {

      int amountOfJets = 2;
      vector<double> *masses = new vector<double>[amountOfJets];
      vector<string> *names = new vector<string>[amountOfJets];

      if (!SeparateJets_LOCAL(Nsum, masses, names, aPartDB))
      {
         mParticlesThisEvent.clear();
         return 99;
      }

      vector4 en;
      vector4 *avec0 = new vector4[masses[0].size()];
      vector4 *avec1 = new vector4[masses[1].size()];
      // get total momentum
      TotEnergy = Etot; //include mass of the particles in the range

      //set starting values to distribute
      en[0] = TotEnergy / 4.;
      en[1] = 0.0;
      en[2] = 0.0;
      en[3] = 0.0; //0 - energy, 1- px, 2-py, 3-px
      if (!ReggaeNegativeEnergyCheck_MINIJETS(Nsum, masses, en, avec0, avec1))
      {
         return 99;
      }
      SaveAllParticles_MINIJETS_REGGAE(masses, names, avec0, avec1, TotEnergy, XYZrand, aPartDB, aParticles);
      delete[] masses;
      delete[] names;
      delete[] avec0;
      delete[] avec1;
      break;
   }
   }
   delete[] Nrand;
   for(int i = 0; i < mNpart; i++)
      delete[] Npart[i];
   delete[] Npart;
   for(int i = 0; i < Nsum; i++)
      delete[] XYZrand[i];
   delete[] XYZrand;
   mParticlesThisEvent.clear();
   return 0;
}
