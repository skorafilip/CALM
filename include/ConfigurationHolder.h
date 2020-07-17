#ifndef _TH2_CONFIGURATIONHOLDER_H_
  #define _TH2_CONFIGURATIONHOLDER_H_

#include <iostream>
#include <string>
#include "Configurator.h"

//FILIPS: nowa klasa, ma na celu przechowywanie parametrów wczytanych w plikach konfiguracyjnych
//zakładam że gdy będą dochodzić nowe paramtery trzeba będzie dodawać nowe zmienne
//wszystko wczytuje się raz na początku działania programu, potem jest tylko odczyt z pól klasy, wydaje mi się, że to optymalne rozwiązanie
using namespace std;
class ConfigurationHolder {
	
	public:
        ConfigurationHolder();
        ConfigurationHolder(Configurator *sMainConfig);

        int pythiaMult;
        string pionsMultDistr;
        int pionsMultDistr_xMin;
        int pionsMultDistr_xMax;
        string kaonsMultDistr;
        int kaonsMultDistr_xMin;
        int kaonsMultDistr_xMax;
        string nucleonsMultDistr;
        int nucleonsMultDistr_xMin;
        int nucleonsMultDistr_xMax;
        string lambdasMultDistr;
        int lambdasMultDistr_xMin;
        int lambdasMultDistr_xMax;

        double* Nmean;
        double RapidityInterval;
        double* XYZ;

        string singleEnergyDistr;
        int MinEtot;
        double* divideEn;

        //int EventsPerFile;

};

#endif
