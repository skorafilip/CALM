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
        double pionsMultDistr_xMin;
        double pionsMultDistr_xMax;
        string kaonsMultDistr;
        double kaonsMultDistr_xMin;
        double kaonsMultDistr_xMax;
        string nucleonsMultDistr;
        double nucleonsMultDistr_xMin;
        double nucleonsMultDistr_xMax;
        string lambdasMultDistr;
        double lambdasMultDistr_xMin;
        double lambdasMultDistr_xMax;

        double* Nmean;
        double RapidityInterval;
        double* XYZ;

        string singleEnergyDistr;
        double singleEnergyDistr_xMin;
        double singleEnergyDistr_xMax;

        int MinEtot;
        double* divideEn;

        //int EventsPerFile;

};

#endif
