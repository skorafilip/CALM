#ifndef _TH2_CONFIGURATIONHOLDER_H_
  #define _TH2_CONFIGURATIONHOLDER_H_

#include <iostream>
#include <string>
#include "Configurator.h"

//FILIPS: nowa klasa, ma na celu przechowywanie parametrów wczytanych w plikach konfiguracyjnych
//zakładam że gdy będą dochodzić nowe paramtery trzeba będzie dodawać nowe zmienne
//wszystko wczytuje się raz na początku działania programu, potem jest tylko odczyt z pól klasy, wydaje mi się, że to optymalne rozwiązanie
class ConfigurationHolder {
	
	public:
        ConfigurationHolder();
        ConfigurationHolder(Configurator *sMainConfig);

        int pythiaMult;
        std::string singleEnergyDistr;
        int MinEtot;
        double* divideEn;
        //int EventsPerFile;

};

#endif
