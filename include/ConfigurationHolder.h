#ifndef _TH2_CONFIGURATIONHOLDER_H_
  #define _TH2_CONFIGURATIONHOLDER_H_

#include <iostream>
#include <string>
#include "Configurator.h"

class ConfigurationHolder {
	
	public:
        ConfigurationHolder();
        ConfigurationHolder(Configurator *sMainConfig);

        int pythiaMult;
        std::string singleEnergyDistr;
        int Etot;
        double* divideEn;
        //int EventsPerFile;

};

#endif
