#include "ConfigurationHolder.h"

ConfigurationHolder::ConfigurationHolder()
: pythiaMult(0), singleEnergyDistr("0.922477*(TMath::Power(x+2.15717,-1.57383)-1.40499e-05)"), 
  Etot(7000), divideEn(new double[2]{1,1})//, EventsPerFile(10000)
{
}

ConfigurationHolder::ConfigurationHolder(Configurator *sMainConfig){
    pythiaMult = sMainConfig->GetParameter("pythiaMult").Atoi();
    singleEnergyDistr = sMainConfig->GetParameter("singleEnergyDistr");
    Etot = sMainConfig->GetParameter("Etot").Atoi();
    divideEn = new double[2]{
        sMainConfig->GetParameter("divideEn_energyofparticles").Atoi(),
        sMainConfig->GetParameter("divideEn_boostenergy").Atoi()
    };
    //EventsPerFile = sMainConfig->GetParameter("EventsPerFile").Atoi();

}
