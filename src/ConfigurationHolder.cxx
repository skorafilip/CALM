#include "ConfigurationHolder.h"
#include <THGlobal.h>
#include <TDatime.h>


ConfigurationHolder::ConfigurationHolder()
: pythiaMult(0), singleEnergyDistr("0.922477*(TMath::Power(x+2.15717,-1.57383)-1.40499e-05)"), 
  MinEtot(7000), divideEn(new double[2]{1,1})//, EventsPerFile(10000)
{
}

ConfigurationHolder::ConfigurationHolder(Configurator *sMainConfig){
    TDatime tDate;
    
    tDate.Set();
    PRINT_MESSAGE("["<<tDate.AsSQLString()<<"]\tLodaing configuration.");

    try{
      pythiaMult = sMainConfig->GetParameter("pythiaMult").Atoi();
      singleEnergyDistr = sMainConfig->GetParameter("singleEnergyDistr");
      MinEtot = sMainConfig->GetParameter("MinEtot").Atoi();
      divideEn = new double[2]{
          sMainConfig->GetParameter("divideEn_energyofparticles").Atoi(),
          sMainConfig->GetParameter("divideEn_boostenergy").Atoi()
      };
      //EventsPerFile = sMainConfig->GetParameter("EventsPerFile").Atoi();
    }
    catch(TString s){
      PRINT_MESSAGE("Lack of parameters ("<<s<<"), program aborted!");
      exit(_ERROR_CONFIG_PARAMETER_NOT_FOUND_);
    }

}
