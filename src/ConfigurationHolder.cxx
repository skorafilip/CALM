#include "ConfigurationHolder.h"
#include <THGlobal.h>
#include <TDatime.h>

//string tokenizer, used to separate values from string arrays defined in config.ini file
vector<string> SplitString(string strarray){
  vector<string> vresult;
  int index;
  while((index = strarray.find(",")) != string::npos){
    string strvalue = strarray.substr(0, index);
    strarray.erase(0, index + 1);
    vresult.push_back(strvalue);
  }
  vresult.push_back(strarray);

  return vresult;
}


/*ConfigurationHolder::ConfigurationHolder()
: pythiaMult(0), singleEnergyDistr("0.922477*(TMath::Power(x+2.15717,-1.57383)-1.40499e-05)"), 
  MinEtot(7000), divideEn(new double[2]{1,1})//, EventsPerFile(10000)
{
}*/

ConfigurationHolder::ConfigurationHolder(Configurator *sMainConfig){
    TDatime tDate;
    
    tDate.Set();
    PRINT_MESSAGE("["<<tDate.AsSQLString()<<"]\tLodaing configuration");

    try{
      string s = sMainConfig->GetParameter("Nmean").Data();
      vector<string> vNmean = SplitString(sMainConfig->GetParameter("Nmean").Data());
      if(vNmean.size()!=4){
        throw *(new TString("Nmean"));
      }
      Nmean = new double[4]{
          stod(vNmean[0]),
          stod(vNmean[1]),
          stod(vNmean[2]),
          stod(vNmean[3])

      };
      RapidityInterval = sMainConfig->GetParameter("RapidityInterval").Atoi();
      vector<string> vXYZ = SplitString(sMainConfig->GetParameter("XYZ").Data());
      if(vXYZ.size()!=3){
        throw *(new TString("XYZ"));
      }
      XYZ = new double[3]{
          stod(vXYZ[0]),
          stod(vXYZ[1]),
          stod(vXYZ[2])
      };

      pythiaMult = sMainConfig->GetParameter("pythiaMult").Atoi();
      if(pythiaMult == 1)
        PRINT_MESSAGE("["<<tDate.AsSQLString()<<"]\tUsing Pythia particle distribution functions");

      pionsMultDistr = sMainConfig->GetParameter("pionsMultDistr");
      pionsMultDistr_xMin = stod(sMainConfig->GetParameter("pionsMultDistr_xMin").Data());
      pionsMultDistr_xMax = stod(sMainConfig->GetParameter("pionsMultDistr_xMax").Data());

      kaonsMultDistr = sMainConfig->GetParameter("kaonsMultDistr");
      kaonsMultDistr_xMin = stod(sMainConfig->GetParameter("kaonsMultDistr_xMin").Data());
      kaonsMultDistr_xMax = stod(sMainConfig->GetParameter("kaonsMultDistr_xMax").Data());

      nucleonsMultDistr = sMainConfig->GetParameter("nucleonsMultDistr");
      nucleonsMultDistr_xMin = stod(sMainConfig->GetParameter("nucleonsMultDistr_xMin").Data());
      nucleonsMultDistr_xMax = stod(sMainConfig->GetParameter("nucleonsMultDistr_xMax").Data());

      lambdasMultDistr = sMainConfig->GetParameter("lambdasMultDistr");
      lambdasMultDistr_xMin = stod(sMainConfig->GetParameter("lambdasMultDistr_xMin").Data());
      lambdasMultDistr_xMax = stod(sMainConfig->GetParameter("lambdasMultDistr_xMax").Data());


      singleEnergyDistr = sMainConfig->GetParameter("singleEnergyDistr");
      singleEnergyDistr_xMin = stod(sMainConfig->GetParameter("singleEnergyDistr_xMin").Data());
      singleEnergyDistr_xMax = stod(sMainConfig->GetParameter("singleEnergyDistr_xMax").Data());

      EtotMax = sMainConfig->GetParameter("EtotMax").Atoi();
      
      vector<string> vdivideEn = SplitString(sMainConfig->GetParameter("divideEn").Data());
      if(vdivideEn.size()!=2){
        throw *(new TString("divideEn"));
      }
      divideEn = new double[2]{
          stoi(vdivideEn[0]),
          stoi(vdivideEn[1])
      };
    }
    catch(TString s){
      PRINT_MESSAGE("Lack of parameters or invalid format ("<<s<<"), program aborted!");
      exit(_ERROR_CONFIG_PARAMETER_NOT_FOUND_);
    }

}
