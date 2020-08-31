#include "ConfigurationHolder.h"
#include <THGlobal.h>
#include <TDatime.h>

//string tokenizer, used to separate values from string arrays defined in config.ini file
vector<string> ConfigurationHolder::SplitString(string strarray, char token){
  vector<string> vresult;
  int index;
  while((index = strarray.find(token)) != string::npos){
    string strvalue = strarray.substr(0, index);
    strarray.erase(0, index + 1);
    vresult.push_back(strvalue);
  }
  vresult.push_back(strarray);

  return vresult;
}


ConfigurationHolder::ConfigurationHolder()
: Nmean(new double[4]{1.493,0.183,0.083,0.048}), RapidityInterval(5), XYZ(new double[3]{5.,5.,5.}), pythiaMult(0),
  pionsMultDistr("0.334508*TMath::Gaus(x,56.8221,23.5326)*(5.97354e-07*x*x*x-8.88401e-05*x*x+0.00434252*x-0.0274243)"), pionsMultDistr_xMin(8), pionsMultDistr_xMax(150),
  kaonsMultDistr("0.731705*TMath::Gaus(x,15.5239,8.95871)*(1.1963e-05*x*x*x-0.000584791*x*x+0.010377*x-0.00451733)"), kaonsMultDistr_xMin(1), kaonsMultDistr_xMax(50),
  nucleonsMultDistr("1.37498*TMath::Gaus(x,8.86527,6.11529)*(2.80839e-05*x*x*x-0.00103983*x*x+0.0134321*x-0.0026861)"), nucleonsMultDistr_xMin(1), nucleonsMultDistr_xMax(50),
  lambdasMultDistr("3.21961*TMath::Gaus(x,1.14316,1.58606)*(0.00221997*x*x*x-0.0125258*x*x+0.00457262*x+0.118927)"), lambdasMultDistr_xMin(0), lambdasMultDistr_xMax(15),
  singleEnergyDistr("0.922477*(TMath::Power(x+2.15717,-1.57383)-1.40499e-05)"), singleEnergyDistr_xMin(0.4), singleEnergyDistr_xMax(1100),
  /*EtotMax(7000),*/ divideEn(new double[2]{1,1})
{
    TDatime tDate;
    
    tDate.Set();
    PRINT_MESSAGE("["<<tDate.AsSQLString()<<"]\tDefault configuration loaded");
}

ConfigurationHolder::ConfigurationHolder(Configurator *config){
    TDatime tDate;
    
    tDate.Set();
    PRINT_MESSAGE("["<<tDate.AsSQLString()<<"]\tLodaing configuration");

    try{
      string s = config->GetParameter("Nmean").Data();
      vector<string> vNmean = SplitString(config->GetParameter("Nmean").Data(), ',');
      if(vNmean.size()!=4){
        throw *(new TString("Nmean"));
      }
      Nmean = new double[4]{
          stod(vNmean[0]),
          stod(vNmean[1]),
          stod(vNmean[2]),
          stod(vNmean[3])

      };
      RapidityInterval = config->GetParameter("RapidityInterval").Atoi();
      vector<string> vXYZ = SplitString(config->GetParameter("XYZ").Data(), ',');
      if(vXYZ.size()!=3){
        throw *(new TString("XYZ"));
      }
      XYZ = new double[3]{
          stod(vXYZ[0]),
          stod(vXYZ[1]),
          stod(vXYZ[2])
      };

      pythiaMult = config->GetParameter("pythiaMult").Atoi();
      if(pythiaMult == 1)
        PRINT_MESSAGE("["<<tDate.AsSQLString()<<"]\tUsing Pythia particle distribution functions");

      pionsMultDistr = config->GetParameter("pionsMultDistr");
      pionsMultDistr_xMin = stod(config->GetParameter("pionsMultDistr_xMin").Data());
      pionsMultDistr_xMax = stod(config->GetParameter("pionsMultDistr_xMax").Data());

      kaonsMultDistr = config->GetParameter("kaonsMultDistr");
      kaonsMultDistr_xMin = stod(config->GetParameter("kaonsMultDistr_xMin").Data());
      kaonsMultDistr_xMax = stod(config->GetParameter("kaonsMultDistr_xMax").Data());

      nucleonsMultDistr = config->GetParameter("nucleonsMultDistr");
      nucleonsMultDistr_xMin = stod(config->GetParameter("nucleonsMultDistr_xMin").Data());
      nucleonsMultDistr_xMax = stod(config->GetParameter("nucleonsMultDistr_xMax").Data());

      lambdasMultDistr = config->GetParameter("lambdasMultDistr");
      lambdasMultDistr_xMin = stod(config->GetParameter("lambdasMultDistr_xMin").Data());
      lambdasMultDistr_xMax = stod(config->GetParameter("lambdasMultDistr_xMax").Data());


      singleEnergyDistr = config->GetParameter("singleEnergyDistr");
      singleEnergyDistr_xMin = stod(config->GetParameter("singleEnergyDistr_xMin").Data());
      singleEnergyDistr_xMax = stod(config->GetParameter("singleEnergyDistr_xMax").Data());

      //EtotMax = config->GetParameter("EtotMax").Atoi();
      
      vector<string> vdivideEn = SplitString(config->GetParameter("divideEn").Data(), ',');
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

ConfigurationHolder::~ConfigurationHolder(){
  delete [] Nmean;
  delete [] XYZ;
  delete [] divideEn;
}

