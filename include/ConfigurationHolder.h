#ifndef _TH2_CONFIGURATIONHOLDER_H_
  #define _TH2_CONFIGURATIONHOLDER_H_

#include <iostream>
#include <string>
#include "Configurator.h"

using namespace std;
class ConfigurationHolder {
	
	public:
        ConfigurationHolder();
        ConfigurationHolder(Configurator *config);
        ~ConfigurationHolder();

        int customMult; ///< %Parameter that determines if custom distribution will be used 
        string pionsMultDistr; ///< Custom distribution function for pions 
        double pionsMultDistr_xMin; ///< Custom distribution function minimum range for pions 
        double pionsMultDistr_xMax; ///< Custom distribution function maximum range for pions 
        string kaonsMultDistr; ///< Custom distribution function for kaons 
        double kaonsMultDistr_xMin; ///< Custom distribution function minimum range for kaons 
        double kaonsMultDistr_xMax; ///< Custom distribution function maximum range for kaons 
        string nucleonsMultDistr; ///< Custom distribution function for nucleons 
        double nucleonsMultDistr_xMin; ///< Custom distribution function minimum range for nucleons 
        double nucleonsMultDistr_xMax; ///< Custom distribution function maximum range for nucleons 
        string lambdasMultDistr; ///< Custom distribution function for lambdas 
        double lambdasMultDistr_xMin; ///< Custom distribution function minimum range for lambdas 
        double lambdasMultDistr_xMax; ///< Custom distribution function maximum range for lambdas 

        double* Nmean; ///< Charged particle yields per rapidity unit from 900 GeV 
        double RapidityInterval; ///< Interval of rapidity 
        double* XYZ; ///< Sigma of Gaus distribution for 3D 

        string singleEnergyDistr; ///< Energy distribution function for all particles 
        double singleEnergyDistr_xMin; ///< Energy distribution function minimum range for for all particles 
        double singleEnergyDistr_xMax;///< Energy distribution function maximum range for for all particles 

        //int EtotMax; ///< Total energy of all particles
        double* divideEn; ///< Energy divider - divideEn[0]: energy of particles, divideEn[1]: boostenergy

  private:
        vector<string> SplitString(string strarray, char token);

};

#endif

/*! @file ConfigurationHolder.h
 * @brief Definition of ConfigurationHolder class. Holds values from configuration file.
 */
/*! @class ConfigurationHolder
 * @brief Holds values of parameters specified in configuration file. CALM sets values once, on the beggining, and they are being read many times during simulation.
 * 
 *
 * @fn ConfigurationHolder::ConfigurationHolder()
 * @brief Default constructor.
 *
 * @fn ConfigurationHolder::ConfigurationHolder(Configurator *config)
 * @brief Reads specified parameters from given Configurator.
 * @param [in] config pointer to Configurator
 *
 * If parameter was not found the parameter keyword will be thrown.
 * @param [in] aKeyword keyword
 * @exception TString requested keyword
 * 
 *  
 * @fn ConfigurationHolder::~ConfigurationHolder()
 * @brief Destructor.
 *
 * @fn vector<string> ConfigurationHolder::SplitString()
 * @brief String splitter, helps to read arrays from string parameters. Returnes vector of strings separated from strarray with token
 * @param [in] strarray array string to separate
 * @param [in] token which separates values
 * 
 * 
 */