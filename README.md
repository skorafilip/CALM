# CALM
ConservAtion Laws Model

Authors: M. A. Janik, A. Zaborowska, P. Modzelewski, F. Skóra

Latest version: 1.0

Helpers
------------------------------
   * CALM-manual.pdf - manual with more specific description of how CALM works.
   * Thesis_PL.pdf - thesis describing latest modifications and publication of CALM (Polish language).
   * Documentation directory - you can find CALM documentation written in HTML (using doxygen tool) inside this directory, please open the Documentation/html/hierarhy.html file in your browser to see full page.



Quick start manual
------------------------------
Depedencies:

   * C++ compiler
   * Cmake
   * ROOT (http://root.cern.ch/)

Installation:

      mkdir build
      cd build
      cmake ..
      make
   
   
   If cmake fails to find ROOT package, you can specify its path via:
   
      cmake -DROOTSYS='your_path_to_root_installation_e_g_/opt/root' ..
      
Capabilities
------------------------------
CALM simulates proton-proton collisions and saves their results into files. It produces 13 types of hadrons, we can distinguish them into 4 kinds:
   1. Pions (π<sup>+</sup>, π<sup>0</sup>, π<sup>-</sup>)
   2. Kaons (K<sup>+</sup>, K<sup>0</sup>, anti-K<sup>0</sup> , K<sup>-</sup>)
   3. Nukleons (p, anti-p, n, anti-n)
   4. Lambdas (Λ, anti-Λ)

Below you can see exemple distribution of multiplicity of all particles types in one event:
![MultPID](https://raw.githubusercontent.com/skorafilip/CALM/master/images/GLOBAL_REGGAE_hevmultPIDPipPip.png)

CALM uses two Monte Carlo generators (Reggae and GENBOD). They can simulate regular collisions (Global) or you may order them to simulate minijet production process with of without (MinijetsLocal or Minijets) conservation laws check for each jet. You can see the correlation functions for all 6 option below.
![Correlation](https://raw.githubusercontent.com/skorafilip/CALM/master/images/CorrelationsPipPim.png)

You can analyze CALM results by your own or by using [tpi program](https://github.com/majanik/tpi_CALM).
