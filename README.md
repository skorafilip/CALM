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
CALM simulates proton-proton collisions and saves their results into files. You can analyze them by your own or by using [tpi program](https://github.com/majanik/tpi_CALM).
