1D solution to the PNP equations for a binary electrolyte between parallel electrodes placed at z = 0 and z = H. The left electrode is powered with potential phi0 (dc and ac options are available.)

How to cite:

  Oscillating Electric Fields in Liquids Create a Long-Range Steady Field  
  Aref Hashemi, Scott C. Bukosky, Sean P. Rader, William D. Ristenpart, and Gregory H. Miller  
  Phys. Rev. Lett. 121, 185504, 2018

  Asymmetric rectified electric fields between parallel electrodes: Numerical and scaling analyses  
  Aref Hashemi, Gregory H. Miller, and William D. Ristenpart  
  Phys. Rev. E 99, 062603, 2019
  
  Net Responses in Nonlinear Dynamical Systems  
  Aref Hashemi  
  PhD dissertation, University of California Davis, 2021

Most important parameters are given in userINPUT.m:

     (parameter) description  

     (phi0)     applied potential on the electrodes at z = 0 (kBT/e)
     (ac)       [0 for dc, 1 for ac]
     (freq)     frequency (Hz) [will be ignored for dc]
     (x_M)      electrolyte strength (M)
     (D1)       cation diffusivity (m^2/s)
     (D2)       anion diffusivity (m^2/s)
     (q1)       cation charge number (-)
     (q2)       anion charge number (-)
     (epsinf)   dielectric constant
     (H)        gap (m)
     (h_S)      thickness of the Stern layer (m)
     (epsRatio) Stern to diffuse layer permittivity ratio (-)

Other parameters are available in INPUT.m.

The ouputs are V_phi, V_n1, V_n2, E

    V_phi, V_n1, V_n2 are (PM, JM_f) matrices (cell centered nodes)  
    E is a (PM,JM_f+1) matrix (face centered nodes)

######## dc

For dc case, the simulation runs untill an equilibrium solution is achieved. 

The script dcPlot_E_n called in the end of the MAIN file (commented by default) plots the last rows of E, V_n1, V_n2.

######## ac

For ac case, the simulation runs for cn_max cycles of the applied potential until a harmonic solution is achieved.

The script acPlot_E_n called in the end of the MAIN file (commented by default) plots the time-average E (AREF), along with 3D plots for E, V_n1, V_n2.

Scalings of the raw output:

    V_n1, V_n2 are scaled by n_inf (number concentration of the electrolyte).  
    V_phi is scaled by kBT/e.  
    All lengths including z_f (vector cell centered nodes), z_e (vector of face centered nodes) are scaled by h=H/2.  
    Electric field is scaled by k_BT/eh.  

Output files will be stored in directory ../matResults/DIR/, where DIR is

       ac=$ac-del=$delta-D1=$D1-phi0=$phi0-f(kHz)=$freq/1000-H(nm)=$H*1e9-c=$x_M


Correspondence should be addressed to:  
	       i) Aref Hashemi (aref@cims.nyu.edu)  
	       ii) Gregory H. Miller (grgmiller@ucdavis.edu)  
	       iii) William D. Ristenpart (wdristenpart@ucdavis.edu) 