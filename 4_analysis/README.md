# Explanation of the data folders:

pbe_plus_u:    xTB MD with the PBE+U (and TD-DFT) on top of it - this is what is used in the manuscript

pbe:  dftb+ MD with PBE nacs (and TD-DFT) on top of it - underestimated gaps, this is presented here just as an example



1. Copy the results of step3 (NACs and TD-DFPT calculations) from step3_nacs folder here:

   * all_logfiles.tar.bz2 
   * all_pdosfiles.tar.bz2
   * res.tar.bz2
   * res_mb_sp.tar.bz2

   Unpack the files, e.g.:

   tar -xf all_logfiles.tar.bz2


2. Run the `analysis.ipynb` notebook, it will produce:

  * pDOS (pdosa.png, pdos_alp.txt, pdos_bet.txt)
  * absorption spectrum ( spectrum_divac.png - don't worry, the file is called this way in both folders, but that's just a name + absorption__spctrum.txt)
  * evolution of state energies and averaged NACs map ( energies_nacs-original.png )
  * dephasing times and the plot of the dephasing times for all pairs of states ( dephasing_times.png )
  * averaged energy gaps and energy gap distribution function ( gap_distributions.png )
  * gap fluctuation ACF and influence spectrum ( acf-IFS.png )


  This script is rather generic - so, likely you'll only need to change the titles of some figures generated


3. Additionally run the `check_NBRA.ipynb` notebook, it will produce:

  * the plots of the ground and excited state energies along the MD path ( nbra.png )

