# Explanation of the data folders:

pbe_plus_u:    xTB MD with the PBE+U (and TD-DFT) on top of it - this one is used in the manuscrit

pbe:  dftb+ MD with PBE nacs (and TD-DFT) on top of it - underestimated gaps, presented here only as an example


# 1. First, run the SP and TD-DFPT calculations 

## 1.1. Prepare files

### 1.1.1. Edit run.py

   * trajectory file name
   * istep, fstep, njobs
   * cp2k template name

### 1.1.2. Edit submit_template.slm

   * number of hours, memory
   * results directory (absolute path)
   * minband, maxband, homo index - should be consistent with the cp2k template
   * the number of steps in the "python" script

### 1.1.3. Edit cp2k input template

   * make sure you create enough CUBE files
   * make sure you have enough added MOs
   * make sure the number of excited states is consistent with the one in the submit_template


## 1.2. Run the calculations

   * make sure to activate libra environment

   * run the job submission script  `python run.py`



# 2. Second, we want to compute the MB couplings

## 2.1. Option 1: Run the Jupyter notebook `nacs.ipynb` on a local computer

   It may take a while, but not too long (could be up to few hours). 

   It will produce embedded graphics at the end


## 2.2. Option 2: Run as a Python script on a cluster

   * Save the `nacs.ipynb` notebook as a Python file `nacs.py`

   * Run it on the cluster using the `submit_nacs.slm` script:

       `sbatch submit_nacs.slm`

   This option will not produce the embedded graphics and may also give
   an error about plotting the end figure, but this should not be a problem - we'll 
   get it in the next stage.


   In both cases, the file `step3_many_body.log` will be produced that contains some 
   information about the SD basis discovered

   
   In both cases, the folder `res_mb_sp` will be generated, with all files needed for the
   next steps.



