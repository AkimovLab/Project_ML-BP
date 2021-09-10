# Do the NA-MD calculations

## 1. Prepare files

### 1.1. Copy the data `res_mb_sp.tar.bz2` from the folder `3_nacs` and unpack the folder

    tar -xjf res_mb_sp.tar.bz2


### 1.2. Edit the `run_namd.py` file

  The number of files, initial timestep (reflected in the names of the files in `res_mb_sp` directory), the number states and the active space

    params["nfiles"]         = 1999
    params["init_times"]     = [5000]
    params["nstates"]        = 11 # total number of electronic states
    params["active_space"]   = list(range(11)) # indexing is from 0!

  Make sure to use the correct number of states here:

    zero = MATRIX(11, 11)

  Make sure to use the correct number of files (timesteps) in the dephasing time calculations:

    #================== COMPUTE DEPHASING AND ENERGY GAPS  ===============
    params["init_times"] = [0]
    params["nsteps"] = 1999


  Change the control and model parameters as needed:

    #================== SET UP THE DYNAMICS AND DISTRIBUTED COMPUTING SCHEME  ===============
    params_nbra = { "nsteps":1999, "dt":1.0*units.fs2au, 
                    "ntraj":10, "x0":[-4.0], "p0":[4.0], "masses":[2000.0], "k":[0.01],                  
                    "nstates":11, "istate":[1, 1],
                    "which_adi_states":range(11), "which_dia_states":range(11),
                    "rep_ham":1, "tsh_method":0, 
                    "force_method":0, "nac_update_method":0,
                    "hop_acceptance_algo":31, "momenta_rescaling_algo":0,
                    "time_overlap_method":1,
                    "mem_output_level":-1,
                    "txt_output_level":3,
                    "properties_to_save": ['timestep', 'time', 'SH_pop', 'SH_pop_raw'],
                    "state_tracking_algo":2, "convergence":0,  "max_number_attempts":100, 
                    "min_probability_reordering":0.01,  
                    "decoherence_algo":0,
                    "decoherence_rates":rates,
                    "ave_gaps":gaps
                  }
    model_params_nbra = {"model":2, "nstates":11, "filename":None }


  In particular, pay attention to:

  - `"nsteps":1999` - should be consistent with the number of files
  - `"dt":1.0*units.fs2au` - should be what you used in the MD
  - `"nstates":11` (in both dictionaries) - should be consistent with the number of states you have elsewhere
  - `"which_adi_states":range(11)` and  `"which_dia_states":range(11)` - should be consistent with the number of states

  - `"txt_output_level":3` and `"properties_to_save": ['timestep', 'time', 'SH_pop', 'SH_pop_raw']` - is the minimal output you really needed
     you may want to request more properties to be saved, but beware of the disk requirements, since these files may become large, especially if
     you have many trajectories

  - `"ntraj":10` - this is one doesn't matter, since it is redefined later - be mindful of this number since it affects the speed of calculations  
     (coordinate it with the computing resources requested in the `submit.slm`)
 
  - `"istate":[1, 1]` - doesn't matter since we'll loop over it


  Finally, this is where the most interesting parameters are changed:

    nthreads = 100
    params_nbra["Temperature"]   = 300.0
    params_nbra["ntraj"]         = 2000
    methods = {0:"FSSH", 1:"IDA", 2:"mSDM", 3:"DISH", 21:"mSDM2", 31:"DISH2" }

    init_states = [1, 2]
    tsh_methods = [0, 1, 2, 3, 21, 31]
    batches = list(range(25))

  Again, be mindful how many threads you request - with too many of them, there may be some delays in scheduling.
  

### 1.3. Edit the SLURM `submit.slm` file

  Mainly pay attention to these parameters:

    #SBATCH --time=20:00:00
    #SBATCH --cpus-per-task=10
    #SBATCH --mem=100000

 With many threads requested, having some good chunk of memory allocated is a good idea.


## 2. Run the calculation

  This is typically done on cluster, but for light jobs can be easily run on a local computer

  For a cluster:

    sbatch submit.slm 

  On a local computer:

    python run_namd.py


## 3. Plot the results with `plot_dynamics.ipynb`

  This notebook will make figures and will store them in the `res` directory. 

  You may need to create it manually or delete it on the first run, if you get some errors
