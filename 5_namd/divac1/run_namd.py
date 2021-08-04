
#================== BASIC SETUPS AND COMMON FUNCTIONS ===============
import os
import sys
import time
import math
import multiprocessing as mp
import numpy as np
import matplotlib.pyplot as plt


# Fisrt, we add the location of the library to test to the PYTHON path
from liblibra_core import *
from libra_py import units
import util.libutil as comn
from libra_py import data_visualize, dynamics_plotting
from libra_py import data_conv, data_stat, data_outs, data_read
import libra_py.workflows.nbra.decoherence_times as decoherence_times
import libra_py.workflows.nbra.step4 as step4
import libra_py.dynamics.tsh.compute as tsh_dynamics
import libra_py.dynamics.tsh.plot as tsh_dynamics_plot

HVIB_CI, ST_CI, HAM_CI, NAC_CI = None, None, None, None  # Global variables to speed up things

def get_avegage_energies(Hvib):
        
    nsteps = len(Hvib[0])
    nstates   = Hvib[0][0].num_of_cols

    # Make a list for the SD energies and populate it
    dE = [0.0]*nstates
        
    for step in range( nsteps ):
        for i in range( nstates ):
            dE[i] += Hvib[0][step].get( i, i ).real - Hvib[0][step].get( 0, 0 ).real 
            
    for i in range( nstates ):
        dE[i] = dE[i] / nsteps
        
    return dE

class tmp:
    pass


def compute_model_nbra(q, params, full_id):
    """   
    Read in the vibronic Hamiltonians along the trajectories    

    Args: 
        q ( MATRIX(1,1) ): coordinates of the particle, ndof, but they do not really affect anything
        params ( dictionary ): model parameters

            * **params["timestep"]** ( int ):  [ index of the file to read ]
            * **params["prefix"]**   ( string ):  [ the directory where the hdf5 file is located ]
            * **params["filename"]** ( string ):  [ the name of the HDF5 file ]
        
    Returns:       
        PyObject: obj, with the members:

            * obj.hvib_adi ( CMATRIX(n,n) ): adiabatic vibronic Hamiltonian 
            
    """
                                    
    hvib_adi, basis_transform, time_overlap_adi = None, None, None
  
    Id = Cpp2Py(full_id)
    indx = Id[-1]    
    timestep = params["timestep"]    
    #filename = params["filename"]
    nadi = params["nstates"]
    
    #x = 5000 + timestep
    
                            
    #============ Vibronic Hamiltonian ===========        
    #filename_re =  F"res_mb_sp/Hvib_ci_{x}_re"
    #filename_im =  F"res_mb_sp/Hvib_ci_{x}_im"
    #hvib_adi = data_read.get_matrix(nadi, nadi, filename_re, filename_im, list(range(nadi)), 1, 1)

    #=========== Basis transform, if available =====
    basis_transform = CMATRIX(nadi, nadi)            
    basis_transform.identity()
                                        
    #========= Time-overlap matrices ===================
    #filename_re =  F"res_mb_sp/St_ci_{x}_re"
    #filename_im =  F"res_mb_sp/St_ci_{x}_im"
    #time_overlap_adi = data_read.get_matrix(nadi, nadi, filename_re, filename_im, list(range(nadi)), 1, 1)
                
    obj = tmp()
    obj.ham_adi = HAM_CI[timestep]
    obj.nac_ci = NAC_CI[timestep]
    obj.hvib_adi = HVIB_CI[timestep] #hvib_adi
    obj.basis_transform = basis_transform
    obj.time_overlap_adi = ST_CI[timestep] #time_overlap_adi
    
    return obj


def compute_model(q, params, full_id):
    model = params["model"]
    res = None    
    if model==1:        
        pass
    elif model==2:
        res = compute_model_nbra(q, params, full_id)
    return res


def run_tsh(common_params, model_params, prefix):    

    params = dict(common_params)            

    # Random numbers generator object
    rnd = Random()
    
    #============ Initialize dynamical variables ==================
    x0, p0, masses, k0 = params["x0"], params["p0"], params["masses"], params["k"]
    ntraj, nstates = params["ntraj"], params["nstates"]
    
    # Nuclear
    init_nucl = {"init_type":2, "force_constant":k0, "ntraj":ntraj}
    #init_nucl = {"init_type":3, "force_constant":k0, "ntraj":ntraj}
    q, p, iM = tsh_dynamics.init_nuclear_dyn_var(x0, p0, masses, init_nucl, rnd)
    
    # Electronic
    istate = params["istate"]
    istates = []
    for i in range(nstates):
        istates.append(0.0)
    istates[ istate[1] ] = 1.0    
    _init_elec = { "init_type":3, "nstates":nstates, "istates":istates, "rep":istate[0],  "ntraj":ntraj   }
    
    #============= Dynamical variables ==============
    dyn_params = dict(common_params)
    
    # This should update only the properties that aren't defined, but not override the existing values!
    critical_params = [  ]     
    default_params = { "prefix":prefix, "mem_output_level":3 }     
    comn.check_input(dyn_params, default_params, critical_params)
                    
    _model_params = dict(model_params)
    _model_params.update({"model0": model_params["model"] })
    
    start = time.time()                               
    res = tsh_dynamics.generic_recipe(q, p, iM, dyn_params, compute_model, _model_params,_init_elec, rnd)
    end = time.time()    
    print(F"Calculation time = {end - start} seconds")


#================== READ IN THE DATA ===============

absolute_path = os.getcwd()
params = {}
params["data_set_paths"] = []
params["data_set_paths"].append(absolute_path+"/res_mb_sp/")

params["Hvib_re_prefix"] = "Hvib_ci_"; params["Hvib_re_suffix"] = "_re"
params["Hvib_im_prefix"] = "Hvib_ci_"; params["Hvib_im_suffix"] = "_im"
params["nfiles"]         = 1999
params["init_times"]     = [5000]
params["nstates"]        = 11 # total number of electronic states
params["active_space"]   = list(range(11)) # indexing is from 0!

# Include HOMO and up to the last electronic state
hvib = step4.get_Hvib2(params)
HVIB_CI = hvib[0]
print ("Length of hvib is: ", len(hvib[0]))

HAM_CI, NAC_CI= [], []
zero = MATRIX(11, 11)
for h in HVIB_CI:
    HAM_CI.append( CMATRIX(h.real(), zero) )
    NAC_CI.append( CMATRIX(zero, -1.0 * h.imag()) )


params["Hvib_re_prefix"] = "St_ci_"; params["Hvib_re_suffix"] = "_re"
params["Hvib_im_prefix"] = "St_ci_"; params["Hvib_im_suffix"] = "_im"
ST_CI = step4.get_Hvib2(params)[0]



#================== COMPUTE DEPHASING AND ENERGY GAPS  ===============

params["init_times"] = [0]
params["nsteps"] = 1999
# Compute energy gaps and decoherence times over part of the trajectory we want
tau, rates = decoherence_times.decoherence_times_ave(hvib, [params["init_times"][0]], params["nsteps"], 0)
dE         = decoherence_times.energy_gaps_ave(      hvib, [params["init_times"][0]], params["nsteps"])
avg_deco   = tau/units.fs2au


dE = get_avegage_energies(hvib)

print("Average energies above the GS are")
for indx, de in enumerate(dE):
    print(F"State {indx} energy {de * units.au2ev} eV")



nst = len(dE)
gaps = MATRIX(nst, nst)
for istate1 in range(nst):
    for istate2 in range(nst):        
        gaps.set(istate1, istate2,  math.fabs(dE[istate1] - dE[istate2]) )





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

#run_tsh(params_nbra, model_params_nbra, "nbra_fssh")



nthreads = 100
params_nbra["Temperature"]   = 300.0
params_nbra["ntraj"]         = 2000
methods = {0:"FSSH", 1:"IDA", 2:"mSDM", 3:"DISH", 21:"mSDM2", 31:"DISH2" }

init_states = [1, 2]
tsh_methods = [0, 1, 2, 3, 21, 31]
batches = list(range(25))

rnd = Random()

var_pool = []
for istate in init_states:  # initial states   
    for method in tsh_methods:  # decoherence method: FSSH, IDA, mSDM, DISH        
        for batch in batches:        
            wait = rnd.uniform(0.0, 25.0)
            var_pool.append( (istate, method, batch, wait) )



prefix = "namd_res"
os.mkdir(prefix)
os.system(F"cp run_namd.py {prefix}/run_namd_params")


def run_dynamics( istate, method, batch, wait_time ):

    prms = dict(params_nbra)

    prms["istate"]    = [1, istate] # Recall index from 0

    if method==0:  # FSSH
        prms["tsh_method"] = 0 # FSSH
        prms["decoherence_algo"] = -1  # no decoherence
        prms["dephasing_informed"] = 0

    elif method==1:  # IDA
        prms["tsh_method"] = 0 # FSSH
        prms["decoherence_algo"] = 1  # ID
        prms["dephasing_informed"] = 0

    elif method==2:  # mSDM
        prms["tsh_method"] = 0 # FSSH
        prms["decoherence_algo"] = 0  # mSDM
        prms["dephasing_informed"] = 0

    elif method==3:  # DISH
        prms["tsh_method"] = 3 # DISH
        prms["decoherence_algo"] = -1  # no other decoherence
        prms["dephasing_informed"] = 0

    elif method==21:  # mSDM
        prms["tsh_method"] = 0 # FSSH
        prms["decoherence_algo"] = 0  # mSDM
        prms["dephasing_informed"] = 1

    elif method==31:  # DISH
        prms["tsh_method"] = 3 # DISH
        prms["decoherence_algo"] = -1  # no other decoherence
        prms["dephasing_informed"] = 1


    name = methods[method]   
    time.sleep( int(wait_time) )

    run_tsh(prms, model_params_nbra, F"{prefix}/start_s{istate}_{name}_batch{batch}" )
    #step4.run(hvib, prms)



#================== DO THE CALCULATIONS  ===============
                        
t1 = time.time()

pool = mp.Pool( nthreads )
pool.starmap( run_dynamics, var_pool )
pool.join()
pool.close()
#pool.join()

t2 = time.time()

print(F"Total time {t2 - t1}")


