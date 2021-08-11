from math import exp, sqrt
from tools import *

# ************************************* Parameters *************************************************

# note that kT = 1
# parameters needed for rate matrices
# values are taken from Hartrich, Seifert, 2016 paper on sensor capacity

# defining parameters for subsystem state spaces
lmax, Nr1, Nm, Nr2 = 1, 3, 5, 3 

# ss_mins_LRM = [-1, 0, 0, 0]
# ss_maxes_LRM = [lmax, Nr1, Nm, Nr2]
# ss_steps_LRM = [0.25, 1, 1, 1]
# dl = ss_steps_LRM[0]

ss_mins_LRM = [0, 0, -1, 0]
ss_maxes_LRM = [Nm, Nr1, lmax, Nr2]
ss_steps_LRM = [1, 1, 0.25, 1]
dl = ss_steps_LRM[2]

# energy in kT required to bind a ligand to a receptor
# energy in kT required to process an ATP
konoff_ratio = 1.0
s_0 = 1/konoff_ratio
DeltaF = lambda l: l + np.log(konoff_ratio*s_0)
# DeltaF_ATP = 4

# rates
w_l, v_r = 1, 10
w_r = v_r*w_l/2

D_l = 0.1
D_r1 = 4*w_r/Nr1
B_r1 = D_r1/D_l

# print(sqrt(1 + ((v_r)**2)/(B_r1)))

# v_mminus = sqrt(1 + ((v_r)**2)/(B_r1))
# print(v_mminus)
v_mminus, v_mplus = .1, .2
kminus, kplus = .3, .2 
# v_mplus = ((kminus*v_mminus)/(kplus)) * exp(DeltaF_ATP)
# print(v_mplus)

# *********************************** Subsystems and Units *****************************************

# subsystem name : index for each subsystem in system
# suborder_LRM = {"l": 0, "r1": 1, "m": 2, "r2": 3}
suborder_LRM = {"m": 0, "r1": 1, "l": 2, "r2": 3}

# list of units in considered unit structure
units_LRM = [["l"],["l", "r1", "m"],["l", "r2"], ["l", "r1"], ["l", "r1", "r2"]]
# units_LRM2 = [["l"],["l", "r1", "m"],["l", "r2"],["l", "r1"]]

# pet names for units : index in unit structure
uind = {'phi': 0, 'omega': 1, 'alpha': 2, 'beta': 3, 'psi': 4}
# uind2 = {'phi': 0, 'omega': 1, 'alpha': 2, 'beta': 3}

# ******************* State Spaces for Subsystems, Units, and Total System *************************

ss_iter = range(len(suborder_LRM))
all_ss_LRM = [get_subsystem_states(ss_mins_LRM[i], 
									ss_maxes_LRM[i], 
									ss_steps_LRM[i], 5) for i in ss_iter]
bin_edges_LRM = [get_bin_edges(state_space, subsystem_step) 
					for state_space, subsystem_step in zip(all_ss_LRM, ss_steps_LRM)]

# generate joint state space for full system
jss_LRM = get_joint_states(all_ss_LRM)

# generate state spaces for units in the unit structure
uss_LRM = [get_marginal_joint_states(u, all_ss_LRM, suborder_LRM) for u in units_LRM]
# uss_LRM2 = [get_marginal_joint_states(u, all_ss_LRM, suborder_LRM) for u in units_LRM2]

# print('num joint states', len(jss_LRM))

# ******************************** Time Evolution Parameters ***************************************

# times to evolve the prob dists over
start, end, timestep = 0, 0.055000, 0.000050
times_LRM = get_evolution_times(start, end, timestep, 8)
# print(times_LRM)

# ********************************** Rate Matrix Functions *****************************************

# rate matrix functions for each of the subsystems
	# each of the functions take as input the joint state values, js
	# and return as output the value of the rate matrix element for the transition from
	#	js --> [s + key*ss_step(s) for s in js if this rmf corresponds to s else s]
	# question - should the DeltaF be multiplied by js["l"] or (js["l"] - js["r1"])?
rmf_l = {	1: lambda js: (D_l/(dl**2)) * exp(-((w_l * js["l"])/(2*D_l)) * dl),
    		-1: lambda js: (D_l/(dl**2)) * exp(+((w_l * js["l"])/(2*D_l)) * dl)}
rmf_r1 = {	1: lambda js: (w_r*exp(DeltaF(js["l"]))) * (Nr1 - js["r1"]),
    		-1: lambda js: (w_r * (js["r1"]))}
rmf_m = {	1: lambda js: (js["r1"]*kplus + v_mminus) * (Nm - js["m"]),
    		-1: lambda js: (js["r1"]*kminus + v_mplus) * (js["m"])}
rmf_r2 = {	1: lambda js: (w_r*exp(DeltaF(js["l"]))) * (Nr2 - js["r2"]),
    		-1: lambda js: (w_r * (js["r2"]))}
rmfs_LRM = [rmf_m, rmf_r1, rmf_l, rmf_r2]

