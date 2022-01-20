import numpy as np
from itertools import product
# from matplotlib import pyplot as plt

def get_subsystem_states(minstateval, maxstateval, statestepval, precision):
	'''
		inputs: 
			minstateval: (float) minimum value that the subsystem state can assume,
			maxstateval: (float) maximum value that the subsystem state can assume,
			statestepval: (float) step between adjacent values in the subsystem state space,
			precision: (int) float precision with which to record the values in the state space
		outputs: 
			all_states, (nparray) representing the subsystem's entire discrete state space
	'''
	dim = int((maxstateval - minstateval)/statestepval) + 1
	all_states = np.around(np.linspace(minstateval, maxstateval, dim), precision)
	return all_states

def get_bin_edges(ss, ss_step):
	'''
		inputs: 
			ss: (list) representing a discrete state space,
			ss_step: (float) statestepval
		outputs: 
			bin_edges: (nparray) 
	'''
	bin_edges = np.zeros(ss.shape[0] + 1)

	# centered bins
	for i, state in enumerate(ss):
		bin_edges[i] = state - ss_step/2

	bin_edges[-1] = ss[-1] + ss_step/2

	return bin_edges

def get_joint_states(list_of_sub_states):
	'''
		inputs: 
			list_of_sub_states: (list) of nparrays representing state spaces for a set of subsystems
		outputs: 
			jss: (list) of lists representing the joint state space (jss) for that set of subsystems
	'''
	joint_states = product(*list_of_sub_states)
	jss = [list(state) for state in joint_states] 
	return jss

def get_marginal_joint_states(subset, all_sub_states, suborder):
	'''
		inputs: 
			subset: (list) of subsystems (strings) for which the joint state space is desired,
				listed in an order not necessarily respecting suborder
			all_sub_states: (list) of nparrays representing state spaces for *all* subsystems, 
				listed in an order respecting suborder
			suborder: (dict) mapping each subsystem (string) to an integer, representing the 
				order that will be used for all variables and lists throughout the code
		outputs: 
			mjss: (list) of lists representing the marginal joint state space (mjss) for the 
				subsystems in subset
	'''
	if len(subset) == 1:
		mjss = [[i] for i in all_sub_states[suborder[subset[0]]]]
	else:
		mjss = get_joint_states([all_sub_states[suborder[s]] for s in subset])

	return mjss

def get_evolution_times(starttime, endtime, tstep, precision):
	'''
		inputs: 
			starttime: (float) start time,
			endtime: (float) end time,
			tstep: (float) timestep,
			precision: (int) float precision with which to record the time values
		outputs: 
			(nparray) representing the evolution times spanning starttime to endtime
	'''
	return get_subsystem_states(starttime, endtime, tstep, precision)