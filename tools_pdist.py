import numpy as np
# from LRM_properties import*
from scipy.stats import norm, binom, randint, uniform
from scipy.linalg import expm

def discretize_pdf(bin_edges, samples):
	'''
		inputs:
			bin_edges: (nparray) of edges for binning samples,
			samples: (nparray) of samples taken from the analogous continuous distribution
		outputs: 
			probs (nparray) representing the distribution over the discrete state space
	'''
	hist, edges = np.histogram(samples, bins=list(bin_edges))
	if sum(hist) == 0:
		probs = hist
	else:
		probs = hist/sum(hist)

	return probs

def discretize_subsystem_pdfs(all_bin_edges, all_samples):
	'''
		inputs: 
			all_bin_edges: (list) of edges for binning samples for a set of subsystems,
			all_samples: (list) of list of samples from each subsystem's continuous state spaces
		outputs: 
			all_dists: (list) of discrete distributions for the set of subsystems
	'''
	all_dists = []

	for bin_edges, samples in zip(all_bin_edges, all_samples):
		all_dists.append(discretize_pdf(bin_edges, samples))

	return all_dists

def get_cond_cdist(dist_type, dist_vars, state):
	'''
		inputs:
			dist_types: (dict) mapping each subsystem (string) to a distribution function,
			dist_vars: (dict) mapping each subsystem (string) to a set of arguments 
				for its distribution function,
			state: (list) of values representing one instance of the joint state
		outputs: 
			(function) representing a conditional continuous distribution for a single subsystem 
				evaluated for the particular value of the joint state
	'''
	return dist_type(*dist_vars(state))

def get_all_cond_cdists(suborder, dist_types, dist_vars, state):
	'''
		inputs: 
			suborder: (dict) mapping each subsystem (string) to an integer, representing the 
				order that will be used for all variables and lists throughout the code,
			dist_types: (dict) mapping each subsystem (string) to a distribution function,
			dist_vars: (dict) mapping each subsystem (string) to a set of arguments 
				for its distribution function,
			state: (list) of values representing one instance of the joint state
		outputs: 
			(list) of conditional continuous distributions for the set of subsystems 
				evaluated for the particular value of the joint state
	'''
	return [get_cond_cdist(dist_types[s], dist_vars[s], state) for s in suborder.keys()]

def get_all_cond_ddists(num_samps, all_bin_edges, suborder, dist_types, dist_vars, state):
	'''
		inputs: 
			num_samps: (int) number of samples to draw, used to construct discrete distributions,
			all_bin_edges: (list) of edges for binning samples for a set of subsystems,
			suborder: (dict) mapping each subsystem (string) to an integer, representing the 
				order that will be used for all variables and lists throughout the code,
			dist_types: (dict) mapping each subsystem (string) to a distribution function,
			dist_vars: (dict) mapping each subsystem (string) to a set of arguments 
				for its distribution function,
			state: (list) of values representing one instance of the joint state
		outputs: 
			all_probs: (list) of conditional discrete distributions for the set of subsystems 
				evaluated for the particular value of the joint state
	'''
	state_dict = {s: state[i] for s, i in suborder.items()}
	all_cond_cdists = get_all_cond_cdists(suborder, dist_types, dist_vars, state_dict)
	all_samples = [dist.rvs(size=num_samps) for dist in all_cond_cdists]
	all_probs = discretize_subsystem_pdfs(all_bin_edges, all_samples)

	return all_probs

def construct_joint_prob_dist(jss, all_ss, all_bin_edges, num_samps, 
								suborder, dist_types, dist_vars, precision):
	'''
		inputs: 
			jss: (list) of lists representing the joint state space for a set of subsystems,
			all_ss: (list) of state spaces for the set of subsystems,
			all_bin_edges: (list) of edges for binning samples for a set of subsystems,
			num_samps: (int) number of samples to draw, used to construct discrete distributions,
			suborder: (dict) mapping each subsystem (string) to an integer, representing the 
				order that will be used for all variables and lists throughout the code,
			dist_types: (dict) mapping each subsystem (string) to a distribution function,
			dist_vars: (dict) mapping each subsystem (string) to a set of arguments 
				for its distribution function,
			precision: (int) float precision with which to record the probabilities
		outputs: 
			jpd: (nparray) representing the probability distribution over the joint state space
	'''
	num_joint_states = len(jss)
	jpd = np.zeros(num_joint_states) # joint probability distribution
	
	# calculate and fill the probability for every state in the jss
	for i in range(num_joint_states):

		jpd[i] = 1 # start with 1 so you can multiply

		# multiply each subsystem's conditional probability contribution to the joint probability
		# this computes, for example, p(x, r1, m, r2) = p(x)*p(r1|x)*p(m|r1)*p(r2|x)
		js = jss[i] #joint state
		sub_cond_dists = get_all_cond_ddists(num_samps, all_bin_edges, 
												suborder, dist_types, dist_vars, js)
		for sub_state, sub_poss_states, s in zip(js, all_ss, range(len(suborder))):
			cp = sub_cond_dists[s][list(sub_poss_states).index(sub_state)]
			jpd[i] = jpd[i]*cp

		jpd[i] = np.round(jpd[i], precision)

	# since the discrete distributions were constructed from sampling, 
	# the sum total might be off in the joint state, so as long as it is within
	# a certain margin of error, renormalize it to 1
	a = sum(jpd)
	if 0.99 < a < 1.01:
		jpd = jpd/a
	else:
		print("Error: sum, %f, of probability distribution is too far from 1", a)
		jpd = jpd/a

	return jpd

def get_marginal_distribution(subset, mjss, jpd, jss, suborder):
	'''
		inputs: 
			subset: (list) of subsystems for which the marginal distribution is desired
			mjss: (list) of lists representing the marginal joint state space (mjss) for the 
				subsystems in subset
			jpd: (nparray) the probability distribution over the entire joint state space
			jss: (list) of lists representing the joint state space for the set of subsystems,
			suborder: (dict) mapping each subsystem (string) to an integer, representing the 
				order that will be used for all variables and lists throughout the code,
		outputs: 
			mpd: (nparray) probability distribution for over the marginal joint state space
	'''
	mpd = np.zeros(len(mjss))
	for i, jms in enumerate(mjss): # one particular joint marginal state (jms)

		# find all states in jss that matches this jms
		for j, js in enumerate(jss):

			# fill in probabilities for joint marginal states by summing over all other states
			match = jms == [js[suborder[sub]] for sub in subset]
			if match:
				mpd[i] += jpd[j]

	return mpd

def get_jpds_over_time(jpd_0, K, ts):
	'''
		inputs: 
			jpd_0: (nparray) initial (t = 0) probability distribution over the joint state space,
			K: (nparray) the rate matrix describing the evolution of the joint distribution,
			ts: (nparray) of times at which the distribution should be calculated
		outputs: 
			jpds: (nparray) of joint distributions at the requested times
	'''
	jpds = np.zeros((len(ts), len(jpd_0)))
	jpds[0] = jpd_0
	for t, time in enumerate(ts):
		if t > 0:
			jpds[t] = np.dot(jpd_0, expm(K*time))
		# print('t', t)

	return jpds

def evolve_whole_system_prob_dists(jpd_0, K, ts, units, uss, jss, suborder, all_ss):
	'''
		inputs: 
			jpd_0: (nparray) initial (t = 0) probability distribution over the joint state space,
			K: (nparray) the rate matrix describing the evolution of the joint distribution,
			ts: (nparray) of times at which the distribution should be calculated,
			units: (list) of units,
			uss: (list) of joint state spaces for the subsystems in each unit in units,
			jss: (list) of lists representing the joint state space for the set of subsystems,
			suborder: (dict) mapping each subsystem (string) to an integer, representing the 
				order that will be used for all variables and lists throughout the code,
			all_ss: (list) of state spaces for each of the subsystems,
		outputs: 
			jpds: (nparray) of joint distributions at the requested times
			all_upds, (list) of nparrays of unit distributions evolved over time, one for each unit
			all_spds, (list) of nparrays of subsystem distributions evolved over time, 
				one for each subsystem

	'''
	# set up the containers for upds and spds
	all_upds = [[get_marginal_distribution(units[u], uss[u], jpd_0, jss, suborder)] 
																	for u in range(len(units))]
	all_spds = [[get_marginal_distribution([sub], all_ss[s], jpd_0, jss, suborder)] 
																	for sub, s in suborder.items()]
	# evolve jpds over time
	jpds = get_jpds_over_time(jpd_0, K, ts)

	# use jpds evolved over time to get the marginal distributions for units and subsystems
	for t, time in enumerate(ts):
		jpds_t = jpds[t]
		if t > 0:
			for u, upds in enumerate(all_upds):
				upds.append(get_marginal_distribution(units[u], uss[u], jpds_t, jss, suborder))
			for s, spds in zip(suborder.items(), all_spds):
				spds.append(get_marginal_distribution([s[0]], all_ss[s[1]], jpds_t, jss, suborder))

	return jpds, all_upds, all_spds

def evolve_unit_prob_dist(upd_0, K_u, ts):
	'''
		inputs: 
			upd_0: (nparray) initial (t = 0) probability distribution over the unit state space,
			K_u, rate matrix to be used to evolve the unit probability distribution,
			ts: (nparray) of times at which the distribution should be calculated
		outputs: 
			upds: (nparray) of unit distributions at the requested times

	'''
	upds = get_jpds_over_time(upd_0, K_u, ts)

	return upds

def comp_probdists(p1, p2, precision):
	'''
		inputs: 
			p1: (nparray) first probability distribution,
			p2: (nparray) second probability distribution,
			precision: (float) desired precision for comparison (to combat floating point errors)
		outputs: 
			isequal: (boolean) representing if p1 == p2 within desired precision
	'''
	isequal = len(p1) == len(p2)
	if isequal == True:
		for a, b in zip(p1, p2):
			isequal = isequal and (a - precision < b <= a + precision)
	return isequal


# ************************ Unsolved Problems With Below ********************************************

def get_all_cond_ddists_directly(suborder, all_ss, dist_types, dist_vars, state):
	'''
		inputs: 
			suborder: (dict) mapping each subsystem (string) to an integer, representing the 
				order that will be used for all variables and lists throughout the code,
			all_ss: (list) of state spaces for the set of subsystems,
			dist_types: (list) mapping each subsystem (string) to a distribution function,
			dist_vars: (list) mapping each subsystem (string) to a set of arguments 
				for its distribution function,
			state: (list) of values representing one instance of the joint state
		outputs: 
			all_probs: (list) of conditional discrete distributions for the set of subsystems 
				evaluated for the particular value of the joint state
	'''
	state_dict = {s: state[i] for s, i in suborder.items()}
	# for dist_type, dist_var in zip(dist_types.items(), dist_vars.items()):
	# 	print(dist_type, dist_vars) 

	# print(state_dict)
	# print([(a[1], b[1](state_dict)) for a, b in zip(dist_types.items(), dist_vars.items())])
	all_cond_ddists = [dist_type[1](*dist_var[1](state_dict)) 
							for dist_type, dist_var in zip(dist_types.items(), dist_vars.items())]

	all_probs = [ddist.pmf(ss) if sum(ddist.pmf(ss)) == 1 else ddist.pmf(ss/dl) 
										for ddist, ss in zip(all_cond_ddists, all_ss)]
	# for probs in all_probs:
	# 	print(probs, sum(probs))
	return all_probs

def construct_joint_prob_dist_directly(jss, all_ss, suborder, dist_types, dist_vars, precision):
	'''
		inputs: 
			jss: (list) of lists representing the joint state space for a set of subsystems,
			all_ss: (list) of state spaces for the set of subsystems,
			suborder: (dict) mapping each subsystem (string) to an integer, representing the 
				order that will be used for all variables and lists throughout the code,
			dist_types: (dict) mapping each subsystem (string) to a distribution function,
			dist_vars: (dict) mapping each subsystem (string) to a set of arguments 
				for its distribution function,
			precision: (int) float precision with which to record the probabilities
		outputs: 
			jpd: (nparray) representing the probability distribution over the joint state space
	'''
	num_joint_states = len(jss)
	jpd = np.zeros(num_joint_states) # joint probability distribution
	
	# calculate and fill the probability for every state in the jss
	for i in range(num_joint_states):

		jpd[i] = 1 # start with 1 so you can multiply

		# multiply each subsystem's conditional probability contribution to the joint probability
		# this computes, for example, p(x, r1, m, r2) = p(x)*p(r1|x)*p(m|r1)*p(r2|x)
		js = jss[i] #joint state
		sub_cond_dists = get_all_cond_ddists_directly(suborder, all_ss, dist_types, dist_vars, js)
		for sub_state, sub_poss_states, s in zip(js, all_ss, range(len(suborder))):
			jpd[i] = jpd[i]*sub_cond_dists[s][list(sub_poss_states).index(sub_state)]

		jpd[i] = np.round(jpd[i], precision)
	print(sum(jpd))

	return jpd
