import numpy as np

def construct_rate_matrix(all_rmfs, jss, ss_steps, suborder, precision):
	'''
		inputs: 
			all_rmfs: (list) of rate matrix functions, which are dictionaries of functions 
				providing the rates at which state transitions happen, for a set of subsystems,
			jss: (list) of lists representing the joint state space for the set of subsystems,
			ss_steps, (list) of steps between adjacent values in the discrete state space 
				for the set of subsystems,
			suborder: (dict) mapping each subsystem (string) to an integer, representing the 
				order that will be used for all variables and lists throughout the code
		outputs: 
			K: (nparray) representing the rate matrix for the set of subsystems,
				K[i][j] is the rate of transitions between joint state i and joint state j
	'''
	num_joint_states = len(jss)
	K = np.zeros((num_joint_states, num_joint_states))

	# fill in all off-diagonal elements
	for i in range(num_joint_states):

		before = list(np.copy(jss[i])) # get state before transition

		# construct all states reachable after a single transition
		# since only one subsystem can change in a single transition, loop through subsystems 
		for sub, rmfs, step in zip(range(len(before)), all_rmfs, ss_steps): 

			# construct the state after a single transition
			for mult in rmfs.keys():
				after = list(np.copy(jss[i])) 
				after[sub] = before[sub] + mult*step 

				# if the state after the transition is a valid state
				if after in jss: 
					# calculate the rate from the rate matrix function
					rate = rmfs[mult]({s : before[suborder[s]] for s in suborder.keys()})
					# and set the appropriate element in the rate matrix
					K[i][jss.index(after)] = np.round(rate, precision)

					#sanity check: this should print nothing otherwise
					#you are violating the multipartite assumption!
					# check = np.array([a-b for a,b in zip(after, before)])
					# if np.count_nonzero(check) > 1:
					# 	print (after, before)
		
	#all diagonal elements should be filled last to ensure sum_x K[x'][x] = 0	
	for i, row in enumerate(K):
		K[i][i] = -1*sum(row)

	return K

def get_subsystem_rate_matrix(K, sub, jss, suborder):
	'''
		inputs: 
			K: (nparray) full joint state rate matrix,
			sub: (string) subsystem for which rate matrix K(s) is desired
			jss: (list) of lists representing the joint state space for the set of subsystems,
			suborder: (dict) mapping each subsystem (string) to an integer, representing the 
				order that will be used for all variables and lists throughout the code
		outputs: 
			K_s: (nparray) rate matrix for sub
	'''
	
	K_s = np.zeros(K.shape) # K(s) must match the shape of K
	num_joint_states = len(jss)
	s = suborder[sub] # index of sub in all lists
	not_s = [suborder[i] for i in suborder if i != sub] # list of indices of not(sub)

	# loop through all entries
	for i in range(num_joint_states):
		for j in range(num_joint_states):

			# only set elements for which ONLY sub changes state, & other subsystem states are fixed
			before, after = jss[i], jss[j]
			fill = (before[s] != after[s])
			for ns in not_s:
				fill = fill and (before[ns] == after[ns])
			if fill == True:
				K_s[i][j] = K[i][j]

	#update the diagonal entries to ensure sum_x K_s[x'][x] = 0
	for i, row in enumerate(K_s):
		K_s[i][i] = -1*sum(row)

	return K_s

def get_unit_rate_matrices(all_rmfs, uss, units, ss_steps, suborder, precision):
	'''
		inputs: 
			all_rmfs: (list) of rate matrix functions, which are dictionaries of functions 
				providing the rates at which state transitions happen, for a set of subsystems,
			uss: (list) of joint state spaces for the subsystems in each unit in units,
			units: (list) of units,
			ss_steps: (list) of statestepval for each subsystem in the entire system,
			suborder: (dict) mapping each subsystem (string) to an integer, representing the 
				order that will be used for all variables and lists throughout the code
		outputs: 
			urmats: (list) of rate matrices for each unit in units, sized correctly to match each
				unit's state space in uss
	'''
	urmats = [construct_rate_matrix([all_rmfs[suborder[a]] for a in u], 
									uss[units.index(u)], 
									[ss_steps[suborder[s]] for s in u], 
									{u[i] : i for i in range(len(u))}, precision) for u in units]
	return urmats

def get_subset_rate_matrix(subset, all_Ks, K, suborder):
	'''
		inputs: 
			subset: (list) of subsystems for which the combined rate matrix is desired,
			all_Ks: (list) of all subsystem rate matrices K_s, each of which is an nparray
				matching the shape of K,
			K: (nparray) system or unit rate matrix,
			suborder: (dict) mapping each subsystem (string) to an integer, representing the 
				order that will be used for all variables and lists throughout the code
		outputs: 
			K_set: (np array) rate matrix for subset
	'''
	K_set = np.zeros(K.shape)
	subs = [suborder[i] for i in suborder if i in subset]

	for sub in subs:
		K_set = K_set + all_Ks[sub]

	return K_set

def comp_rmats(K1, K2, precision):
	'''
		inputs: 
			K1: (nparray) first rate matrix,
			K2: (nparray) second rate matrix,
			precision: (float) desired precision of similarity (to combat floating point errors)
		outputs: 
			isequal: (boolean) representing if K1 == K2 within the desired precision
	'''
	isequal = K1.shape == K2.shape
	if isequal == True:
		for i in range(len(K1[0])):
			for j in range(len(K1[0])):
				isequal = isequal and (K1[i][j] - precision < K2[i][j] <= K1[i][j] + precision)
	return isequal