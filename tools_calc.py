import numpy as np
from scipy.integrate import simps
from scipy.linalg import expm

def L_t(pds_t):
	'''
		inputs: 
			pds_t: (nparray) of probability distributions at each timestep
		outputs: 
			Ls_t: (nparray) of total variation distance from the initial distribution to the 
				distribution at each timestep
	'''
	pd_0 = pds_t[0] # initial probability distribution
	Ls_t = []

	for pd in pds_t: # at each time point

		# compute the total variation distance between the distribution at this time point
		# and the initial distribution (this is a measure of how much the distribution has changed)
		# L(t) = sum_x |p_x(0) - p_x(t)|
		L = 0
		for i, prob in enumerate(pd_0): 
			L += abs(prob - pd[i])
		Ls_t.append(L)

	return np.array(Ls_t)

def Ldot_t(pds_t, K, dt):
	'''
		inputs: 
			pds_t: (nparray) of probability distributions at each timestep,
			K: (nparray) the (potentially counterfactual) rate matrix describing instantaneous 
				evolution of the distribution,
			dt: (float) infinitesimal time over which to perform the derivative approximation
		outputs: 
			Ldots_t: (nparray) of total variation distance from the initial distribution to the 
				distribution at each timestep
	'''

	Ldots_t = []

	for pd in pds_t: # at each time point

		# compute the instantaneous change (during time dt) in the distribution
		Ldot = 0
		pd_tdt = np.dot(pd, expm(K*dt))
		for i, prob in enumerate(pd): 
			# do abs on each term or once after total sum?
			Ldot += abs(prob - pd_tdt[i])/dt
		Ldots_t.append(Ldot)

	return np.array(Ldots_t)

def A_t(pds_t, K):
	'''
		inputs: 
			pds_t: (nparray) of probability distributions at each timestep,
			K: (nparray) the rate matrix describing the evolution of the distribution
		outputs: l
			As_t: (nparray) of the dynamical activity in the state space at each timestep
	'''
	As_t = []

	for pd in pds_t: # at each time point

		# A(t) = sum_{x' neq x} K^{x'}_{x} * p_{x'}(t)
		# remember Dot(p_x(t)) = sum_x' K^{x'}_{x} * p_{x'}(t)
		# so A(t) = sum([(Dot(p_x(t)) - K^{x}_{x} * p_{x}(t)) for x in possible_states])
		# the following performs this calculation

		A = np.dot(pd, K)
		for state, prob in enumerate(pd):
			A[state] = A[state] - K[state][state]*prob

		# collapse into a sum to get dynamical activity
		As_t.append(sum(A))

	# note: the following code does the same calculation, but is more time costly
	# as_t = []
	# num_states = K.shape[0]
	# for pd in pds_t:
	# 	a = 0
	# 	for before in range(num_states):
	# 		for after in range(num_states):
	# 			if before != after:
	# 				a += K[before][after] * pd[before]
	# 	as_t.append(a)

	# return np.array(as_t)

	# print('dyn act old', As_t)
	# print('dyn act new', as_t)
	# print("")

	return np.array(As_t)

def EPrate_t(pds_t, K):
	'''
		inputs: 
			pds_t: (nparray) of probability distributions at each timestep,
			K: (nparray) the rate matrix describing the evolution of the distribution
		outputs: 
			EPrates: (nparray) of the entropy production rate in the universe as brought about 
				by the state space represented by K and pds_t at each timestep
	''' 
	num_states = K.shape[0]
	EPrates = []

	for pd in pds_t: # at each time point

		# langle Dot{sigma}(t) rangle = 
		#	sum_{x', x} K^{x'}_{x} * p_{x'}(t) * ln[frac{K^{x'}_{x}*p_{x'}(t)}{K^{x}_{x'}*p_{x}(t)}]
		EPrate = 0
		for i in range(num_states): # x' state
			for j in range(num_states): # x state
				a = K[i][j]*pd[i]
				b = K[j][i]*pd[j]
				# print(a, b)
				if a != 0 and b != 0: # this filters for transitions that are actually possible
					EPrate += a*np.log(a/b)
		EPrates.append(EPrate)

	return np.array(EPrates)

def S_t(pds_t):
	'''
		inputs:
			pds_t: (nparray) of probability distributions at each timestep,
		outputs:
			Ss: (nparray) of the increase in entropy in the system as brought about by the 
				state space represented by K and pds_t at at each timestep
	'''
	
	num_states = len(pds_t[0])
	Ss = []

	for pd in pds_t: # at each time point

		# S = sum_{x} -p_{x}(t) * ln(p_{x}(t))
		S = 0
		for prob in pd: # x state
			if prob > 0:
				S += -1*prob*np.log(prob)
		Ss.append(S)

	return np.array(Ss)

def get_L_A_EPrate_S_ts(pds_t, K):
	'''
		inputs: 
			pds_t: (nparray) of probability distributions at each timestep,
			K: (nparray) the rate matrix describing the evolution of the distribution
		outputs: 
			Ls_t: (nparray) of total variation distance from the initial distribution to the 
				distribution at each timestep,
			As: (nparray) of the dynamical activity in the state space at each timestep,
			EPrates: (nparray) of the entropy production rate in the universe as brought about 
				by the state space represented by K and pds_t at each timestep,
			Ss: (nparray) of the increase in entropy in the system as brought about by the 
				state space represented by K and pds_t at at each timestep
	''' 
	num_states = len(pds_t[0])
	pd_0 = pds_t[0] # initial probability distribution
	Ls, As, EPrates, Ss = [], [], [], []

	for pd in pds_t: # at each time point

		L, A, EPrate, S = 0, np.dot(pd, K), 0, 0
		for i in range(num_states): # x' state

			before_prob = pd[i]

			L += abs(pd_0[i] - before_prob)

			A[i] = A[i] - K[i][i]*before_prob

			if before_prob > 0:
				S += -1*before_prob*np.log(before_prob)

			for j in range(num_states): # x state
				after_prob = pd[j]

				a = K[i][j]*before_prob
				b = K[j][i]*after_prob
				# print(a, b)
				if a != 0 and b != 0: # this filters for transitions that are actually possible
					EPrate += a*np.log(a/b)
		
		Ls.append(L)
		As.append(sum(A))
		EPrates.append(EPrate)
		Ss.append(S)

	return np.array([np.array(Ls), np.array(As), np.array(EPrates), np.array(Ss)])

def DotS_t(pds_t, K):
	'''
		inputs:
			pds_t: (nparray) of probability distributions at each timestep,
			K: (nparray) the rate matrix describing the evolution of the distribution
		outputs:
			DotSs: (nparray) of the increase in entropy in the system as brought about by the 
				state space represented by K and pds_t at at each timestep
	'''
	
	num_states = K.shape[0]
	DotSs = []

	for pd in pds_t: # at each time point

		# Dot_S = sum_{x', x} -K^{x'}_{x} * p_{x'}(t) * ln(p_{x}(t))
		DotS = 0
		for i in range(num_states): # x' state
			for j in range(num_states): # x state
				a = pd[j]
				if a > 0:
					DotS += -1*K[i][j]*pd[i]*np.log(a)
		DotSs.append(S)

	return np.array(DotSs)

def group_vals_t(group, unit_vals_t):
	'''
		inputs:
			group: (list) of indices of units that are in the group
			unit_vals: (nparray) of values of each unit at each timestep
		outputs:
			group_sums: (nparray) of the sum of values of the units in group at each timestep
	'''
	group_sums = []

	for unit_vals in unit_vals_t.T: # at each time point, transpose to get [time][unit] order
		group_sum = 0

		for unit_ind in group:
			group_sum += unit_vals[unit_ind]

		group_sums.append(group_sum)

	return np.array(group_sums)

def InExInfo_t(Ss, unit_Ss_t, odd_evens):
	'''
		inputs:
			Ss: (nparray) of S values for the system at each timestep
			unit_Ss_t: (list) of nparrays of S values for each unit at each timestep
			odd_evens: (list) of two lists of indices of odd units and indices of even units
		outputs:
			InExInfos: (nparray) of the InExInfo at each timestep
	'''
	odd_units, even_units = odd_evens[0], odd_evens[1]
	odd_Ss, even_Ss = group_vals_t(odd_units, unit_Ss_t), group_vals_t(even_units, unit_Ss_t)
	InExInfos = [odd_S - even_S - S for odd_S, even_S, S in zip(odd_Ss, even_Ss, Ss)]
	
	return np.array(InExInfos)

def InExf_t(unit_fs_t, odd_evens):
	'''
		inputs:
			unit_fs_t: (list) of nparrays of values of some functional f for each unit 
				at each timestep 
			odd_evens: (list) of two lists of indices of odd units and indices of even units)
		outputs:
			InExfs: (nparray) of the InEx sum of the functional f at each timestep
	'''
	odd_units, even_units = odd_evens[0], odd_evens[1]
	odd_fs, even_fs = group_vals_t(odd_units, unit_fs_t), group_vals_t(even_units, unit_fs_t)
	InExfs = [odd_f - even_f for odd_f, even_f in zip(odd_fs, even_fs)]
	
	return np.array(InExfs)

def Deltaf_t(fs_t):
	'''
		inputs:
			fs_t: (nparray) of f(t) over time
		outputs:
			Deltafs: (nparray) of the Delta of f over time, i.e. f(t) - f(0) at each timestep
	'''
	f_0 = fs_t[0]
	Deltafs = [f - f_0 for f in fs_t]

	return np.array(Deltafs)

def integrate(times, quants):
	'''
		inputs:
			times: (nparray) of times over which to integrate,
			quants: (list) of quantities at each of those times, this will be integrated
		outputs:
			integ_results: (nparray) of integration results by each time step
	'''
	integ_results = []
	for t in range(len(times)):
		if t == 0:
			integ_results.append(0.0)
		else: 
			integ_results.append(simps(quants[:t+1], times[:t+1]))

	return np.array(integ_results)

def time_avg(times, quants):
	'''
		inputs:
			times: (nparray) of times over which to average,
			quants: (list) of quantities at each of those times
		outputs:
			avg: (nparray) of time averages of quants by each time step
	'''
	avg = []
	vals = integrate(times, quants)
	for t, val in enumerate(vals):
		if t == 0:
			avg.append(quants[t])
		else:
			avg.append(val/(times[t]))

	return np.array(avg)

def classical_EP_bound(L, A, tau):
	return (L**2)/(2*A*tau)

def classical_tau_bound(L, A, EP):
	return (L**2)/(2*A*EP)

def classical_bound(L, A, EP_or_t):
	return (L**2)/(2*A*EP_or_t)

def mpp_tau_bound(Lsyst, Ais, EPis):
	# print (Lsyst, Ais, EPis)
	denom = 2*(sum([np.sqrt(a*e) for a,e in zip(Ais, EPis)]))**2
	return Lsyst**2/denom

def MPP_RegSum_bound(Ls, As, x):
	'''
		inputs:
			Ls: (list) of total variation distances at time tau for each of the terms in the bound,
			As: (list) of dynamical activities at time tau for each of the terms in the bound,
			tau: (float) time for which EP bound is being calculated
		outputs:
			bound: (float) bound in MPP at a given time, Eq. (),
	'''
	bound = 0
	for l, a in zip(Ls, As):
		bound += classical_bound(l, a, x)
	return bound

def MPP_RegSum_contribs(Ls, As, x):
	'''
		inputs:
			Ls: (list) of total variation distances at time tau for each of the terms in the bound,
			As: (list) of dynamical activities at time tau for each of the terms in the bound,
			tau: (float) time for which EP bound is being calculated
		outputs:
			contribs: (nparray) of contributions to the bound from each separate term
	'''
	contribs = []
	for l, a in zip(Ls, As):
		contribs.append(classical_bound(l, a, x))
	return np.array(contribs)

def MPP_InExSum_bound(Lus, Aus, EPs_or_ts, odd_evens):
	'''
		inputs:
			Lus: (nparray) of Ls for each unit for each timestep,
			Aus: (nparray) of As for each unit for each timestep,
			EPs_or_ts: (nparray) of system EPs or time at each timestep,
			odd_evens: (list) of two lists of indices of odd units and indices of even units)
		outputs:
			InEx_bounds: (nparray) of in-ex sum bounds at each timestep, Eq. ()
	'''
	InEx_bounds = []
	odds, evens = odd_evens[0], odd_evens[1]

	# transpose to get [time][unit value] order instead of [unit][time value] order
	for Ls, As, EP_or_t in zip(Lus.T, Aus.T, EPs_or_ts):
		if EP_or_t > 0: # ignore the trivial t=0 bound
			odd_vals = [classical_bound(Ls[i], As[i], EP_or_t) for i in odds]
			even_vals = [classical_bound(Ls[i], As[i], EP_or_t) for i in evens]
			InEx_bounds.append(sum(odd_vals) - sum(even_vals))

	return np.array(InEx_bounds)

def MPP_ReArr_InEx_EP_bound(Lus, Aus, ts, EPus, DeltaInExInfo, odd_evens):
	'''
		inputs:
			Lus: (nparray) of Ls for each unit for each timestep,
			Aus: (nparray) of As for each unit for each timestep,
			ts: (nparray) of time at each timestep,
			EPus: (nparray) of EPs for each unit for each timestep,
			DeltaInExInfo: (nparray) of the change in in-ex info each timestep (i(t) - i(0)),
			odd_evens: (list) of two lists of indices of odd units and indices of even units)
		outputs:
			ReArr_InEx_bounds: (nparray) of rearranged in-ex sum EP bounds at each timestep, Eq. ()
	'''
	ReArr_InEx_bounds = []
	odds, evens = odd_evens[0], odd_evens[1]

	# transpose to get [time][unit value] order instead of [unit][time value] order
	for Ls, As, EPs, dinexinfo, t in zip(Lus.T, Aus.T, EPus.T, DeltaInExInfo, ts):
		if t > 0: # ignore the trivial t=0 bound
			odd_vals = [classical_bound(Ls[i], As[i], t) for i in odds]
			even_vals = [EPs[i] for i in evens]
			ReArr_InEx_bounds.append(sum(odd_vals) - sum(even_vals) - dinexinfo)

	return np.array(ReArr_InEx_bounds)

def MPP_ReArr_InEx_tau_bound(Lus, Aus, ts, EPus, EPsys_t, DeltaInExInfo, odd_evens):
	'''
		inputs:
			Lus: (nparray) of Ls for each unit for each timestep,
			Aus: (nparray) of As for each unit for each timestep,
			ts: (nparray) of time at each timestep,
			EPus: (nparray) of EPs for each unit for each timestep,
			EPsys_t: (nparray) of system EP at each timestep,
			DeltaInExInfo: (nparray) of the change in in-ex info each timestep (i(t) - i(0)),
			odd_evens: (list) of two lists of indices of odd units and indices of even units)
		outputs:
			ReArr_InEx_bounds: (nparray) of rearranged in-ex sum tau bounds at each timestep, Eq. ()
	'''
	ReArr_InEx_bounds = []
	odds, evens = odd_evens[0], odd_evens[1]

	# transpose to get [time][unit value] order instead of [unit][time value] order
	for Ls, As, EPs, epsys, dinexinfo, t in zip(Lus.T, Aus.T, EPus.T, EPsys_t, DeltaInExInfo, ts):
		if t > 0: # ignore the trivial t=0 bound
			weird_EP_sum = sum([EPs[i] for i in evens]) + epsys + dinexinfo
			ReArr_InEx_bounds.append(sum([classical_bound(Ls[i], As[i], weird_EP_sum) 
																				for i in odds]))

	return np.array(ReArr_InEx_bounds)