from tools import *
from tools_pdist import *
from LRM_properties import *
from LRM_rmat import *
# import time

# **************************** Initial Probability Distribution ************************************

# set the parameters for the initial conditional distributions
# sets the distribution types for each of the subsystem initial conditional distributions
dist_types_LRM = {"l": norm,"r1": norm, "m": norm, "r2": norm}
# takes a joint state and constructs the parameters for the initial conditional distributions
dist_vars_LRM = {"l": lambda js: [0, 0.5],
                "r1": lambda js: [exp(js["l"]-lmax)*Nr1, 2],
                "m": lambda js: [exp(js["l"]-lmax)*Nm, 3],
                "r2": lambda js: [exp(js["l"]-lmax)*Nr2, 3]}

# constructs the t=0 joint probability distribution given the subsystem conditional distributions
jpd_LRM_0 = construct_joint_prob_dist(jss_LRM, all_ss_LRM, bin_edges_LRM, 10000,
													suborder_LRM, dist_types_LRM, dist_vars_LRM, 10)

print('Constructed initial probability distributions')

# ************************** Probability Distributions, All Times **********************************

# generate joint, unit, and subsystem probability distributions under actual rate matrix
jpds_LRM, all_upds_LRM, all_spds_LRM = evolve_whole_system_prob_dists(jpd_LRM_0, K_LRM, times_LRM,
											units_LRM, uss_LRM, jss_LRM, suborder_LRM, all_ss_LRM)

print('Evolved probability distributions over time')


# *********************************** Sanity Checks ************************************************

# later, figure out which is faster, getting the unit marginals from evolving the system prob dist,
# or evolving the marginal unit prob dists directly themselves

# come back to trying binom and uniform later

# dist_types_LRM = {"l": randint,"r1": binom, "m": binom, "r2": binom}
# # takes a joint state and constructs the parameters for the initial conditional distributions
# dist_vars_LRM = {"l": lambda js: [int(ss_mins_LRM[0]/dl), int(ss_maxes_LRM[0]/dl)],
# 				"r1": lambda js: [Nr1, js["l"]],
# 				"m": lambda js: [Nm, js["r1"]/Nr1],
# 				"r2": lambda js: [Nr2, js["l"]]}

# jpd_LRM_0 = construct_joint_prob_dist_directly(jss_LRM, all_ss_LRM, suborder_LRM,
# 													dist_types_LRM, dist_vars_LRM, 10)


# for spd in all_spds_LRM:
# 	print(spd[0]) # initial prob dist for each subsystem
# 	print(spd[-1]) # final prob dist for each subsystem
# 	print("")

# print(jpd_LRM_0)
# print(sum(jpd_LRM_0))

# for t, time in enumerate(times_LRM):
# 	jpds_LRM_t = jpds_LRM[t]
# 	if t > 0:
# 		for u, upd in enumerate(upds_LRM):
# 			upd.append(get_marginal_distribution(units_LRM[u], jpds_LRM_t, jss_LRM, all_ss_LRM, suborder_LRM))
# 		for s, spd in zip(suborder_LRM.keys(), spds_LRM):
# 			spd.append(get_marginal_distribution([s], jpds_LRM_t, jss_LRM, all_ss_LRM, suborder_LRM))

# 	plt.plot(jpds_LRM_t)
# 	print(sum(jpds_LRM_t))
# plt.show()

# for u in range(len(units_LRM)):
# 	upd = upds_LRM[u]
# 	for t, time in enumerate(times_LRM):
# 		upd_t = upd[t]
# 		plt.plot(upd_t)
# 		print(sum(upd_t))
# 	plt.show()

# color = plt.cm.rainbow(np.linspace(0,1,len(times_LRM)))

# ignore this, write another function for plotting the distributions
# for s in range(len(suborder_LRM.keys())):
# 	spd = spds_LRM[s]
# 	for t, time in enumerate(times_LRM):
# 		spd_t = spd[t]
# 		# print(all_ss_LRM[s])
# 		# print(spd_t)
# 		plt.vlines(all_ss_LRM[s], 0, spd_t, colors = color[t], lw = 20)
# 		# print(sum(spd_t))
# 	plt.show()

# sanity check if CTMC master equation over unit omega holds

# fuck to the yessss it works
# test = get_jpds_over_time(upds_LRM[1][0], Ku_LRM[1], times_LRM)
# print(comp_probdists(upds_LRM[1][1], test[1], 0.00001))


# tries to come back to

# def get_all_cond_ddists(dist_types, dist_vars, all_ss, state):
# 	return [dist_type(*dist_var(state)).pmf(ss)
# 				for dist_type, dist_var, ss in zip(dist_types, dist_vars, all_ss)]

# def construct_joint_prob_dist():
# 	pass

# dist_types_LRM2 = {"l": uniform,"r1": binom, "m": binom, "r2": binom}
# # takes a joint state and constructs the parameters for the initial conditional distributions
# dist_vars_LRM2 = {"l": lambda js: [int(ss_mins_LRM[0]/dl), int(ss_maxes_LRM[0]/dl)],
# 				"r1": lambda js: [Nr1, js["l"]],
# 				"m": lambda js: [Nm, js["r1"]/Nr1],
# 				"r2": lambda js: [Nr2, js["l"]]}


# print(get_all_cond_ddists(dist_types_LRM2, dist_vars_LRM2, all_ss_LRM, {"l": 0.6, "r1": 2, "m": 3, "r2": 2}))

# # ********* Counterfactual Probability Distributions, All Times, Unit Structure 1 ******************

# # generate joint, unit, and subsystem probability distributions under counterfactual rate matrices
# jpds_N_omega, all_upds_N_omega, all_spds_N_omega = evolve_whole_system_prob_dists(jpd_LRM_0,
# 						K_N_omega, times_LRM, units_LRM, uss_LRM, jss_LRM, suborder_LRM, all_ss_LRM)

# jpds_N_alpha, all_upds_N_alpha, all_spds_N_alpha = evolve_whole_system_prob_dists(jpd_LRM_0,
# 						K_N_alpha, times_LRM, units_LRM, uss_LRM, jss_LRM, suborder_LRM, all_ss_LRM)

# omegapd_0 = all_upds_LRM[uind['omega']][0]
# uss_omega = uss_LRM[uind['omega']]
# omegapds_omega_phi = evolve_unit_prob_dist(omegapd_0, K_omega_phi, times_LRM)

# alphapd_0 = all_upds_LRM[uind['alpha']][0]
# uss_omega = uss_LRM[uind['alpha']]
# alphapds_alpha_phi = evolve_unit_prob_dist(alphapd_0, K_alpha_phi, times_LRM)

# print('Evolved counterfactual probability distributions over time, unit structure 1')

# # ********* Counterfactual Probability Distributions, All Times, Unit Structure 2 ******************

# # only evolve the non-repetitive ones because this is time-costly
# # generate joint, unit, and subsystem probability distributions under counterfactual rate matrices
# jpds_N_omega, all_upds_N_omega, all_spds_N_omega = jpds_N_omega, all_upds_N_omega, all_spds_N_omega

# jpds_N_beta, all_upds_N_beta, all_spds_N_beta = evolve_whole_system_prob_dists(jpd_LRM_0,
# 					K_N_beta, times_LRM, units_LRM2, uss_LRM2, jss_LRM, suborder_LRM, all_ss_LRM)

# omegapd_0 = all_upds_LRM2[uind2['omega']][0]
# uss_omega = uss_LRM2[uind2['omega']]
# omegapds_omega_beta = evolve_unit_prob_dist(omegapd_0, K_omega_beta, times_LRM)

# betapd_0 = all_upds_LRM2[uind2['beta']][0]
# uss_beta = uss_LRM2[uind2['beta']]
# betapds_beta_phi = evolve_unit_prob_dist(betapd_0, K_beta_phi, times_LRM)

# print('Evolved counterfactual probability distributions over time, unit structure 2')
