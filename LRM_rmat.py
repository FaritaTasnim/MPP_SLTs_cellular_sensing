from tools import *
from tools_rmat import *
from LRM_properties import *

# *********************************** K Full System ************************************************

# construct rate matrices for evolving distributions over the entire joint state space
K_LRM = construct_rate_matrix(rmfs_LRM, jss_LRM, ss_steps_LRM, suborder_LRM, 16)
Ks_LRM_X = [get_subsystem_rate_matrix(K_LRM, s, jss_LRM, suborder_LRM) for s in suborder_LRM.keys()]
# Ku_LRM_X = [get_subset_rate_matrix(unit, Ks_LRM_X, K_LRM, suborder_LRM) for unit in units_LRM]

# ************************************** K Units ***************************************************

# construct rate matrices for evolving distributions over the unit state spaces
Ku_LRM = get_unit_rate_matrices(rmfs_LRM, uss_LRM, units_LRM, ss_steps_LRM, suborder_LRM, 10)
Ks_LRM = [[] for u in units_LRM]
for i, u in enumerate(units_LRM):
	Ks_LRM[i] = [get_subsystem_rate_matrix(Ku_LRM[i], s, uss_LRM[i],
											{u[i] : i for i in range(len(u))}) for s in u]

# # ***************************** K Units, Unit Structure 2 ******************************************

# # construct rate matrices for evolving distributions over the unit state spaces
# Ku_LRM2 = get_unit_rate_matrices(rmfs_LRM, uss_LRM2, units_LRM2, ss_steps_LRM, suborder_LRM, 10)
# Ks_LRM2 = [[] for u in units_LRM2]
# for i, u in enumerate(units_LRM2):
# 	Ks_LRM2[i] = [get_subsystem_rate_matrix(Ku_LRM2[i], s, uss_LRM2[i],
# 											{u[i] : i for i in range(len(u))}) for s in u]

# # ************************ K Counterfactuals, Unit Structure 1 *************************************

# # subsets of subsystems needed for counterfactual rate matrices over the entire joint state space
# N_phi = [i for i in suborder_LRM.keys() if i not in units_LRM[uind['phi']]]
# N_omega = [i for i in suborder_LRM.keys() if i not in units_LRM[uind['omega']]]
# N_alpha = [i for i in suborder_LRM.keys() if i not in units_LRM[uind['alpha']]]

# # counterfactual rate matrix K_N_A, i.e. rate matrix for the entire system if only the subsystems
# # in N_A (= N\A) are allowed to change state
# K_N_phi = get_subset_rate_matrix(N_phi, Ks_LRM_X, K_LRM, suborder_LRM)
# K_N_omega = get_subset_rate_matrix(N_omega, Ks_LRM_X, K_LRM, suborder_LRM)
# K_N_alpha = get_subset_rate_matrix(N_alpha, Ks_LRM_X, K_LRM, suborder_LRM)

# # subsets of subsystems needed for counterfactual rate matrices over a unit state space
# omega, alpha = units_LRM[uind['omega']], units_LRM[uind['alpha']]
# omega_phi = [i for i in units_LRM[uind['omega']] if i not in units_LRM[uind['phi']]]
# alpha_phi = [i for i in units_LRM[uind['alpha']] if i not in units_LRM[uind['phi']]]

# # counterfactual rate matrix K_W_W', i.e. rate matrix for the unit W if only the subsystems
# # in W_W' (= W\W') are allowed to change state
# K_omega_phi = get_subset_rate_matrix(omega_phi, Ks_LRM[uind['omega']],
# 								Ku_LRM[uind['omega']], {omega[i] : i for i in range(len(omega))})
# K_alpha_phi = get_subset_rate_matrix(alpha_phi, Ks_LRM[uind['alpha']],
# 								Ku_LRM[uind['alpha']], {alpha[i] : i for i in range(len(alpha))})

# # ************************ K Counterfactuals, Unit Structure 1 *************************************

# # subsets of subsystems needed for counterfactual rate matrices over the entire joint state space
# N_phi = [i for i in suborder_LRM.keys() if i not in units_LRM2[uind2['phi']]]
# N_omega = [i for i in suborder_LRM.keys() if i not in units_LRM2[uind2['omega']]]
# N_alpha = [i for i in suborder_LRM.keys() if i not in units_LRM2[uind2['alpha']]]
# N_beta = [i for i in suborder_LRM.keys() if i not in units_LRM2[uind2['beta']]]

# # counterfactual rate matrix K_N_A, i.e. rate matrix for the entire system if only the subsystems
# # in N_A (= N\A) are allowed to change state
# K_N_phi = get_subset_rate_matrix(N_phi, Ks_LRM_X, K_LRM, suborder_LRM)
# K_N_omega = get_subset_rate_matrix(N_omega, Ks_LRM_X, K_LRM, suborder_LRM)
# K_N_alpha = get_subset_rate_matrix(N_alpha, Ks_LRM_X, K_LRM, suborder_LRM) # to use
# K_N_beta = get_subset_rate_matrix(N_beta, Ks_LRM_X, K_LRM, suborder_LRM) # to use

# # subsets of subsystems needed for counterfactual rate matrices over a unit state space
# omega = units_LRM2[uind2['omega']]
# alpha = units_LRM2[uind2['alpha']]
# beta = units_LRM2[uind2['beta']]

# omega_phi = [i for i in units_LRM2[uind2['omega']] if i not in units_LRM2[uind2['phi']]]
# omega_alpha = [i for i in units_LRM2[uind2['omega']] if i not in units_LRM2[uind2['alpha']]]
# omega_beta = [i for i in units_LRM2[uind2['omega']] if i not in units_LRM2[uind2['beta']]]

# alpha_phi = [i for i in units_LRM2[uind2['alpha']] if i not in units_LRM2[uind2['phi']]]
# beta_phi = [i for i in units_LRM2[uind2['beta']] if i not in units_LRM2[uind2['phi']]]

# # counterfactual rate matrix K_W_W', i.e. rate matrix for the unit W if only the subsystems
# # in W_W' (= W\W') are allowed to change state
# K_omega_phi = get_subset_rate_matrix(omega_phi, Ks_LRM2[uind2['omega']],
# 						Ku_LRM2[uind2['omega']], {omega[i] : i for i in range(len(omega))})
# K_omega_alpha = get_subset_rate_matrix(omega_alpha, Ks_LRM2[uind2['omega']],
# 						Ku_LRM2[uind2['omega']], {omega[i] : i for i in range(len(omega))})
# K_omega_beta = get_subset_rate_matrix(omega_beta, Ks_LRM2[uind2['omega']], # to use
# 						Ku_LRM2[uind2['omega']], {omega[i] : i for i in range(len(omega))})

# K_alpha_phi = get_subset_rate_matrix(alpha_phi, Ks_LRM2[uind2['alpha']],
# 						Ku_LRM2[uind2['alpha']], {alpha[i] : i for i in range(len(alpha))})
# K_beta_phi = get_subset_rate_matrix(beta_phi, Ks_LRM2[uind2['beta']], # to use
# 						Ku_LRM2[uind2['beta']], {beta[i] : i for i in range(len(beta))})

print('Constructed all rate matrices')

# *********************************** Sanity Checks ************************************************

# print(np.argwhere(K_LRM < 0)) # should be only the diagonal terms

# print([sum(i) for i in Ku_LRM[1]])
# print(Ku_LRM[1])

# print([np.count_nonzero(k) for k in Ks_LRM_X])

# print('K constructed')
# print(np.count_nonzero(K_LRM))

# print('Ku constructed')
# print([np.count_nonzero(Ku) for Ku in Ku_LRM])
# print([Ku.shape for Ku in Ku_LRM])

# for Ks in Ks_LRM:
# 	print([K.shape for K in Ks])

# K_omega = Ks_LRM[1][0] + Ks_LRM[1][1] + Ks_LRM[1][2]
# print('isequal', comp_rmats(Ku_LRM[1], K_omega, 0.00001))

# K_sum = K_l + K_r1 + K_m + K_r2
# print('isequal', comp_rmats(K_LRM, K_sum, 0.00001))
# print('isequal', comp_rmats(K_omega, K_l + K_r1 + K_m, 0.00001))

# print(N_phi, N_omega, N_alpha)
# print(omega_phi, alpha_phi)
