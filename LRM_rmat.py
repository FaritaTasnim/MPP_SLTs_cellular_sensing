from tools import *
from tools_rmat import *
from LRM_properties import *

# *********************************** K Full System ************************************************

# construct rate matrices for evolving distributions over the entire joint state space
K_LRM = construct_rate_matrix(rmfs_LRM, jss_LRM, ss_steps_LRM, suborder_LRM, 10)
Ks_LRM_X = [get_subsystem_rate_matrix(K_LRM, s, jss_LRM, suborder_LRM) for s in suborder_LRM.keys()]

# ************************************** K Units ***************************************************

# construct rate matrices for evolving distributions over the unit state spaces
Ku_LRM = get_unit_rate_matrices(rmfs_LRM, uss_LRM, units_LRM, ss_steps_LRM, suborder_LRM, 10)
Ks_LRM = [[] for u in units_LRM]
for i, u in enumerate(units_LRM):
	Ks_LRM[i] = [get_subsystem_rate_matrix(Ku_LRM[i], s, uss_LRM[i], 
											{u[i] : i for i in range(len(u))}) for s in u]

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