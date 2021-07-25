from tools_calc import *
from tools_unit import *
from LRM_properties import *
from LRM_rmat import *
from LRM_pdist import *
# from matplotlib import pyplot as plt

# ************************ Thermodynamic Quantities for System *************************************

Lsys, Asyst, EPratesys, Ssys = get_L_A_EPrate_S_ts(jpds_LRM, K_LRM)
Asys = time_avg(times_LRM, Asyst)
EPsys = integrate(times_LRM, EPratesys)


# ********************** Thermodynamic Quantities for Subsystems ***********************************

# ss_iter, the iterable over all subsystems, was defined in LRM_properties
sub_calcs = np.array([get_L_A_EPrate_S_ts(jpds_LRM, Ks_LRM_X[s]) for s in ss_iter])
Lx, Ait, EPratei, Si = sub_calcs[:,0,:], sub_calcs[:,1,:], sub_calcs[:,2,:], sub_calcs[:,3,:]
Ai = np.array([time_avg(times_LRM, Ait[s]) for s in ss_iter])
EPi = np.array([integrate(times_LRM, epratei) for epratei in EPratei])
Li = np.array([L_t(all_spds_LRM[s]) for s in ss_iter])


#---quick, very important sanity check to ensure that the subsystem contributions to
#------global and local quantities are the same!

# ss_itero = range(len(units_LRM[uind['omega']]))
# sub_calcso = np.array([get_L_A_EPrate_S_ts(all_upds_LRM[uind['omega']], Ks_LRM[uind['omega']][s]) for s in ss_itero])
# Lxo, Aito, EPrateio, Sio = sub_calcso[:,0,:], sub_calcso[:,1,:], sub_calcso[:,2,:], sub_calcso[:,3,:]
# Aio = np.array([time_avg(times_LRM, aito) for aito in Aito])
# EPio = np.array([integrate(times_LRM, epratei) for epratei in EPrateio])

# for a, b in zip (EPi[0], EPio[2]):
# 	print(np.round(a,5), np.round(b,5))

# ************************** Thermodynamic Quantities for Units ************************************

unit_iter = range(len(units_LRM))
unit_calcs = np.array([get_L_A_EPrate_S_ts(all_upds_LRM[u], Ku_LRM[u]) for u in unit_iter])
Lu, Aut, EPrateu, Su = unit_calcs[:,0,:], unit_calcs[:,1,:], unit_calcs[:,2,:], unit_calcs[:,3,:]
Au = np.array([time_avg(times_LRM, Aut[u]) for u in unit_iter])
EPu = np.array([integrate(times_LRM, eprateu) for eprateu in EPrateu])

def get_subsystem_contribs_to_unit(all_subsystem_contribs, unit):
	return np.array([all_subsystem_contribs[suborder_LRM[s]] for s in units_LRM[uind[unit]]])

Ai_units = np.array([get_subsystem_contribs_to_unit(Ai, u) for u in uind.keys()])
EPi_units = np.array([get_subsystem_contribs_to_unit(EPi, u) for u in uind.keys()])

# Ai_omega, Ai_alpha, Ai_psi = [get_subsystem_contribs_to_unit(Ai, u) for u in ["omega", "alpha", "psi"]]
# EPi_omega, EPi_alpha, EPi_psi = [get_subsystem_contribs_to_unit(EPi, u) for u in ["omega", "alpha", "psi"]]


# ******************************** All Speed Limit Bounds ******************************************

global_tau_bounds = np.array([classical_tau_bound(l, a, ep) 
								for l, a, ep in zip(Lsys[1:], Asys[1:], EPsys[1:])])

print('Calculated global bound')

mpp_tau_bounds = np.array([mpp_tau_bound(l, ais, epis) 
								for l, ais, epis in zip(Lsys[1:], Ai.T[1:], EPi.T[1:])])

print('Calculated mpp bound')

mpp_unit_tau_bounds = np.array([np.array([mpp_tau_bound(l, a, ep) 
								for l, a, ep in zip(Lu[u][1:], Ai_units[u].T[1:], EPi_units[u].T[1:])]) 
									for u in unit_iter])

mpp_phi_tau_bounds, mpp_omega_tau_bounds, mpp_alpha_tau_bounds, mpp_beta_tau_bounds, mpp_psi_tau_bounds = mpp_unit_tau_bounds


print('Calculated mpp unit bounds')


subsystem_local_tau_bounds = np.array([np.array([classical_tau_bound(li, ai, epi) 
								for li, ai, epi in zip(Li[s][1:], Ai[s][1:], EPi[s][1:])]) 
									for s in ss_iter])

sub_1_tau_bounds, sub_2_tau_bounds, sub_3_tau_bounds, sub_4_tau_bounds = subsystem_local_tau_bounds


print('Calculated mpp subsystem-local bounds')


# unit_local_tau_bounds = np.array([np.array([classical_tau_bound(l, a, ep) 
# 								for l, a, ep in zip(Lu[u][1:], Au[u][1:], EPu[u][1:])]) 
# 									for u in unit_iter])

# phi_tau_bounds, omega_tau_bounds, alpha_tau_bounds, beta_tau_bounds, psi_tau_bounds = unit_local_tau_bounds


# print('Calculated unit-local bounds')


print('Calculated all bounds')


# # **************** Thermodynamic Quantities for Units, Unit Structure 2 ****************************

# unit_iter2 = range(len(units_LRM2))
# unit_calcs2 = np.array([get_L_A_EPrate_S_ts(all_upds_LRM2[u], Ku_LRM2[u]) for u in unit_iter2])
# Lu2, Aut2, EPrateu2, Su2=unit_calcs2[:,0,:],unit_calcs2[:,1,:],unit_calcs2[:,2,:],unit_calcs2[:,3,:]
# Au2 = np.array([time_avg(times_LRM, Aut2[u]) for u in unit_iter2])
# EPu2 = np.array([integrate(times_LRM, eprateu) for eprateu in EPrateu2])

# # **************** Thermodynamic Quantities for Odd, Even, and InEx Sums ***************************

# odd_even_LRM = get_odd_even_intersections(units_LRM, suborder_LRM)
# InExEP_LRM = InExf_t(EPu, odd_even_LRM)
# InEx_Info_LRM = InExInfo_t(Ssys, Su, odd_even_LRM)
# DeltaInEx_Info_LRM = Deltaf_t(InEx_Info_LRM)

# odd_even_LRM2 = get_odd_even_intersections(units_LRM2, suborder_LRM)
# InExEP_LRM2 = InExf_t(EPu2, odd_even_LRM2)
# InEx_Info_LRM2 = InExInfo_t(Ssys, Su2, odd_even_LRM2)
# DeltaInEx_Info_LRM2 = Deltaf_t(InEx_Info_LRM2)

# plus_minus_LRM = get_plus_minus(odd_even_LRM, units_LRM) # not repeating for LRM2 because same

# # ******************************** Classical Bounds ************************************************

# # entire system
# class_EP_bounds_LRM = np.array([classical_EP_bound(l, a, t) 
# 								for l, a, t in zip(Lsys, Asys, times_LRM) if t > 0])
# class_tau_bounds_LRM = np.array([classical_tau_bound(l, a, ep) 
# 								for l, a, ep in zip(Lsys, Asys, EPsys) if ep > 0])

# # units in unit structure 1
# class_EPu_bounds_LRM = np.array([np.array([classical_EP_bound(l, a, t) 
# 									for l, a, t in zip(Lu[u], Au[u], times_LRM) if t > 0]) 
# 										for u in unit_iter])
# class_tauu_bounds_LRM = np.array([np.array([classical_tau_bound(l, a, ep) 
# 									for l, a, ep in zip(Lu[u], Au[u], EPu[u]) if ep > 0]) 
# 										for u in unit_iter])

# # units in unit structure 2
# class_EPu_bounds_LRM2 = np.array([np.array([classical_EP_bound(l, a, t) 
# 									for l, a, t in zip(Lu2[u], Au2[u], times_LRM) if t > 0]) 
# 										for u in unit_iter2])
# class_tauu_bounds_LRM2 = np.array([np.array([classical_tau_bound(l, a, ep) 
# 									for l, a, ep in zip(Lu2[u], Au2[u], EPu2[u]) if ep > 0]) 
# 										for u in unit_iter2])

# # Note that bounds for unit phi will not be a lower bound unless time granularity is fine enough

# for l, l1, l2, l3, a, a1, a2, a3 in zip(Lsys, Lu[1], Lu[2], Lu[0], Asys, Au[1], Au[2], Au[0]):
# 	print(l**2/a - l1**2/a1, l**2/a - l2**2/a2, l**2/a - l3**2/a3)

# # for sys, u1, u2, u3 in zip(class_tau_bounds_LRM, *class_tauu_bounds_LRM):
# # 	print(sys, u1, u2, u3)

# # ******************************** In-Ex Sum Bounds ************************************************

# MPP_InExSum_EP_bounds_LRM = MPP_InExSum_bound(Lu, Au, times_LRM, odd_even_LRM)
# MPP_InExSum_tau_bounds_LRM = MPP_InExSum_bound(Lu, Au, EPsys, odd_even_LRM)

# MPP_InExSum_EP_bounds_LRM2 = MPP_InExSum_bound(Lu2, Au2, times_LRM, odd_even_LRM2)
# MPP_InExSum_tau_bounds_LRM2 = MPP_InExSum_bound(Lu2, Au2, EPsys, odd_even_LRM2)

# # ************************** Rearranged In-Ex Sum Bounds *******************************************

# MPP_ReArr_InEx_EP_bounds_LRM = MPP_ReArr_InEx_EP_bound(Lu, Au, times_LRM, EPu, 
# 																DeltaInEx_Info_LRM, odd_even_LRM)
# MPP_ReArr_InEx_tau_bounds_LRM = MPP_ReArr_InEx_tau_bound(Lu, Au, times_LRM, EPu, EPsys, 
# 																DeltaInEx_Info_LRM, odd_even_LRM)

# MPP_ReArr_InEx_EP_bounds_LRM2 = MPP_ReArr_InEx_EP_bound(Lu2, Au2, times_LRM, EPu2, 
# 																DeltaInEx_Info_LRM2, odd_even_LRM2)
# MPP_ReArr_InEx_tau_bounds_LRM2 = MPP_ReArr_InEx_tau_bound(Lu2, Au2, times_LRM, EPu2, EPsys, 
# 																DeltaInEx_Info_LRM2, odd_even_LRM2)

# MPP_PlusMinus_InEx_EP_bounds_LRM = MPP_ReArr_InEx_EP_bound(Lu, Au, times_LRM, EPu, 
# 																DeltaInEx_Info_LRM, plus_minus_LRM)
# MPP_PlusMinus_InEx_tau_bounds_LRM = MPP_ReArr_InEx_tau_bound(Lu, Au, times_LRM, EPu, EPsys, 
# 																DeltaInEx_Info_LRM, plus_minus_LRM)

# # GET CONTRIBUTIONS FROM UNITS FOR THIS

# # ****************************** VDFT Regular Sum Bounds *******************************************

# LtermsB1 = [[l1, l2] for l1, l2 in zip(Lu[uind['omega']], Lu[uind['alpha']])]
# AtermsB1 = [[a1, a2] for a1, a2 in zip(Au[uind['omega']], Au[uind['alpha']])]
# EPtermsB1 = [ep1 + ep2 for ep1, ep2 in zip(EPu[uind['omega']], EPu[uind['alpha']])]

# MPP_B1_EP_bounds_LRM = np.array([MPP_RegSum_bound(ls, acts, t) 
# 								for ls, acts, t in zip(LtermsB1, AtermsB1, times_LRM) if t > 0])
# MPP_B1_tau_bounds_LRM = np.array([MPP_RegSum_bound(ls, acts, ep) 
# 								for ls, acts, ep in zip(LtermsB1, AtermsB1, EPtermsB1) if ep > 0])

# MPP_B1_EP_contribs = np.array([MPP_RegSum_contribs(ls, acts, t) 
# 								for ls, acts, t in zip(LtermsB1, AtermsB1, times_LRM) if t > 0]).T
# MPP_B1_tau_contribs = np.array([MPP_RegSum_contribs(ls, acts, ep) 
# 								for ls, acts, ep in zip(LtermsB1, AtermsB1, EPtermsB1) if ep > 0]).T

# LtermsB2 = [[l1, l2] for l1, l2 in zip(Lu2[uind2['beta']], Lu2[uind2['alpha']])]
# AtermsB2 = [[a1, a2] for a1, a2 in zip(Au2[uind2['beta']], Au2[uind2['alpha']])]
# EPtermsB2 = [ep1 + ep2 for ep1, ep2 in zip(EPu2[uind2['beta']], EPu2[uind2['alpha']])]

# MPP_B2_EP_bounds_LRM = np.array([MPP_RegSum_bound(ls, acts, t) 
# 								for ls, acts, t in zip(LtermsB2, AtermsB2, times_LRM) if t > 0])
# MPP_B2_tau_bounds_LRM = np.array([MPP_RegSum_bound(ls, acts, ep) 
# 								for ls, acts, ep in zip(LtermsB2, AtermsB2, EPtermsB2) if ep > 0])

# MPP_B2_EP_contribs = np.array([MPP_RegSum_contribs(ls, acts, t) 
# 								for ls, acts, t in zip(LtermsB2, AtermsB2, times_LRM) if t > 0]).T
# MPP_B2_tau_contribs = np.array([MPP_RegSum_contribs(ls, acts, ep) 
# 								for ls, acts, ep in zip(LtermsB2, AtermsB2, EPtermsB2) if ep > 0]).T

# # ************************************* CFT Bounds *************************************************
# MPP_CFT_EP_alphaphi = class_EPu_bounds_LRM[uind['phi']] - DeltaInEx_Info_LRM[1:]

# LtermsCFT, AtermsCFT = Lu[uind['phi']], Au[uind['phi']]
# EPpItermsCFT = EPu[uind['alpha']] + DeltaInEx_Info_LRM
# MPP_CFT_tau_alphaphi = np.array([classical_bound(l, a, x) 
# 				for l, a, x, t in zip(LtermsCFT, AtermsCFT, EPpItermsCFT, times_LRM) if t > 0])

# # *********************************** MPP EPR Bounds ***********************************************
# deltat = timestep/10

# Ldot_phi = Ldot_t(all_upds_LRM[uind['phi']], Ku_LRM[uind['phi']], deltat)
# Ldot_omega_phi = Ldot_t(all_upds_LRM[uind['omega']], K_omega_phi, deltat)
# Ldot_N_omega = Ldot_t(jpds_LRM, K_N_omega, deltat)

# At_omega_phi = A_t(all_upds_LRM[uind['omega']], K_omega_phi)
# At_N_omega = A_t(jpds_LRM, K_N_omega)

# Ldot_ID1 = [[l1, l2, l3] for l1, l2, l3 in zip(Ldot_phi, Ldot_omega_phi, Ldot_N_omega)]
# At_ID1 = [[a1, a2, a3] for a1, a2, a3 in zip(Aut[uind['phi']], At_omega_phi, At_N_omega)]

# MPP_EPR_ID1 = np.array([MPP_RegSum_bound(ls, acts, 1) 
# 								for ls, acts, t in zip(Ldot_ID1, At_ID1, times_LRM) if t > 0]).T


# Ldot_beta_phi = Ldot_t(all_upds_LRM2[uind2['beta']], K_beta_phi, deltat)
# Ldot_omega_beta = Ldot_t(all_upds_LRM2[uind2['omega']], K_omega_beta, deltat)

# At_beta_phi = A_t(all_upds_LRM2[uind2['beta']], K_beta_phi)
# At_omega_beta = A_t(all_upds_LRM2[uind2['omega']], K_omega_beta)

# Ldot_ID2 = [[l1, l2, l3, l4] for l1, l2, l3, l4 in zip(Ldot_phi, Ldot_beta_phi, Ldot_omega_beta, Ldot_N_omega)]
# At_ID2 = [[a1, a2, a3, a4] for a1, a2, a3, a4 in zip(Aut[uind['phi']], At_beta_phi, At_omega_beta, At_N_omega)]

# MPP_EPR_ID2 = np.array([MPP_RegSum_bound(ls, acts, 1) 
# 								for ls, acts, t in zip(Ldot_ID2, At_ID2, times_LRM) if t > 0]).T

















# *********************************** Sanity Checks ************************************************

# odd_EP, even_EP = group_vals_t(odd_even_LRM[0], EPu), group_vals_t(odd_even_LRM[1], EPu)

# to check if iterative decomposition holds, these should be equal within some desired precision
# for ep, epno, epoopr, epopr in zip(EPsys, EP_N_omega, EP_omega_phi, EPu[0]):
# 	print(ep, epno+epoopr+epopr)

# Lsys = L_t(jpds_LRM)
# Asyst = A_t(jpds_LRM, K_LRM)
# EPratesys = EPrate_t(jpds_LRM, K_LRM)
# Ssys = S_t(jpds_LRM)

# Lu = np.array([L_t(all_upds_LRM[u]) for u in unit_iter])
# Aut = np.array([A_t(all_upds_LRM[u], Ku_LRM[u]) for u in unit_iter])
# EPrateu = np.array([EPrate_t(all_upds_LRM[u], Ku_LRM[u]) for u in unit_iter])
# Su = np.array([S_t(all_upds_LRM[u]) for u in unit_iter])

# L_N_omega = L_t(jpds_N_omega)
# A_N_omegat = A_t(jpds_N_omega, K_N_omega)
# EPrate_N_omega = EPrate_t(jpds_N_omega, K_N_omega)

# L_omega_phi = L_t(omegapds_omega_phi)
# A_omega_phit = A_t(omegapds_omega_phi, K_omega_phi)
# EPrate_omega_phi = EPrate_t(omegapds_omega_phi, K_omega_phi)

# to check if in-ex sums hold these should be equal within some desired precision
# for ep, inexep, alphainexinfo in zip(EPsys, InEx_EP, DeltaInEx_Info):
# 	print(ep, inexep - alphainexinfo)

# for lu, au in zip(Lu, Au):
# 	print(np.array([l**2 for l in lu]))
# 	print(au)
# 	print("")

# for a, epu in zip(class_EPu_bounds_LRM, EPu):
# 	for bd, ep in zip(a, epu[1:]):
# 		print(bd, ep)
# 	print("")

# for a in class_tauu_bounds_LRM:
# 	for bd, t in zip(a, times_LRM[1:]):
# 		print(bd, t)
# 	print("")

# for bd, ep in zip(class_EPu_bounds_LRM[0], EPu[0][1:]):
# 	print(bd, ep)
# print("")


# for bd, t in zip(class_tauu_bounds_LRM[0], times_LRM[1:]):
# 	print(bd, t)
# print("")

# print('L_omega_phi')
# for l in L_omega_phi:
# 	print(l)
# print('')

# print('L_omega')
# for l in Lu[uind['omega']]:
# 	print(l)
# print('')

# print('Lsys')
# for l in Lsys:
# 	print(l)
# print('')

# print('A_omega_phi')
# for l in A_omega_phi:
# 	print(l)
# print('')

# print('A_omega')
# for l in Au[uind['omega']]:
# 	print(l)
# print('')

# print('Asys')
# for l in Asys:
# 	print(l)
# print('')

# print('Asyst')
# for a in Asyst:
# 	print(a)
# print('')

# print('<Asys>')
# for a in Asys:
# 	print(a)
# print('')

# print('EPsys')
# for ep in EPsys:
# 	print(ep)

# print('EPratesys')
# for epr, ep in zip(EPratesys, EPsys):
# 	print(epr, ep)

# print('EP integ vs. EP direct')
# for ep, epr in zip(EPsys, EPratesys):
# 	print(ep, epr)

# print('******************')

# print('L_N_omega')
# for l in L_N_omega:
# 	print(l)
# print('')

# print('A_N_omegat')
# for a in A_N_omegat:
# 	print(a)
# print('')

# print('<A_N_omega>')
# for a in A_N_omega:
# 	print(a)
# print('')

# print('EP_N_omega')
# for ep in EP_N_omega:
# 	print(ep)

# for u in range(numunits):
# 	print('******************')

# 	print('L', units_LRM[u])
# 	for l in Lu[u]:
# 		print(l)
# 	print('')

# 	print('A', units_LRM[u])
# 	for a in Au[u]:
# 		print(a)
# 	print('')

# 	print('EP', units_LRM[u])
# 	for ep in EPu[u]:
# 		print(ep)

	# print('EPrate', units_LRM[u])
	# for epr in EPrateu[u]:
	# 	print(epr)

# print('EP bounds')
# for b, mppb, ep in zip(class_EP_bounds_LRM, MPP_EP_bounds_LRM, EPsys[1:]):
# 	print(np.round(b, 8), np.round(mppb, 8), np.round(ep, 8))
# print('')

# print('tau bounds')
# for b, mppb, t in zip(class_tau_bounds_LRM, MPP_tau_bounds_LRM, times_LRM[1:]):
# 	print(np.round(b, 8), np.round(mppb, 8), np.round(t, 8))
# print('')

# for classbd, epbound, ep in zip(class_EP_bounds_LRM, MPP_InExSum_EP_bounds_LRM, EPsys[1:]):
# 	print(classbd, epbound, ep)

# print("")

# for classbd, mppbound, t in zip(class_tau_bounds_LRM, MPP_InExSum_tau_bounds_LRM, times_LRM[1:]):
# 	print(classbd, mppbound, t)

# # *********************** Thermodynamic Quantities Under K(N\w) ************************************ 

# L_N_omega, A_N_omegat, EPrate_N_omega, S_N_omega = get_L_A_EPrate_S_ts(jpds_N_omega, K_N_omega)
# A_N_omega = time_avg(times_LRM, A_N_omegat)
# EP_N_omega = integrate(times_LRM, EPrate_N_omega)

# # *********************** Thermodynamic Quantities Under K(w\w') ***********************************

# L_omega_phi, A_omega_phit, EPrate_omega_phi, S_omega_phi = get_L_A_EPrate_S_ts(
# 															omegapds_omega_phi, K_omega_phi)
# A_omega_phi = time_avg(times_LRM, A_omega_phit)
# EP_omega_phi = integrate(times_LRM, EPrate_omega_phi)

# # *********************** Thermodynamic Quantities Under K(N\a) ************************************ 

# L_N_alpha, A_N_alphat, EPrate_N_alpha, S_N_alpha = get_L_A_EPrate_S_ts(jpds_N_alpha, K_N_alpha)
# A_N_alpha = time_avg(times_LRM, A_N_alphat)
# EP_N_alpha = integrate(times_LRM, EPrate_N_alpha)

# # *********************** Thermodynamic Quantities Under K(a\w') ***********************************

# L_alpha_phi, A_alpha_phit, EPrate_alpha_phi, S_alpha_phi = get_L_A_EPrate_S_ts(
# 															alphapds_alpha_phi, K_alpha_phi)
# A_alpha_phi = time_avg(times_LRM, A_alpha_phit)
# EP_alpha_phi = integrate(times_LRM, EPrate_alpha_phi)

# # *********************** Thermodynamic Quantities Under K(N\g) ************************************ 

# # this is costly so I'm not going to repeat these calculations
# L_N_omega, A_N_omegat, EPrate_N_omega, S_N_omega = L_N_omega, A_N_omegat, EPrate_N_omega, S_N_omega
# A_N_omega, EP_N_omega = A_N_omega, EP_N_omega

# # *********************** Thermodynamic Quantities Under K(N\e) ************************************ 

# L_N_beta, A_N_betat, EPrate_N_beta, S_N_beta = get_L_A_EPrate_S_ts(jpds_N_beta, 
# 																					K_N_beta)
# A_N_beta = time_avg(times_LRM, A_N_betat)
# EP_N_beta= integrate(times_LRM, EPrate_N_beta)

# # *********************** Thermodynamic Quantities Under K(g\e) ***********************************

# L_omega_beta, A_omega_betat, EPrate_omega_beta, S_omega_beta = get_L_A_EPrate_S_ts(
# 															omegapds_omega_beta, K_omega_beta)
# A_omega_beta = time_avg(times_LRM, A_omega_betat)
# EP_omega_beta = integrate(times_LRM, EPrate_omega_beta)

# # *********************** Thermodynamic Quantities Under K(e\b) ***********************************

# L_beta_phi, A_beta_phit, EPrate_beta_phi, S_beta_phi = get_L_A_EPrate_S_ts(
# 															betapds_beta_phi, K_beta_phi)
# A_beta_phi = time_avg(times_LRM, A_beta_phit)
# EP_beta_phi = integrate(times_LRM, EPrate_beta_phi)

# # ****************** Iterative Decomposition Bounds, Unit Structure 1 ******************************

# Lterms1a = [[l1, l2, l3] for l1, l2, l3 in zip(Lu[uind['phi']], L_omega_phi, L_N_omega)]
# Aterms1a = [[a1, a2, a3] for a1, a2, a3 in zip(Au[uind['phi']], A_omega_phi, A_N_omega)]

# MPP_ID1a_EP_bounds_LRM = np.array([MPP_IterDecomp_EP_bound(ls, acts, t) 
# 								for ls, acts, t in zip(Lterms1a, Aterms1a, times_LRM) if t > 0])
# MPP_ID1a_tau_bounds_LRM = np.array([MPP_IterDecomp_tau_bound(ls, acts, ep) 
# 								for ls, acts, ep in zip(Lterms1a, Aterms1a, EPsys) if ep > 0])

# MPP_ID1a_EP_contribs = np.array([MPP_IterDecomp_EP_contribs(ls, acts, t) 
# 								for ls, acts, t in zip(Lterms1a, Aterms1a, times_LRM) if t > 0]).T
# MPP_ID1a_tau_contribs = np.array([MPP_IterDecomp_EP_contribs(ls, acts, ep) 
# 								for ls, acts, ep in zip(Lterms1a, Aterms1a, EPsys) if ep > 0]).T

# Lterms1b = [[l1, l2, l3] for l1, l2, l3 in zip(Lu[uind['phi']], L_alpha_phi, L_N_alpha)]
# Aterms1b = [[a1, a2, a3] for a1, a2, a3 in zip(Au[uind['phi']], A_alpha_phi, A_N_alpha)]

# MPP_ID1b_EP_bounds_LRM = np.array([MPP_IterDecomp_EP_bound(ls, acts, t) 
# 								for ls, acts, t in zip(Lterms1b, Aterms1b, times_LRM) if t > 0])
# MPP_ID1b_tau_bounds_LRM = np.array([MPP_IterDecomp_tau_bound(ls, acts, ep) 
# 								for ls, acts, ep in zip(Lterms1b, Aterms1b, EPsys) if ep > 0])

# MPP_ID1b_EP_contribs = np.array([MPP_IterDecomp_EP_contribs(ls, acts, t) 
# 								for ls, acts, t in zip(Lterms1b, Aterms1b, times_LRM) if t > 0]).T
# MPP_ID1b_tau_contribs = np.array([MPP_IterDecomp_EP_contribs(ls, acts, ep) 
# 								for ls, acts, ep in zip(Lterms1b, Aterms1b, EPsys) if ep > 0]).T

# # ****************** Iterative Decomposition Bounds, Unit Structure 2 ******************************

# Lterms2a = [[l1, l2, l3] for l1, l2, l3 in zip(Lu2[uind2['beta']], L_omega_beta, L_N_omega)]
# Aterms2a = [[a1, a2, a3] for a1, a2, a3 in zip(Au2[uind2['beta']], A_omega_beta, A_N_omega)]

# MPP_ID2a_EP_bounds_LRM = np.array([MPP_IterDecomp_EP_bound(ls, acts, t) 
# 								for ls, acts, t in zip(Lterms2a, Aterms2a, times_LRM) if t > 0])
# MPP_ID2a_tau_bounds_LRM = np.array([MPP_IterDecomp_tau_bound(ls, acts, ep) 
# 								for ls, acts, ep in zip(Lterms2a, Aterms2a, EPsys) if ep > 0])

# MPP_ID2a_EP_contribs = np.array([MPP_IterDecomp_EP_contribs(ls, acts, t) 
# 								for ls, acts, t in zip(Lterms2a, Aterms2a, times_LRM) if t > 0]).T
# MPP_ID2a_tau_contribs = np.array([MPP_IterDecomp_EP_contribs(ls, acts, ep) 
# 								for ls, acts, ep in zip(Lterms2a, Aterms2a, EPsys) if ep > 0]).T

# Lterms2b = [[l1, l2, l3, l4] for l1, l2, l3, l4 in zip(Lu2[uind2['phi']], L_beta_phi, 
# 																		L_omega_beta, L_N_omega)]
# Aterms2b = [[a1, a2, a3, a4] for a1, a2, a3, a4 in zip(Au2[uind2['phi']], A_beta_phi, 
# 																		A_omega_beta, A_N_omega)]

# MPP_ID2b_EP_bounds_LRM = np.array([MPP_IterDecomp_EP_bound(ls, acts, t) 
# 								for ls, acts, t in zip(Lterms2b, Aterms2b, times_LRM) if t > 0])
# MPP_ID2b_tau_bounds_LRM = np.array([MPP_IterDecomp_tau_bound(ls, acts, ep) 
# 								for ls, acts, ep in zip(Lterms2b, Aterms2b, EPsys) if ep > 0])

# MPP_ID2b_EP_contribs = np.array([MPP_IterDecomp_EP_contribs(ls, acts, t) 
# 								for ls, acts, t in zip(Lterms2b, Aterms2b, times_LRM) if t > 0]).T
# MPP_ID2b_tau_contribs = np.array([MPP_IterDecomp_EP_contribs(ls, acts, ep) 
# 								for ls, acts, ep in zip(Lterms2b, Aterms2b, EPsys) if ep > 0]).T

# Lterms2c = [[l1, l2, l3] for l1, l2, l3 in zip(Lu2[uind2['phi']], L_beta_phi, L_N_beta)]
# Aterms2c = [[a1, a2, a3] for a1, a2, a3 in zip(Au2[uind2['phi']], A_beta_phi, A_N_beta)]

# MPP_ID2c_EP_bounds_LRM = np.array([MPP_IterDecomp_EP_bound(ls, acts, t) 
# 								for ls, acts, t in zip(Lterms2c, Aterms2c, times_LRM) if t > 0])
# MPP_ID2c_tau_bounds_LRM = np.array([MPP_IterDecomp_tau_bound(ls, acts, ep) 
# 								for ls, acts, ep in zip(Lterms2c, Aterms2c, EPsys) if ep > 0])

# MPP_ID2c_EP_contribs = np.array([MPP_IterDecomp_EP_contribs(ls, acts, t) 
# 								for ls, acts, t in zip(Lterms2c, Aterms2c, times_LRM) if t > 0]).T
# MPP_ID2c_tau_contribs = np.array([MPP_IterDecomp_EP_contribs(ls, acts, ep) 
# 								for ls, acts, ep in zip(Lterms2c, Aterms2c, EPsys) if ep > 0]).T


# for bound, ep, t in zip(class_EP_bounds_LRM, EPsys[1:], times_LRM):
# 	print (t, bound*t, ep, t*bound/ep)

# print("")

# numunits = len(EPu2)

# for u in range(numunits):
# 	print(u)
# 	for ubound, uep, t in zip(class_EPu_bounds_LRM2[u], EPu2[u][1:], times_LRM):
# 		print(t, ubound*t, uep, t*ubound/uep)
# 	print("")

# for bound, t in zip(class_tau_bounds_LRM, times_LRM[1:]):
# 	print (t, bound)

# print("")

# numunits = len(EPu2)

# for u in range(numunits):
# 	print(u)
# 	for ubound, t in zip(class_tauu_bounds_LRM2[u], times_LRM[1:]):
# 		print(t, ubound)
# 	print("")

# # should not be equal
# print("L")
# for l, lphi, lomega, lalpha in zip (Lsys, Lu[0], Lu[1], Lu[2]):
# 	print(l, lomega + lalpha - lphi)

# # should be equal
# print("A")
# for a, aphi, aomega, aalpha in zip (Asys, Au[0], Au[1], Au[2]):
# 	print(a, aomega + aalpha - aphi)