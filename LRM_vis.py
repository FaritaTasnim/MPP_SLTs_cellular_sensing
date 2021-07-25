from tools_vis import *
from LRM_calc import *
from LRM_pdist import *

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Source Sans Pro']})

# actuals = [EPsys[1:], times[1:]*1000]
# actual_names = [r'$\langle \sigma \rangle$', r'$\tau$ ' + '(ms)']

actual = times_LRM[1:]*1000
actual_name = r'$\tau$ ' + '(ms)'

# colors = ['c', 'y', 'm']
palettes = ['rocket', 'mako', 'flare', 'crest', 'bright', 'hls']
names = [r'$x$', r'$n_{b,1}$', r'$n_y$', r'$n_{b,2}$']

all_bounds = {	'actual': times_LRM[1:]*1000,
			'global': global_tau_bounds*1000, 
			'mpp': mpp_tau_bounds*1000,
			'sub_1': sub_1_tau_bounds*1000,
			'sub_2': sub_2_tau_bounds*1000,
			'sub_3': sub_3_tau_bounds*1000,
			'sub_4': sub_4_tau_bounds*1000,
			# 'local_phi': phi_tau_bounds*1000,
			# 'local_omega': omega_tau_bounds*1000,
			# 'local_alpha': alpha_tau_bounds*1000,
			# 'local_beta': beta_tau_bounds*1000,
			# 'local_psi': psi_tau_bounds*1000,
			'mpp_phi': mpp_phi_tau_bounds*1000,
			'mpp_omega': mpp_omega_tau_bounds*1000,
			'mpp_alpha': mpp_alpha_tau_bounds*1000,
			'mpp_beta': mpp_beta_tau_bounds*1000,
			'mpp_psi': mpp_psi_tau_bounds*1000}

all_names = {	'actual': r'$\tau$ ' + '(ms)',
				'global': 'global CSL',
				'mpp': 'Eq. (22)',
				'sub_1': 'Eq. (24), i=1',
				'sub_2': 'Eq. (24), i=2',
				'sub_3': 'Eq. (24), i=3',
				'sub_4': 'Eq. (24), i=4',
				# 'local_phi': r'$\phi$' + '-local CSL',
				# 'local_omega': r'$\omega$' + '-local CSL',
				# 'local_alpha': r'$\alpha$' + '-local CSL',
				# 'local_beta': r'$\beta$' + '-local CSL',
				# 'local_psi': r'$\psi$' + '-local CSL',
				'mpp_phi': 'Eq. (23), ' + r'$\phi$',
				'mpp_omega': 'Eq. (23), ' + r'$\omega$', 
				'mpp_alpha': 'Eq. (23), ' + r'$\alpha$', 
				'mpp_beta': 'Eq. (23), ' + r'$\beta$', 
				'mpp_psi': 'Eq. (23), ' + r'$\psi$'}

colors = {	'actual': '#2f4f4f',
			'global': '#228b22', 
			'mpp': '#ff1493', 
			'sub_1': '#ffd700',
			'sub_2': '#00bfff',
			'sub_3': '#ff8c00',
			'sub_4': '#00ff00',
			# 'local_phi': '#42d4f4',
			# 'local_omega': '#ff00ff',
			# 'local_alpha': '#4169e1',
			# 'local_beta': '#911eb4', 
			# 'local_psi': '#adff2f',
			'mpp_phi': '#911eb4',
			'mpp_omega': '#dda0dd',
			'mpp_alpha': '#7fffd4',
			'mpp_beta': '#4169e1',
			'mpp_psi': '#808000'}

subplot_1_keys = ['global', 'mpp', 'sub_1', 'sub_2', 'sub_3', 'sub_4']
subplot_2_keys = ['global', 'mpp_omega', 'sub_1', 'sub_2', 'sub_3']
subplot_3_keys = ['global', 'mpp_alpha', 'sub_3', 'sub_4']
subplot_4_keys = ['global', 'mpp_phi', 'sub_3']
subplot_5_keys = ['global', 'mpp_beta', 'sub_2', 'sub_3']
subplot_6_keys = ['global', 'mpp_psi', 'sub_2', 'sub_3', 'sub_4']
subplot_7_keys = ['global', 'mpp', 'mpp_phi', 'mpp_omega', 'mpp_alpha', 'mpp_beta', 'mpp_psi']

first_plot_keys = [subplot_1_keys, subplot_7_keys]
first_plot_tags = ['subs', 'mpps']

first_bounds_to_plot = []
first_bound_names_to_plot = []
first_bound_colors_to_plot = []

for plot_keys, plot_tag in zip(first_plot_keys, first_plot_tags):
	first_bounds_to_plot.append([all_bounds[x] for x in plot_keys])
	first_bound_names_to_plot.append([all_names[x] for x in plot_keys])
	first_bound_colors_to_plot.append([colors[x] for x in ['actual'] + plot_keys])
	# bounds = [[all_bounds[x] for x in plot_keys]]
	# names = [[all_names[x] for x in plot_keys]]
	# plot_colors = [colors[x] for x in ['actual'] + plot_keys]
	# plot_bounds(times_LRM*1000, actuals, actual_names, bounds, names, plot_colors, plot_tag)

plot_bounds(times_LRM*1000, actual, actual_name, first_bounds_to_plot, first_bound_names_to_plot, first_bound_colors_to_plot, 'main', (10,4))



second_plot_keys = [subplot_2_keys, subplot_3_keys, subplot_4_keys, subplot_5_keys, subplot_6_keys]
second_plot_tags = ['omega', 'alpha', 'phi', 'beta', 'psi']

second_bounds_to_plot = []
second_bound_names_to_plot = []
second_bound_colors_to_plot = []

for plot_keys, plot_tag in zip(second_plot_keys, second_plot_tags):
	second_bounds_to_plot.append([all_bounds[x] for x in plot_keys])
	second_bound_names_to_plot.append([all_names[x] for x in plot_keys])
	second_bound_colors_to_plot.append([colors[x] for x in ['actual'] + plot_keys])
	# bounds = [[all_bounds[x] for x in plot_keys]]
	# names = [[all_names[x] for x in plot_keys]]
	# plot_colors = [colors[x] for x in ['actual'] + plot_keys]
	# plot_bounds(times_LRM*1000, actuals, actual_names, bounds, names, plot_colors, plot_tag)

plot_bounds(times_LRM*1000, actual, actual_name, second_bounds_to_plot, second_bound_names_to_plot, second_bound_colors_to_plot, 'units', (25,4))

# ************************************************************
# plot_1_keys = ['global', 'mpp', 'sub_1', 'sub_2', 'sub_3', 'sub_4']
# plot_2_keys = ['global', 'local_phi', 'local_omega', 'local_alpha', 'local_beta', 'local_psi']
# plot_3_keys = ['global', 'local_phi', 'mpp_phi', 'sub_3']
# plot_4_keys = ['global', 'local_omega', 'mpp_omega', 'sub_1', 'sub_2', 'sub_3']
# plot_5_keys = ['global', 'local_alpha', 'mpp_alpha', 'sub_3', 'sub_4']
# plot_6_keys = ['global', 'local_beta', 'mpp_beta', 'sub_2', 'sub_3']
# plot_7_keys = ['global', 'local_psi', 'mpp_psi', 'sub_2', 'sub_3', 'sub_4']
# plot_8_keys = ['global', 'mpp', 'mpp_phi', 'mpp_omega', 'mpp_alpha', 'mpp_beta', 'mpp_psi']
# all_plot_keys = [plot_1_keys, plot_2_keys, plot_3_keys, plot_4_keys, 
# 					plot_5_keys, plot_6_keys, plot_7_keys, plot_8_keys]
# all_plot_tags = ['N', 'globalvslocal', 'phi', 'omega', 'alpha', 'beta', 'psi', 'mpps']



# for plot_keys, plot_tag in zip(all_plot_keys, all_plot_tags):
# 	bounds = [[all_bounds[x] for x in plot_keys]]
# 	names = [[all_names[x] for x in plot_keys]]
# 	plot_colors = [colors[x] for x in ['actual'] + plot_keys]
# 	plot_bounds(times_LRM*1000, actuals, actual_names, bounds, names, plot_colors, plot_tag)



#***************************************************************
# N_bds = 	[[global_tau_bounds*1000,
# 						mpp_tau_bounds*1000,
# 						sub_1_tau_bounds*1000,
# 						sub_2_tau_bounds*1000,
# 						sub_3_tau_bounds*1000,
# 						sub_4_tau_bounds*1000]]

# N_bd_names = [['global CSL',
# 							'Eq. (40)',
# 							'Eq. (49), i=1',
# 							'Eq. (49), i=2',
# 							'Eq. (49), i=3',
# 							'Eq. (49), i=4']]

# plot_bounds(times_LRM*1000, actuals, actual_names, N_bds, N_bd_names, palettes[4], 'N')



# globalvslocal_bds = 	[[global_tau_bounds*1000,
# 						omega_tau_bounds*1000,
# 						alpha_tau_bounds*1000,
# 						beta_tau_bounds*1000,
# 						psi_tau_bounds*1000]]

# globalvslocal_bd_names = [['global CSL',
# 							r'$\omega$' + '-local CSL',
# 							r'$\alpha$' + '-local CSL',
# 							r'$\beta$' + '-local CSL',
# 							r'$\psi$' + '-local CSL']]

# plot_bounds(times_LRM*1000, actuals, actual_names, globalvslocal_bds, globalvslocal_bd_names, palettes[4], 'globalvslocal')



# omega_bds = 	[[omega_tau_bounds*1000,
# 						mpp_omega_tau_bounds*1000,
# 						sub_1_tau_bounds*1000,
# 						sub_2_tau_bounds*1000,
# 						sub_3_tau_bounds*1000]]

# omega_bd_names = [[r'$\omega$' + '-local CSL',
# 							'Eq. (41), ' + r'$\omega$',
# 							'Eq. (49), i=1',
# 							'Eq. (49), i=2',
# 							'Eq. (49), i=3']]

# plot_bounds(times_LRM*1000, actuals, actual_names, omega_bds, omega_bd_names, palettes[4], 'omega')



# alpha_bds = 	[[alpha_tau_bounds*1000,
# 						mpp_alpha_tau_bounds*1000,
# 						sub_3_tau_bounds*1000,
# 						sub_4_tau_bounds*1000]]

# alpha_bd_names = [[r'$\alpha$' + '-local CSL',
# 						'Eq. (41), ' + r'$\alpha$',
# 						'Eq. (49), i=3',
# 						'Eq. (49), i=4']]
# #i values match that in figure in paper


# plot_bounds(times_LRM*1000, actuals, actual_names, alpha_bds, alpha_bd_names, palettes[4], 'alpha')



# beta_bds = 	[[beta_tau_bounds*1000,
# 						mpp_beta_tau_bounds*1000,
# 						sub_2_tau_bounds*1000,
# 						sub_3_tau_bounds*1000]]

# beta_bd_names = [[r'$\beta$' + '-local CSL',
# 						'Eq. (41), ' + r'$\beta$',
# 						'Eq. (49), i=2',
# 						'Eq. (49), i=3']]
# #i values match that in figure in paper


# plot_bounds(times_LRM*1000, actuals, actual_names, beta_bds, beta_bd_names, palettes[4], 'beta')




# psi_bds = 	[[psi_tau_bounds*1000,
# 						mpp_psi_tau_bounds*1000,
# 						sub_2_tau_bounds*1000,
# 						sub_3_tau_bounds*1000,
# 						sub_4_tau_bounds*1000]]

# psi_bd_names = [[r'$\psi$' + '-local CSL',
# 						'Eq. (41), ' + r'$\psi$',
# 						'Eq. (49), i=2',
# 						'Eq. (49), i=3',
# 						'Eq. (49), i=4']]
# #i values match that in figure in paper


# plot_bounds(times_LRM*1000, actuals, actual_names, psi_bds, psi_bd_names, palettes[4], 'psi')









# *******************************************************
# plot_bounds(times_LRM*1000, actuals, actual_names, 
# 	[N_bds, globalvslocal_bds, omega_bds, alpha_bds, beta_bds, psi_bds], 
# 	[N_bd_names, globalvslocal_bd_names, omega_bd_names, alpha_bd_names, beta_bd_names, psi_bd_names], 
# 	palettes[4], 'N')

# VDFT Bounds

# first checking if \sigma \geq \sum_{b_i} \sigma^{b_i}
# a = [EPsys[1:], -1*DeltaInEx_Info[1:]]
# a_names = [r'$\langle \sigma \rangle$', '-I']
# EP_B = [[EPu[uind['omega']][1:] + EPu[uind['alpha']][1:]],
# 		[EPu[uind['phi']][1:]]]
# EP_B_names = [['B1'],[r'$\langle \sigma^{\phi} \rangle$']]
# plot_bounds(times*1000, a, a_names, EP_B, EP_B_names, palettes[4])

# EP_B = [[EPu2[uind2['beta']][1:] + EPu2[uind2['alpha']][1:]],
# 		[EPu[uind['phi']][1:]]]
# EP_B_names = [['B2'],[r'$\langle \sigma^{\phi} \rangle$']]
# plot_bounds(times*1000, a, a_names, EP_B, EP_B_names, palettes[4])

# FOUND OUT ONLY Neitherrr is valid

# B1_labels = [r"$\omega$", r"$\alpha$"]
# B2_labels = [r"$\beta$", r"$\alpha$"]

# bounds_B = 	[[class_EP_bounds, 
# 					# MPP_B1_EP_bounds, 
# 					MPP_B2_EP_bounds],
# 				[class_tau_bounds*1000, 
# 					# MPP_B1_tau_bounds*1000, 
# 					MPP_B2_tau_bounds*1000]]

# bound_names_B = 	[[r'$\langle \sigma \rangle$' + ' bound', 
# 						# r'$\langle \sigma \rangle$' +' MPP B1 bound',
# 						r'$\langle \sigma \rangle$' +' MPP B2 bound'],
# 					[r'$\tau$' + ' bound', 
# 						# r'$\tau$' + ' MPP B1 bound',
# 						r'$\tau$' + ' MPP B2 bound']]

# bounds
# plot_contribs_to_bound(times*1000, MPP_B1_EP_contribs, B1_labels, 'B1 EP', palettes[2])
# plot_contribs_to_bound(times*1000, MPP_B1_tau_contribs*1000, B1_labels, 'B1 tau', palettes[2])

# plot_contribs_to_bound(times*1000, MPP_B2_EP_contribs, B2_labels, 'B2 EP', palettes[2])
# plot_contribs_to_bound(times*1000, MPP_B2_tau_contribs*1000, B2_labels, 'B2 tau', palettes[2])

# plot_bounds(times*1000, actuals, actual_names, bounds_B, bound_names_B, palettes[2])

#InEx, OddEven ReArr InEx, PlusMinus ReArr InEx Sum bounds

# bounds_InEx = 	[[class_EP_bounds, 
# 					MPP_InExSum_EP_bounds,
# 					MPP_InExSum_EP_bounds2],
# 				[class_tau_bounds*1000,
# 					MPP_InExSum_tau_bounds*1000,
# 					MPP_InExSum_tau_bounds2*1000]]

# bound_names_InEx = 	[[r'$\langle \sigma \rangle$' + ' bound',
# 						r'$\langle \sigma \rangle$' +' MPP InEx bound',
# 						r'$\langle \sigma \rangle$' +' MPP InEx2 bound'],
# 					[r'$\tau$' + ' bound',
# 						r'$\tau$' + ' MPP InEx bound',
# 						r'$\tau$' + ' MPP InEx2 bound']]

# bounds_RearrInEx = 	[[class_EP_bounds, 
# 						MPP_ReArr_InEx_EP_bounds,
# 						MPP_ReArr_InEx_EP_bounds2],
# 					[class_tau_bounds*1000,
# 						MPP_ReArr_InEx_tau_bounds*1000,
# 						MPP_ReArr_InEx_tau_bounds2*1000]]

# bound_names_RearrInEx = [[r'$\langle \sigma \rangle$' + ' bound',
# 							r'$\langle \sigma \rangle$' +' MPP ReArr InEx bound',
# 							r'$\langle \sigma \rangle$' +' MPP ReArr InEx2 bound'],
# 						[r'$\tau$' + ' bound',
# 							r'$\tau$' + ' MPP ReArr InEx bound',
# 							r'$\tau$' + ' MPP ReArr InEx2 bound']]

# bounds_PlusMinusInEx = 	[[class_EP_bounds, 
# 						MPP_PlusMinus_InEx_EP_bounds],
# 					[class_tau_bounds*1000,
# 						MPP_PlusMinus_InEx_tau_bounds*1000]]

# bound_names_PlusMinusInEx = [[r'$\langle \sigma \rangle$' + ' bound',
# 							r'$\langle \sigma \rangle$' +' MPP PlusMinus InEx bound'],
# 						[r'$\tau$' + ' bound',
# 							r'$\tau$' + ' MPP PlusMinus InEx bound']]

# plot_bounds(times*1000, actuals, actual_names, bounds_InEx, bound_names_InEx, palettes[4])
# plot_bounds(times*1000, actuals, actual_names, bounds_RearrInEx, bound_names_RearrInEx, palettes[4])
# plot_bounds(times*1000, actuals, actual_names, bounds_PlusMinusInEx, bound_names_PlusMinusInEx, palettes[4])

# B_bds = [[class_EP_bounds, 
# 					EPu[uind['omega']][1:] + EPu[uind['alpha']][1:],
# 					EPu2[uind2['beta']][1:] + EPu2[uind2['alpha']][1:],
# 					MPP_B1_EP_bounds,
# 					MPP_B2_EP_bounds],
# 				[class_tau_bounds*1000,
# 					MPP_B1_tau_bounds*1000,
# 					MPP_B2_tau_bounds*1000]]

# B_bd_names = [[r'$\langle \sigma \rangle$' + ' CSL',
# 						r'$\langle \sigma \rangle$' +' EPsum B1',
# 						r'$\langle \sigma \rangle$' +' EPsum B2',
# 						r'$\langle \sigma \rangle$' +' MPP B1',
# 						r'$\langle \sigma \rangle$' +' MPP B2'],
# 					[r'$\tau$' + ' CSL',
# 						r'$\tau$' + ' MPP B1',
# 						r'$\tau$' + ' MPP B2']]

# plot_bounds(times*1000, actuals, actual_names, B_bds, B_bd_names, palettes[4])


# all_bds = [[class_EP_bounds, 
# 					MPP_ReArr_InEx_EP_bounds,
# 					MPP_ReArr_InEx_EP_bounds2,
# 					MPP_PlusMinus_InEx_EP_bounds,
# 					class_EPu_bounds[0]],
# 				[class_tau_bounds*1000,
# 					MPP_ReArr_InEx_tau_bounds*1000,
# 					MPP_ReArr_InEx_tau_bounds2*1000,
# 					MPP_PlusMinus_InEx_tau_bounds*1000,
# 					class_tauu_bounds[0]*1000]]

# all_bd_names = [[r'$\langle \sigma \rangle$' + ' CSL',
# 						r'$\langle \sigma \rangle$' +' MPP R1',
# 						r'$\langle \sigma \rangle$' +' MPP R2',
# 						r'$\langle \sigma \rangle$' +' MPP PM',
# 						r'$\langle \sigma \rangle^{*}_{\phi}$'],
# 					[r'$\tau$' + ' CSL',
# 						r'$\tau$' + ' MPP R1',
# 						r'$\tau$' + ' MPP R2',
# 						r'$\tau$' + ' MPP PM',
# 						r'$\tau^{*}_{\phi}$']]

# plot_bounds(times*1000, actuals, actual_names, all_bds, all_bd_names, palettes[4])

#**************TODOOOO******************************
# #Plot Contributions to Rearranged InEx tau Bounds 
# more tricky because of multiplicities
#**************TODOOOO******************************

# CFT ReArr Bounds
# actuals_CFT = [EPu[uind['alpha']][1:], times[1:]*1000]
# actual_names_CFT = [r'$\langle \sigma^{\alpha} \rangle$', r'$\tau$ ' + '(ms)']
# bounds_CFT = [[MPP_CFT_EP_alphaphi],[MPP_CFT_tau_alphaphi]]
# bound_names_CFT = [['CFT ' + r'$\phi$ ' + r'$\sigma$'],['CFT ' + r'$\phi$ ' + r'$\tau$']]
# plot_bounds(times*1000, actuals_CFT, actual_names_CFT, bounds_CFT, bound_names_CFT, palettes[4])

# # EP Rate Counterfactual Bounds
# actuals_EPR = [EPratesys[1:]]
# actual_names_EPR = [r'$\langle \dot{\sigma} \rangle$']
# bounds_EPR = [[MPP_EPR_ID1, MPP_EPR_ID2]]
# bound_names_EPR = [['ID1', 'ID2']]
# plot_bounds(times*1000, actuals_EPR, actual_names_EPR, bounds_EPR, bound_names_EPR, palettes[4])


# # Plot Evolution of Subsystem Probability Distributions Over Time
# plot_multiple_pdist_evolution(all_ss, all_spds, names, palettes[2])


# *********************************** Sanity Checks ************************************************

# plot_bounds(times*1000, [EPsys[1:]], [r'$\langle \sigma \rangle$'], 
# 			[[InExEP[1:] - DeltaInEx_Info[1:], InExEP2[1:] - DeltaInEx_Info2[1:]]], 
# 			[[r'$\langle \sigma \rangle$' + 'InEx1', r'$\langle \sigma \rangle$' + 'InEx2']], palettes[2])


# bounds = 	[[class_EP_bounds, 
# 					MPP_ID1a_EP_bounds, 
# 					MPP_ID1b_EP_bounds,
# 					MPP_ID2a_EP_bounds, 
# 					MPP_ID2b_EP_bounds,
# 					MPP_ID2c_EP_bounds, 
# 					MPP_InExSum_EP_bounds,
# 					MPP_ReArr_InEx_EP_bounds,
# 					MPP_InExSum_EP_bounds2,
# 					MPP_ReArr_InEx_EP_bounds2],
# 				[class_tau_bounds*1000, 
# 					MPP_ID1a_tau_bounds*1000, 
# 					MPP_ID1b_tau_bounds*1000,
# 					MPP_ID2a_tau_bounds*1000, 
# 					MPP_ID2b_tau_bounds*1000,
# 					MPP_ID2c_tau_bounds*1000,  
# 					MPP_InExSum_tau_bounds*1000,
# 					MPP_ReArr_InEx_tau_bounds*1000,
# 					MPP_InExSum_tau_bounds2*1000,
# 					MPP_ReArr_InEx_tau_bounds2*1000]]

# bound_names = 	[[r'$\langle \sigma \rangle$' + ' bound', 
# 						r'$\langle \sigma \rangle$' +' MPP ID1a bound',
# 						r'$\langle \sigma \rangle$' +' MPP ID1b bound',
# 						r'$\langle \sigma \rangle$' +' MPP ID2a bound',
# 						r'$\langle \sigma \rangle$' +' MPP ID2b bound',
# 						r'$\langle \sigma \rangle$' +' MPP ID2c bound',
# 						r'$\langle \sigma \rangle$' +' MPP InEx bound',
# 						r'$\langle \sigma \rangle$' +' MPP ReArr InEx bound',
# 						r'$\langle \sigma \rangle$' +' MPP InEx2 bound',
# 						r'$\langle \sigma \rangle$' +' MPP ReArr InEx2 bound'],
# 					[r'$\tau$' + ' bound', 
# 						r'$\tau$' + ' MPP ID1a bound',
# 						r'$\tau$' + ' MPP ID1b bound',
# 						r'$\tau$' + ' MPP ID2a bound',
# 						r'$\tau$' + ' MPP ID2b bound',
# 						r'$\tau$' + ' MPP ID2c bound',
# 						r'$\tau$' + ' MPP InEx bound',
# 						r'$\tau$' + ' MPP ReArr InEx bound',
# 						r'$\tau$' + ' MPP InEx2 bound',
# 						r'$\tau$' + ' MPP ReArr InEx2 bound']]


# bounds_ID1 = 	[[class_EP_bounds, 
# 					MPP_ID1a_EP_bounds, 
# 					MPP_ID1b_EP_bounds],
# 				[class_tau_bounds*1000, 
# 					MPP_ID1a_tau_bounds*1000, 
# 					MPP_ID1b_tau_bounds*1000]]

# bound_names_ID1 = 	[[r'$\langle \sigma \rangle$' + ' bound', 
# 						r'$\langle \sigma \rangle$' +' MPP ID1a bound',
# 						r'$\langle \sigma \rangle$' +' MPP ID1b bound'],
# 					[r'$\tau$' + ' bound', 
# 						r'$\tau$' + ' MPP ID1a bound',
# 						r'$\tau$' + ' MPP ID1b bound']]

# bounds_ID2 = 	[[class_EP_bounds, 
# 					MPP_ID2a_EP_bounds, 
# 					MPP_ID2b_EP_bounds,
# 					MPP_ID2c_EP_bounds],
# 				[class_tau_bounds*1000, 
# 					MPP_ID2a_tau_bounds*1000, 
# 					MPP_ID2b_tau_bounds*1000,
# 					MPP_ID2b_tau_bounds*1000]]

# bound_names_ID2 = 	[[r'$\langle \sigma \rangle$' + ' bound', 
# 						r'$\langle \sigma \rangle$' +' MPP ID2a bound',
# 						r'$\langle \sigma \rangle$' +' MPP ID2b bound',
# 						r'$\langle \sigma \rangle$' +' MPP ID2c bound'],
# 					[r'$\tau$' + ' bound', 
# 						r'$\tau$' + ' MPP ID2a bound',
# 						r'$\tau$' + ' MPP ID2b bound',
# 						r'$\tau$' + ' MPP ID2c bound']]

# ID1a_labels = [r"$\phi$", r"$\omega \backslash \phi$", r"$\mathcal{N} \backslash \omega$"]
# ID1b_labels = [r"$\phi$", r"$\alpha \backslash \phi$", r"$\mathcal{N} \backslash \alpha$"]

# ID2a_labels = [r'$\beta$', r'$\omega \backslash \beta$', r'$\mathcal{N} \backslash \omega$']
# ID2b_labels = [r'$\phi$', r'$\beta \backslash \phi$', r'$\omega \backslash \beta$', r'$\mathcal{N} \backslash \omega$']
# ID2c_labels = [r'$\phi$', r'$\beta \backslash \phi$', r'$\mathcal{N} \backslash \beta$']

# # Plot Contributions to Iterative Decomp Bounds Under Unit Structure 1
# plot_contribs_to_bound(times*1000, MPP_ID1a_EP_contribs, ID1a_labels, 'ID1a EP', palettes[2])
# plot_contribs_to_bound(times*1000, MPP_ID1a_tau_contribs*1000, ID1a_labels, 'ID1a tau', palettes[2])
# plot_contribs_to_bound(times*1000, MPP_ID1b_EP_contribs, ID1b_labels, 'ID1b EP', palettes[2])
# plot_contribs_to_bound(times*1000, MPP_ID1b_tau_contribs*1000, ID1b_labels, 'ID1b tau', palettes[2])

# # Plot Contributions to Iterative Decomp Bounds Under Unit Structure 2
# plot_contribs_to_bound(times*1000, MPP_ID2a_EP_contribs, ID2a_labels, 'ID2a EP', palettes[2])
# plot_contribs_to_bound(times*1000, MPP_ID2a_tau_contribs*1000, ID2a_labels, 'ID2a tau', palettes[2])
# plot_contribs_to_bound(times*1000, MPP_ID2b_EP_contribs, ID2b_labels, 'ID2b EP', palettes[2])
# plot_contribs_to_bound(times*1000, MPP_ID2b_tau_contribs*1000, ID2b_labels, 'ID2b tau', palettes[2])
# plot_contribs_to_bound(times*1000, MPP_ID2c_EP_contribs, ID2c_labels, 'ID2c EP', palettes[2])
# plot_contribs_to_bound(times*1000, MPP_ID2c_tau_contribs*1000, ID2c_labels, 'ID2c tau', palettes[2])
