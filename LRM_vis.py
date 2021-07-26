from tools_vis import *
from LRM_calc import *
from LRM_pdist import *

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Source Sans Pro']})

actual = times_LRM[1:]*1000
actual_name = r'$\tau$ ' + '(ms)'

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

plot_bounds(times_LRM*1000, actual, actual_name, second_bounds_to_plot, second_bound_names_to_plot, second_bound_colors_to_plot, 'units', (25,4))