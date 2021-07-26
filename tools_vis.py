from matplotlib import pyplot as plt
import seaborn as sns

def plot_bounds(ts, actual, actual_name, bounds, bound_names, p, name, size):
	'''
		inputs: 
			ts: (nparray) of times over which to plot,
			actuals: (list) of quantities to plot, each of which is a list of the quantity's value
				at every time in ts,
			actual_names: (list) of names (strings) labeling actuals,
			bounds: (list) of lists of bounds to plot, with bounds[i] a list of bounds for 
				actuals[i] and bounds[i][j] the jth type of bound for actuals[i],
			bound_names: (list) of list of strings labeling the bounds,
			p: (list) of colors to cycle through when plotting each subplot of actual value vs.
				bounds for that value,
			name: name to save the figure as,
			size: desired size of the figure
		outputs: 
			None. This functions plots bounds versus actual values.
	'''

	sns.set_style("whitegrid", {"axes.linecolor": ".9"})
	sns.set_context("notebook", rc={"lines.linewidth": 4})

	num_subplots = len(bounds)
	#create the entire figure with one subplot for every type of actual vs. bounds to plot
	f, axs = plt.subplots(1, num_subplots, sharey=True, figsize=size) 

	if num_subplots > 1:

		for i, ax in enumerate(axs): # populate a subplot
			ax.plot(ts[1:], actual, label = actual_name, color=p[i][0])
			for j, bound in enumerate(bounds[i]):
				ax.plot(ts[1:], bound, label = bound_names[i][j], color=p[i][j+1])

			# ax.set_xlabel('time (ms)')
			# ax.set_ylabel(actual_names[i])
			# ax.legend()

	else:

		axs.plot(ts[1:], actual, label = actual_name)
		for j, bound in enumerate(bounds[0]):
			axs.plot(ts[1:], bound, label = bound_names[0][j])

		# axs.set_xlabel('time (ms)')
		# axs.set_ylabel(actual_names[0])
		# axs.legend()

	# plt.tight_layout()
	plt.savefig(name + '.png', dpi=600)
	# plt.show()

def plot_contribs_to_bound(ts, contribs, lbls, ttle, p):
	'''
		inputs:
		outputs:
	'''
	sns.set('talk', style='whitegrid', palette=sns.color_palette(p, len(contribs)))
	f, axs = plt.subplots(1, 1, sharey=False)
	axs.stackplot(ts[1:], contribs, labels=lbls)
	axs.set_xlabel('time (ms)')
	# axs.set_ylabel(y)
	axs.set_title(ttle)
	axs.legend(loc='upper left')
	plt.show()

def plot_prob_dist_evolution(ss, pds_t, p):
	'''
		inputs:
		outputs:
	'''
	sns.set('talk', style='white' , palette=sns.color_palette(p, len(pds_t)))
	for pd in pds_t:
		plt.plot(ss, pd)
	plt.tight_layout()
	plt.show()

def plot_multiple_pdist_evolution(all_ss, all_pds_t, all_names, p):
	'''
		inputs:
		outputs:
	'''
	sns.set('talk', style='white' , palette=sns.color_palette(p, len(all_pds_t[0])))
	f, axs = plt.subplots(1, len(all_ss), sharey=False) 
	# f.set_figwidth(4) 
	# f.set_figheight(1) 

	for i, ax in enumerate(axs): # populate a subplot
		for pd in all_pds_t[i]:
			ax.plot(all_ss[i], pd)

		ax.set_xticks(all_ss[i])
		ax.set_xlabel(all_names[i])
		ax.set_ylabel('p(' + all_names[i] + ')')

	# plt.tight_layout()
	plt.show()