from itertools import combinations

def inter(l1, l2): 
	'''
		inputs:
			l1: (list) of subsystems,
			l2: (list) of subsystems
		outputs:
			(list) of subsystems at the intersection of l1 and l2
	'''
	return list(set(l1) & set(l2)) 

def intersection(list_of_units): 
	'''
		inputs:
			list_of_units: (list) of units, each a list itself
		outputs:
			intersect: (list) of subsystems representing the intersection of those units
				note that by definition this must also be a unit in the unit structure
	'''
	l = len(list_of_units)
	intersect = list_of_units[0]
	if l > 1:
		intersect = inter(list_of_units[0], list_of_units[1])
		for i in range(l):
			if i > 1:
				intersect = inter(intersect, list_of_units[i])
	return intersect

def ordered_unit(unit, suborder):
	ordered_indices = [suborder[u] for u in unit]
	return [i for i, u in suborder.items() if u in ordered_indices]

def get_odd_even_intersections(unit_struct, suborder):
	'''
		inputs:
			unit_struct: (list) of units in the considered unit structure
		outputs:
			intersections: (list) of lists of even and odd intersection units, referenced by their
				index in unit_struct
	'''
	intersections = [[],[]]

	all_combos = [[],[]]
	l = len(unit_struct)
	for i in range(l):

		#generate combos of units without i units in it:
		k = l - i
		kcombos = combinations(unit_struct, k) 

		#append them respectively to the group of even or odd combos
		for combo in kcombos:
			if k % 2 == 0:
				all_combos[1].append(list(combo))
			else:
				all_combos[0].append(list(combo))

	intersect = []
	for i in range(2):
		for combo in all_combos[i]:
			intersect = intersection(combo)
			if len(intersect) > 0:
				intersections[i].append(unit_struct.index(ordered_unit(intersect, suborder)))

	return intersections

def get_plus_minus(oddeven, unit_struct):
	'''
		inputs:
			oddeven: (list) of lists of even and odd intersection units, referenced by their
				index in unit_struct
		outputs:
			(list) of lists of units that give positive and negative terms upon evaluating the 
				in-ex sum; each unit is referenced by its index in unit_struct
	'''
	odds, evens = oddeven
	plus, minus = [], []
	for u, unit in enumerate(unit_struct):
		i = sum([u == o for o in odds]) - sum([u == e for e in evens])
		if i > 0:
			for j in range(i):
				plus.append(u)
		elif i < 0:
			for j in range(abs(i)):
				minus.append(u)
	return [plus, minus]