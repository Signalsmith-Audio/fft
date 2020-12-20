import os

import matplotlib
from matplotlib import pyplot

def small(*args, **kwargs):
	figure, axes = pyplot.subplots(*args, **kwargs)
	figure.set_size_inches(4.5, 3)
	return figure, axes

def medium(*args, **kwargs):
	figure, axes = pyplot.subplots(*args, **kwargs)
	figure.set_size_inches(6.5, 4)
	return figure, axes

def tall(*args, **kwargs):
	figure, axes = pyplot.subplots(*args, **kwargs)
	figure.set_size_inches(4.5, 5.5)
	return figure, axes

def short(*args, **kwargs):
	figure, axes = pyplot.subplots(*args, **kwargs)
	figure.set_size_inches(7, 3)
	return figure, axes

def wide(*args, **kwargs):
	figure, axes = pyplot.subplots(*args, **kwargs)
	figure.set_size_inches(11, 4)
	return figure, axes

def full(*args, **kwargs):
	figure, axes = pyplot.subplots(*args, **kwargs)
	figure.set_size_inches(16, 10)
	return figure, axes

def save(prefix, figure, legend_loc=0):
	dirname = os.path.dirname(prefix)
	if len(dirname) and not os.path.exists(dirname):
		os.makedirs(dirname)

	if len(figure.get_axes()) > 1:
		figure.set_tight_layout(True)
	for axes in figure.get_axes():
		if len(axes.get_lines()) > 0 and legend_loc != None:
			axes.legend(loc=legend_loc)

	print(prefix)
	figure.savefig(prefix + '.svg', bbox_inches='tight')
