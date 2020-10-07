import os.path
import csv
import numpy
import matplotlib
from matplotlib import pyplot

import article

def plot_if_exists(axis, csvFile):
	axis.set(xlabel="size", ylabel="normalised rate")
	axis.set_xscale('log', basex=4)
	axis.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, pos: ("%i"%x) if (x <= 65536) else ("2^%i"%int(numpy.log2(x)))))
	if not os.path.exists(csvFile):
		return False
	with open(csvFile) as inputCsv:
		reader = csv.reader(inputCsv)
		header = next(reader)
		(_, _, lineName) = header

		sizes = []
		rates = []
		for row in reader:
			(size, rate, scaledRate) = [float(x) for x in row]
			sizes.append(size)
			rates.append(scaledRate)
		axis.plot(sizes, rates, label=lineName)

figure, axis = article.wide()
plot_if_exists(axis, "out/results/signalsmith.csv")
plot_if_exists(axis, "out/results/fftw-estimate.csv")
plot_if_exists(axis, "out/results/fftw-measure.csv")
plot_if_exists(axis, "out/results/kissfft.csv")
article.save("out/comparison", figure)

figure, axis = article.wide()
plot_if_exists(axis, "out/results/signalsmith.csv")
plot_if_exists(axis, "out/results/previous-v4.csv")
plot_if_exists(axis, "out/results/previous-permute.csv")
plot_if_exists(axis, "out/results/previous.csv")
article.save("out/previous", figure)

figure, axis = article.wide()
plot_if_exists(axis, "out/results/signalsmith.csv")
plot_if_exists(axis, "out/results/dev-history-radix23.csv")
plot_if_exists(axis, "out/results/dev-history-factorise.csv")
plot_if_exists(axis, "out/results/dev-history-direct.csv")
article.save("out/history", figure)
