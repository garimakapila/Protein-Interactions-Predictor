# Copyright © Garima Kapila

import matplotlib.pyplot as plt
import pandas as pd, pptx, re, scipy, seaborn, StringIO
seaborn.set(color_codes = True)
seaborn.set_color_codes('pastel')


def column_comparisions(dataframe0, dataframe1, file_name, p_values):
	graphs = []
	p_values_dict = dict(zip(p_values['Column'], p_values['P-Value']))
	for col in p_values['Column']:
		x0, x1 = dataframe0[col], dataframe1[col]
		plt.clf()
		seaborn.violinplot(data = [[], x0, x1])
		plt.boxplot([x0, x1])
		p_value = p_values_dict[col]
		graphs.append(graph_as_image(title = col + ' ' + str(p_value)))
	graph_file_name = 'Graphs_' + file_name[:-4]
	write_images_to_ppt(graphs, graph_file_name)

def line_graph(data, title, x_axis, y_axis):
	label_plot(title, x_axis, y_axis)
	plt.plot(data)
	return graph_as_image()



##################
# Helper Methods #
##################

def write_images_to_ppt(images, title, file=''):
	if file != '':
		ppt = pptx.Presentation(file)
		name = file
	else:
		ppt = pptx.Presentation()
		name = title + '.pptx'
	for image in images:
		slide = ppt.slides.add_slide(ppt.slide_layouts[6])
		inch = pptx.util.Inches(1)
		slide = slide.shapes.add_picture(image, inch, inch) 
	ppt.save(name)

def graph_as_image(title = ''):
	if title != '':
		plt.title(title.replace('_', ' '))
	fig = plt.gcf()
	imgdata = StringIO.StringIO()
	fig.savefig(imgdata, format = 'png')
	imgdata.seek(0)
	return imgdata

def label_plot(title, x_axis, y_axis):
	plt.clf()
	plt.title(title)
	plt.xlabel(x_axis)
	plt.ylabel(y_axis)

"""
Splits data1 into nbins # bins, plot against data2
Plots boxplot-violin plot for 0, 1 for each bucket
"""
def bin_graph(data1, data2, labels, nbins):
	bins, ranges = get_bins_and_ranges(data1, data2, labels, nbins)
	split_bins = []
	for sub_bin in bins:
		sub_bin0, sub_bin1 = split_by_label(sub_bin)
		split_bins.append(sub_bin0)
		split_bins.append(sub_bin1)
	plt.clf()
	seaborn.violinplot(data = [[]] + split_bins)
	x_label = data1.name.replace('_', ' ')
	plt.xlabel(x_label + ', Bins: '+ format_ranges(ranges))
	y_label = data2.name.replace('_', ' ')
	plt.ylabel(y_label)
	plt.boxplot(split_bins)
	plt.xticks([])
	plt.show()

def format_ranges(ranges):
	formatted = ''
	for start, end in ranges:
		formatted += '[' + str(start) + ', ' + str(end) + '], '
	return formatted[:-2]

def get_bins_and_ranges(data1, data2, labels, nbins):
	ranges = get_ranges(data1, nbins)
	bins, new_ranges = [], []
	for start, end in ranges:
		bin_range = get_in_range(data1, data2, labels, start, end)
		bins.insert(0, bin_range)
		new_ranges.insert(0, [int(start), int(end)])
	return bins, new_ranges

def get_ranges(data, nbins):
	min_value = float(data.min())
	max_value = float(data.max())
	bin_size = float(max_value - min_value)/nbins
	index = min_value
	ranges = []
	while index < max_value - float(bin_size)/2:
		new_index = index + bin_size
		ranges.append((index, new_index))
		index = new_index
	return ranges

# inclusive
def get_in_range(data1, data2, labels, start, end):
	values1 = data1.values.tolist()
	values2 = data2.values.tolist()
	labels_list = labels.values.tolist()
	ranges = []
	for value1, value2, label in zip(values1, values2, labels_list):
		if value1 >= start and value1 <= end:
			ranges.append([value2, label])
	return pd.DataFrame(ranges, columns=[data2.name, labels.name])

def get_p_values(dataframe0, dataframe1, file_name):
	p_values = []
	for col in dataframe0.columns:
		x0, x1 = dataframe0[col], dataframe1[col]
		pvalue = scipy.stats.mannwhitneyu(x0, x1)[1]
		p_values.append([col, pvalue])
	p_values = pd.DataFrame(p_values, columns = ['Column', 'P-Value'])
	p_values = p_values.sort_values(by = ['P-Value'])
	p_values.to_csv(file_name, sep = ',', index = False)
	return p_values
