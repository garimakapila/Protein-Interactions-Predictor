from scipy.optimize import curve_fit
import numpy as np
import pandas as pd
from sklearn import preprocessing
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from DataCollector.DataParser import *
from Visualizer import *



# keep/delete certain columns
def filter_cols(features_file, include=[], exclude=[]):
	features = pd.read_csv(features_file)
	cols = [
		col for col in features.columns if
			(any(keyword in col for keyword in include) and
			 all(keyword not in col for keyword in exclude))
	]
	results = features[cols]
	results_file = features_file.replace('Features', 'Reduced')
	results.to_csv(results_file, index = False)
	return results_file

# standardize data then PCA
def dimensionality_reduction(file_name, scree_file_name='', graph=False):
	data, labels = columns_from_data(file_name)
	std_data = standardize_data(data)
	n_cols = len(data.columns)
	n, scree, log_scree = get_n_components(std_data, n_cols)
	scree_plots_name = file_name.replace('Features', 'Scree_Plots')
	write_images_to_ppt([scree, log_scree], scree_plots_name, scree_file_name)
	pca = PCA(n_components = n).fit(std_data)
	reduced = pca.transform(std_data)
	if graph:
		scatter_plots(file_name, reduced)
	reduced_file = 'Reduced_' + file_name
	np.savetxt(reduced_file, reduced, delimiter = ',')
	return n, reduced_file, scree_file_name

# scatter points using PCA and TSNE
def scatter_plots(data, labels):
	scatter_graph(data, labels)

# scale features by standardization
def standardize_data(data):
	scaler = preprocessing.StandardScaler()
	std_scale = scaler.fit(data)
	std_data = std_scale.transform(data)
	return std_data

# get top principal components (# at kink in scree plot)
def get_n_components(data, n_cols):
	pca = PCA(n_components = n_cols).fit(data)
	data = pca.transform(data)
	variance = pca.explained_variance_
	log_variance = np.log(variance)
	scree, log_scree = scree_plot(variance, log_variance)
	n = int((np.abs(log_variance)).argmin() * 1.5)
	return n, scree, log_scree

# visualize kink of explained variance plot
def scree_plot(variance, log_variance):
	title = 'Scree Plot'
	x_axis = 'Principal Component #'
	y_axis = 'Explained Variance'
	scree = line_graph(variance, title, x_axis, y_axis)
	title = 'Log Scree Plot'
	y_axis = 'Log Explained Variance'
	log_scree = line_graph(log_variance, title, x_axis, y_axis)
	return scree, log_scree

# get number columns and labels
def columns_from_data(file_name):
	data = pd.read_csv(file_name, sep = ',')
	data = data.select_dtypes(['number'])
	labels = data['Label']
	data = data.drop(['Label'], axis = 1)
	return data, labels

# For scatter plots

def pca_2D(reduced_pca_data):
	return reduced_pca_data[:, 0], reduced_pca_data[:, 1]

def tsne_2D(reduced_data):
	return TSNE(n_components = 2).fit_transform(reduced_data)
	

def combine_features(file_path):
	features = pd.read_csv(file_path, sep=',')
	features = add_blosum_positives(features)
	features = add_score_averages(features)
	features = join_pair_columns(features)
	features = move_col_to_end(features, 'Database')
	features = move_col_to_end(features, 'Label')
	file_name = file_path.split('/')[-1]
	file_name_prefix = file_name.split('_')[0]
	file_name = file_name.replace(file_name_prefix, 'Features')
	features.fillna(value = 0)
	features.to_csv(file_name, sep = ',', index = False)
	return file_name


# input pandas dataframe, dataframe split by labels
def split_by_label(dataframe, labels_col, labels=[0, 1]):
	label_name = labels_col.name
	dataframe[label_name] = labels_col
	splits = [dataframe[dataframe[label_name] == label] for label in labels]
	return [split.drop(columns=[label_name]) for split in splits]


### Helper methods for add_features_to_interologs ###

# Add blosum match + similar (non-match) to get positives column
def add_blosum_positives(features):
	for col in features.columns.values:
		if '_Match' in col and '_Blosum' in col:
			col_type = column_type(col)
			pr_type = pair_type(col)
			for col2 in features.columns.values:
				if '_Similar' in col2 and pr_type == pair_type(col2) and      \
					col_type == column_type(col2):
					features[col_type + '_Sum_Blosum_Positive_Scores_Pair_'   \
					+ pr_type] = features[col] + features[col2]
	return features

# average # or score divided by length
def add_score_averages(features):
	
	# Put column names that have 'length' in a dictionary along with pair type
	lengths_dict = {}
	for col in features.columns.values:
		if 'Length' in col:
			col_type = column_type(col)
			pr_type = pair_type(col)
			lengths_dict[col_type + pr_type] = col
			# avoid infinite values
			features[col] = features[col].replace(0, 1)

	# Get averages/percents using lengths
	for col in features.columns.values:
		if 'Sum' in col or 'Count' in col or 'Difference' in col:
			length_col = column_type(col) + pair_type(col)
			if length_col in lengths_dict.keys():
				lengths = features[lengths_dict[length_col]]
				factor = 1
				if 'Count' in col:
					factor = 100
				new_col = col.replace('Sum', 'Average_Sum')
				new_col = new_col.replace('Count', 'Percent')
				new_col = new_col.replace('Difference', 'Average_Difference')
				features[new_col] = features[col] * factor / lengths
				features[new_col] = features[new_col].round(3)
	
	# Add % interface residues
	for col in features.columns.values:
		if 'Interface_' in col and 'Length' in col:
			pr_type = pair_type(col)
			fasta_length = lengths_dict['Fasta' + pr_type]
			new_col = col.replace('Length', 'Percent')
			features[new_col] = 100 * features[col] / features[fasta_length]
			features[new_col].round(3)
	
	return features

# Get joint columns for pair 1, pair 2 across features
def join_pair_columns(features):
	for col1 in features.columns:
		if 'Pair_1' in col1:
			column_name = col1.replace('Pair_1', 'Joint')
			data1 = features[col1]
			col2 = col1[:-1] + '2'
			data2 = features[col2]
			joint = data1 * data2
			features[column_name] = joint
	return features

# first part of name is column type, types: Blast, Global, Interface, etc
def column_type(column):
	return column.split('_')[0]

# last part of name is pair type, ex. Pair_1 = 1
def pair_type(column):
	return column.split('_')[-1]

def filter_interologs(features_file):
	data = pd.read_csv(features_file)
	data0 = data[data['Label'] == 0]
	data1 = data[data['Label'] == 1]
	q0_0, q0_1, q1_0, q1_1 = 0, 0, 0, 0
	for x in range(0, 10):
		x /= float(100)
		percent_m1, percent_m2 = [
			'Interface_Matching_Percent_Pair_1',
			'Interface_Matching_Percent_Pair_2',
		]
		x0_0 = data0[percent_m1].quantile(x)
		x0_1 = data0[percent_m2].quantile(x)
		if x0_0 == 0 and x0_1 == 0:
			q0_0 = x0_0
			q0_1 = x0_1
			q1_0 = data1[percent_m1].quantile(x)
			q1_1 = data1[percent_m2].quantile(x)
		else:
			break
