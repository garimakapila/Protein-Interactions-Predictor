from ML.DataParser import *
from ML.Visualizer import *
from ML.Voter import *
from Organism import Organism
import numpy as np, os, shutil, pandas as pd
from sklearn.svm import LinearSVC
from sklearn.mixture import GaussianMixture
from sklearn.ensemble import RandomForestClassifier
from sklearn.cluster import AgglomerativeClustering

class Results():

	# input organism is of type Organism
	def __init__(self, organism1, organism2, option):
		self.organism1 = organism1
		self.organism2 = organism2
		self.option = option
		self.names = organism1.name + '_' + organism2.name

		# directory to move collected files to
		self.directory = self.option + '/Results/' + organism2.name

		# make directories
		if not os.path.exists(self.directory):
			os.makedirs(self.directory)


	def run(self):
		# 1. Combine raw data features and generate graphs
		reduced_file, labels = self.analyze_data()	
		# 2. Make predictions
		self.vote(reduced_file, labels)
		# 3. Move results files into 'Results' directory
		self.clean()


	def analyze_data(self, graph=False):
		interologs_file_path = self.option + '/Data/' + self.organism2.name + \
			'/' + self.organism1.name + '/Interologs_' + self.names + '.csv'
		features_file = combine_features(interologs_file_path)

		data, labels = columns_from_data(features_file)
		data0, data1 = split_by_label(data, labels)

		p_values_file = features_file.replace('Features', 'P-Values')
		p_values = get_p_values(data0, data1, p_values_file)
		
		if graph:
			column_comparisions(data0, data1, features_file, p_values)
		
		#dimensionality_reduction(features_file)
		#n, reduced_file, scree_file = dimensionality_reduction(features_file)
		#cols = reduce_features_by_n(n, p_values_file, features_file)

		reduced_file = filter_cols(features_file, include=['Interface'],
			exclude=['Difference', 'Database', 'Overlapping', 'Gap_Score'])
		filter_interologs(features_file)
		self.standardize(reduced_file)

		#return cols, reduced_file, labels
		return reduced_file, labels


	def write_predictions(self, features_file, predictions):
		features = pd.read_csv(features_file)
		features = features[['A1', 'B1', 'A2', 'B2', 'Label']]
		features['Prediction'] = predictions
		predictions_file = features_file.replace('Features', 'Predictions')
		features.to_csv('Predictions', index = 1)


	def vote(self, reduced_file, labels):
		SVM = LinearSVC(penalty = 'l2', dual = False)
		RFC = RandomForestClassifier(max_depth = 2, random_state = 0)
		n = 5 # number clusters
		GMM = GaussianMixture(n_components = n)
		Ward = AgglomerativeClustering(linkage = 'ward', n_clusters = n)
		supervised_models = [SVM, RFC]
		unsupervised_models = [Ward, GMM]
 
		x = genfromtxt(reduced_file, delimiter = ',')
		y = labels
		predictions = []
		for supervised in supervised_models:
			for unsupervised in unsupervised_models:
				voter = Voter(supervised, unsupervised, x, y, n)
				predictions.append(voter.get_votes())
				print predictions
		predictions = pd.DataFrame(predictions)
		predictions = predictions.mode(axis = 1)
		return predictions


	def standardize(self, results_file):
		results = pd.read_csv(results_file)
		std_data = standardize_data(results)
		std_data = pd.DataFrame(std_data, columns = results.columns)
		std_data.to_csv(results_file, index = False)


	def clean(self):
		move_keywords = ['.csv', '.pptx']
		for file in os.listdir('.'):
			if 'Features_' in file:
				os.remove(file)
			if any(word in file for word in move_keywords):
				shutil.move(file, self.directory + '/' + file)


	def reduce_features_by_n(self, n, p_values_file, features_file):
		p_values = pd.read_csv(p_values_file, sep=',')
		p_values = p_values.tail(n - int(n * 1.5))
		features = pd.read_csv(features_file, sep=',')
		cols = [col for col in p_values['Column']]
		features.drop(columns = cols, axis = 1, inplace = True)
		features.to_csv(features_file, sep = ',', index = False)
		return [col for col in p_values.head(n)]


