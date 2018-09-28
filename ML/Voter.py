# Copyright © Garima Kapila

from numpy import genfromtxt
import pandas as pd
import time
from DataParser import *
from sklearn.metrics import confusion_matrix
import pandas as pd#, tensorflow as tf
from sklearn.model_selection import train_test_split
from sklearn.decomposition import PCA

class Voter():

	def __init__(self, supervised, unsupervised, x, y, n):
		self.supervised = supervised
		self.unsupervised = unsupervised
		self.x = np.nan_to_num(x)
		self.y = np.nan_to_num(y)
		self.n = n

	# take mode over several runs
	def get_predictions(self):
		predictions = []
		for i in range(0, 10):
			pred = self.get_votes()
			predictions.append(pred)
		predictions = pd.DataFrame(predictions)
		predictions = predictions.mode(axis = 1)
		return predictions

	def get_votes(self):
		# Get seed data
		seed_x, seed_y = self.get_seed() #self.get_seed_one_to_one()
		x_train, x_test, y_train, y_test = train_test_split(
			seed_x, seed_y, test_size = 0.33)

		# Supervised
		self.supervised.fit(x_train, y_train)
		predictions_supervised = self.supervised.predict(self.x)
		self.print_accuracy(predictions_supervised)

		# Unsupervised
		pca = PCA(n_components = 5).fit(self.x)
		x = pca.transform(self.x)
		self.unsupervised.fit(x)
		predictions_unsupervised = self.unsupervised.fit_predict(x)
		self.print_accuracy(predictions_unsupervised)

		label_counts = [{} for i in range(self.n)]
		labels = {}
		for i in range(predictions_supervised.shape[0]):
			y = predictions_unsupervised[i]
			if predictions_supervised[i] not in label_counts[y]:
			  label_counts[y][predictions_supervised[i]] = 0
			label_counts[y][predictions_supervised[i]] += 1

		# Label cluster based on majority
		sketch_level = 0
		for i in range(len(label_counts)):
			count = label_counts[i]
			if count:
			  labels[i] = max(count, key = count.get)
			else:
			  labels[i] = 0
			  sketch_level += 1
		predictions_unsupervised = np.vectorize(labels.get)                   \
			(predictions_unsupervised)

		return predictions_supervised

	def get_seed(self):
		seed_x, seed_y = [], []
		random = int(time.time())
		for i, points in enumerate(zip(self.x, self.y)):
			x_point, y_point = points
			if (i + random) % 20 == 0:
				seed_x.append(x_point)
				seed_y.append(y_point)
		return seed_x, seed_y

	def get_seed_one_to_one(self):
		seed_x, seed_y = [], []
		count0, count1 = 0, 0
		random = int(time.time())
		for i, points in enumerate(zip(self.x, self.y)):
			x_point, y_point = points
			if (i + random) % 20 == 0:
				if y_point == 1:
					count1 += 1
					seed_x.append(x_point)
					seed_y.append(1)
				elif count0 < count1:
					seed_x.append(x_point)
					seed_y.append(0)
		return seed_x, seed_y

	def print_accuracy(self, predictions):
		categories = [0, 1]
		for category in categories:
			count = 0
			count_total = 0
			for label, pred in zip(self.y, predictions):
				if label == category and pred == category:
					count += 1
				if label == category:
					count_total += 1
			print count, count_total
			
