import os
import sys
import argparse
from argparse import RawTextHelpFormatter

# TensorFlow and tf.keras
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import tensorflow as tf

# Helper libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def is_valid_file(x):
	if not os.path.exists(x):
		raise argparse.ArgumentTypeError("{0} does not exist".format(x))
	return x

def create_model(opt):
	'''
	This creates and returns a new model
	'''
	model = tf.keras.Sequential([
					tf.keras.layers.Dense(24, input_shape=(24,)),
					tf.keras.layers.Dense(128, activation='relu'),
					tf.keras.layers.Dense(128, activation='relu'),
					tf.keras.layers.Dense(10, activation='softmax')
	])
	model.compile(optimizer = opt,
				  loss=tf.keras.losses.SparseCategoricalCrossentropy(from_logits=True),
				  metrics=['accuracy']
				  )
	return model

if __name__ == '__main__':
	usage = '%s [-opt1, [-opt2, ...]] infile' % __file__
	parser = argparse.ArgumentParser(description='', formatter_class=RawTextHelpFormatter, usage=usage)
	parser.add_argument('infile', type=is_valid_file, help='input file')
	parser.add_argument('-o', '--outfile', action="store", default=sys.stdout, type=argparse.FileType('w'), help='where to write the output [stdout]')
	args = parser.parse_args()

	train = pd.read_csv(sys.argv[1], header=None, sep='\t')
	X = train.iloc[:,2:26]
	Y = train.iloc[:,1].replace({'None':0, 'False':1, 'True':2})

	model = create_model('adam')
	model.fit(X, Y, epochs=10)




