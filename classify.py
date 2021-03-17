import os
import sys
import argparse
from argparse import RawTextHelpFormatter
from statistics import mode

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


def xsmooth(data):
	out = np.zeros_like(data)
	var = np.array([
					data[0::6], 
					data[1::6],
					data[2::6],
					data[3::6],
					data[4::6],
					data[5::6]
					])
	for i in range(var.shape[1]):
		#print(var[:,max(i-2, 0) : i+3])
		counts = np.count_nonzero(var[:,max(i-4, 0) : i+5] == 2, axis=1)
		idx = np.argmax(counts)
		#idx = np.argmax(np.count_nonzero(var[:,max(i-2, 0) : i+3] == 2, axis=1))
		if counts[idx] >= 3:
			out[6*i+idx] = 2
	return out

if __name__ == '__main__':
	usage = '%s [-opt1, [-opt2, ...]] infile' % __file__
	parser = argparse.ArgumentParser(description='', formatter_class=RawTextHelpFormatter, usage=usage)
	parser.add_argument('infile', type=is_valid_file, help='input file')
	parser.add_argument('-o', '--outfile', action="store", default=sys.stdout, type=argparse.FileType('w'), help='where to write the output [stdout]')
	args = parser.parse_args()

	train = pd.read_csv(sys.argv[1], header=None, sep='\t')
	X = train.iloc[:,2:26]
	Y = train.iloc[:,1].replace({'None':0, 'False':1, 'True':2})

	'''
	cp_callback = tf.keras.callbacks.ModelCheckpoint(filepath='cp.ckpt', save_weights_only=True, verbose=1)
	model = create_model('adam')
	model.fit(X, Y, epochs=11, callbacks=[cp_callback])
	'''

	model = create_model('adam')
	model.load_weights('cp.ckpt')

	p = model.predict(X)
	f = np.argmax(p,axis=-1)
	#smooth(f)
	ff = xsmooth(f)
	
	for row in zip(train.iloc[:,0].to_list(), ff):
		if row[1] == 2:
			if row[0] > 0:
				print('     CDS             ', row[0] , '..', row[0]+2, sep='')
			else:
				print('     CDS             complement(', abs(row[0]), '..', abs(row[0])+2, ')', sep='')



