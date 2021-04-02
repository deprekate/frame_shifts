import os
import sys
import argparse
from argparse import RawTextHelpFormatter
from statistics import mode

import make_train as mt

# TensorFlow and tf.keras
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import tensorflow as tf

# Helper libraries
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
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
					#tf.keras.layers.Dense(24, input_shape=(24,)),
					tf.keras.layers.Dense(47, input_shape=(47,)),
					tf.keras.layers.Dense(256, activation='relu'),
					tf.keras.layers.Dense(246, activation='relu'),
					tf.keras.layers.Dense(10, activation='softmax')
	])
	model.compile(optimizer = opt,
				  loss=tf.keras.losses.SparseCategoricalCrossentropy(from_logits=True),
				  metrics=['accuracy']
				  )
	return model

def pack(features, label):
  return tf.stack(list(features.values()), axis=-1), label

def smooth(data):
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
		counts = np.count_nonzero(var[:,max(i-19, 0) : i+20] == 2, axis=1)
		idx = np.argmax(counts)
		#idx = np.argmax(np.count_nonzero(var[:,max(i-2, 0) : i+3] == 2, axis=1))
		if counts[idx] >= 3:
			out[6*i+idx] = 2
	return out

def gen_series():
  i = 0
  while True:
    size = np.random.randint(0, 10)
    yield [i, np.random.normal(size=(size,))]
    i += 1


if __name__ == '__main__':
	usage = '%s [-opt1, [-opt2, ...]] infile' % __file__
	parser = argparse.ArgumentParser(description='', formatter_class=RawTextHelpFormatter, usage=usage)
	parser.add_argument('infile', type=is_valid_file, help='input file')
	parser.add_argument('-o', '--outfile', action="store", default=sys.stdout, type=argparse.FileType('w'), help='where to write output [stdout]')
	args = parser.parse_args()
	'''
	if args.labels: print("\t".join(['ID','TYPE','GC'] + translate.amino_acids))
		exit()
	'''

	model = create_model('adam')
	model.load_weights('small.ckpt')
	
	contigs = mt.read_fasta(args.infile)
	for header in contigs:
		dataset = tf.data.Dataset.from_generator(
								mt.get_windows,
								args=[contigs[header]],
								output_types=tf.float32,
								output_shapes = (41,) 
								).batch(10)
		#for feature in train.take(1):
		#	print( feature )
		#exit()
	
		p = model.predict(dataset)
		Y = np.argmax(p,axis=-1)
	
		Y = smooth(Y)
	
		#for row in zip(train.iloc[:,0].to_list(), Y):
		for i,row in enumerate(Y):
			if row == 2:
				if i%2:
					print('     CDS             complement(', ((i-1)//2)+1, '..', ((i-1)//2)+3, ')', sep='')
				else:
					print('     CDS             ', (i//2)+1 , '..', (i//2)+3, sep='')



