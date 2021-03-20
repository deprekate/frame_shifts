import os
import sys
import argparse
from argparse import RawTextHelpFormatter
from statistics import mode

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

	letters = ['#', '*', '+', 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
	tfiles = tf.data.experimental.make_csv_dataset(
		file_pattern = "test2/*.tsv",
		field_delim='\t',
		header=False,
		column_names= ['ID', 'TYPE','GC'] + [letter for pair in zip([l+'a' for l in letters], [l+'b' for l in letters]) for letter in pair],
		batch_size=8, num_epochs=1,
		num_parallel_reads=20,
		shuffle_buffer_size=10000,
		select_columns=['TYPE','GC'] + [letter for pair in zip([l+'a' for l in letters], [l+'b' for l in letters]) for letter in pair],
		label_name='TYPE'
		)

	pdata = tfiles.map(pack)
	#for feature in tfiles.take(1):
	#	print( feature )

	cp_callback = tf.keras.callbacks.ModelCheckpoint(filepath='cp.ckpt', save_weights_only=True, verbose=1)
	model = create_model('adam')
	model.fit(pdata, epochs=9, callbacks=[cp_callback])

	#model = create_model('adam')
	#model.load_weights('cp.ckpt')
	
	
	train = pd.read_csv(sys.argv[1], header=None, sep='\t')
	X = tf.stack(train.iloc[:,2:])
	#Y = train.iloc[:,1].replace({'None':0, 'False':1, 'True':2})
	p = model.predict(X)
	f = np.argmax(p,axis=-1)
	#smooth(f)
	ff = f #xsmooth(f)
	
	for row in zip(train.iloc[:,0].to_list(), ff):
		if row[1] == 2:
			if row[0] > 0:
				print('     CDS             ', row[0] , '..', row[0]+2, sep='')
			else:
				print('     CDS             complement(', abs(row[0]), '..', abs(row[0])+2, ')', sep='')



