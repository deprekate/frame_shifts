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
					tf.keras.layers.Dense(41, input_shape=(41,)),
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


if __name__ == '__main__':
	usage = '%s [-opt1, [-opt2, ...]] directory' % __file__
	parser = argparse.ArgumentParser(description='', formatter_class=RawTextHelpFormatter, usage=usage)
	parser.add_argument('directory', type=is_valid_file, help='input directory')
	parser.add_argument('-o', '--outfile', action="store", default=sys.stdout, type=argparse.FileType('w'), help='where to write the output [stdout]')
	args = parser.parse_args()

	letters = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
	colnames = ['TYPE','GC'] + [letter for pair in zip([l+'1' for l in letters], [l+'2' for l in letters]) for letter in pair]
	tfiles = tf.data.experimental.make_csv_dataset(
		file_pattern        = args.directory + "/*.tsv",
		field_delim         = '\t',
		header              = False,
		column_names        = colnames,
		column_defaults     = [tf.int32] + [tf.float32] * 41,
		batch_size          = 100,
		num_epochs          = 1,
		num_parallel_reads  = 10,
		shuffle_buffer_size = 10000,
		label_name          = colnames[0],
		select_columns      = colnames
		)
	pdata = tfiles.map(pack)
	#for feature in tfiles.take(1):
	#	print( feature )

	cp_callback = tf.keras.callbacks.ModelCheckpoint(filepath=args.directory + '.ckpt', save_weights_only=True, verbose=1)
	model = create_model('adam')
	model.fit(pdata, epochs=5, callbacks=[cp_callback])



