# DBN with nolearn library
# Written for personal spot AWS EC2

import sys, string, datetime, random, copy, os, commands
import fnmatch, re, argparse, csv, pickle
import numpy as np
from collections import defaultdict

sys.path.append('/home/ec2-user/Kellis/util')
from lib import *

from sklearn.cross_validation import train_test_split
from sklearn.metrics import classification_report
from nolearn.dbn import DBN

def main():
  data_fn = '/home/ec2-user/Kellis/data/bravo.formatted/dat.all.txt'
  blacklist_fn = '/home/ec2-user/Kellis/data/bravo.formatted/dat.blacklist.txt'
  y_fn = '/home/ec2-user/Kellis/data/bravo.formatted/dat.y.txt'

  data = read_delimited_txt(data_fn, '\t')
  blacklist = read_delimited_txt(blacklist_fn, '\t')
  y = read_delimited_txt(y_fn, '\t')
  
  # Get names and remove the first element of each row which is the row number
  names = data[0]
  data = data[1:]
  for i in range(len(data)):
    data[i] = data[i][1:]

  y = y[1:]
  for i in range(len(y)):
    y[i] = y[i][-1]
  y = convert_y_binary(y)

  # Normalizes column-wise so all values are between 0 and 1
  data = normalize_0_1(data)

  # Split into training, testing
  xtrain, xtest, ytrain, ytest = train_test_split(data, y, test_size = 0.2, random_state = 1)

  # Input to 300 node RBM to 2 node output
  dbn = DBN( \
    [xtrain.shape[1], 300, 2], \
    learn_rates = 5, \
    learn_rate_decays = 0.9, \
    epochs = 501, \
    verbose = 1)
  dbn.fit(xtrain, ytrain)

  preds = dbn.predict(xtest)
  print classification_report(ytest, preds)

  out_fn = 'dbn.pickle'
  with open(out_fn, 'w') as f:
    pickle.dump(dbn, out_fn)

  return


def convert_y_binary(y):
  biny = []
  for s in y:
    biny.append(int(float(s) > 1))
  return biny

if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start