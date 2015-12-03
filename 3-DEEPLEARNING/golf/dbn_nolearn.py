# DBN with nolearn library
# Written for personal spot AWS EC2

import sys, string, datetime, random, copy, os, commands
import fnmatch, re, argparse, csv, pickle
import numpy as np
from collections import defaultdict

sys.path.append('/broad/compbio/maxwshen/Kellis/util')
from lib import *

from sklearn.metrics import classification_report
from nolearn.dbn import DBN


def main():
  data_id = 'B'
  data_path = '/broad/compbio/maxwshen/data/1-MAKETRAINTEST/complete/golf/'
  
  print 'train...', datetime.datetime.now()
  train_set = readin(data_id, 'train', data_path)
  print 'valid...', datetime.datetime.now()
  valid_set = readin(data_id, 'valid', data_path)
  print 'test...', datetime.datetime.now()
  test_set = readin(data_id, 'test', data_path)

  # Input to 300 node RBM to 2 node output
  dbn = DBN( \
    [xtrain.shape[1], 300, 2], \
    learn_rates = 5, \
    learn_rate_decays = 0.9, \
    epochs = 31, \
    verbose = 1)
  dbn.fit(dat_train, y_train)

  preds = dbn.predict(dat_test)
  print classification_report(y_test, preds)

  out_fn = 'dbn.pickle'
  with open(out_fn, 'w') as f:
    pickle.dump(dbn, out_fn)

  return


def readin(data_id, typ, data_path):
  dat = read_delimited_txt(data_path + data_id + '.' + typ + '.dat.txt', '\t')
  dat = np.asarray(get_names_trim(dat))
  y = read_delimited_txt(data_path + data_id + '.' + typ + '.y.txt', '\t')
  y = np.asarray(convert_TF_to_num(get_names_trim(y)))
  return (dat, y)


def convert_TF_to_num(nparr):
    for i in range(len(nparr)):
        for j in range(len(nparr[0])):
            if nparr[i][j] == 'FALSE':
                nparr[i][j] = 0
            if nparr[i][j] == 'TRUE':
                nparr[i][j] = 1
    return nparr

def get_names_trim(dat):
  names = dat[0]
  dat = dat[1:]
  for i in range(len(dat)):
    dat[i] = dat[i][1:]
  return dat


if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start