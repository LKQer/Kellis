# Combines features produced by maketraintest.py 
# into a single data matrix text file

import sys, string, datetime, random, copy, os, commands, fnmatch, re, argparse, pdb
import numpy as np
from collections import defaultdict
sys.path.append('/broad/compbio/maxwshen/Kellis/util')
from lib import *

def main():
  global MTT_FOLD
  global INP_PATH
  global OUT_FN
  global OUT_Y_FN
  global ADD_LABELS
  global TYP

  name = sys.argv[1]
  TYP = sys.argv[2]
  # cell_types = ['IMR90', 'GM12878', 'K562']
  cell_types = ['IMR90', 'K562']
  # chros = ['1', '2', '3', '4', '5', '6', '7', '8', \
          # '10', '11', '12', '13', '14', '15', '16', '17', \
          # '18', '19', '20', '21', '22']
  chros = ['1']

  common_fold = '/broad/compbio/maxwshen/data/1-MAKETRAINTEST/'
  MTT_FOLD = common_fold + 'traintest/' + name + '/' + TYP + '/'
  out_path = common_fold + 'combined/' + name + '/' + TYP + '/'
  INP_PATH = common_fold + 'fgbg/' + name + '/' + TYP + '/'
  ensure_dir_exists(out_path)


  for cell_type in cell_types:
    # Ensure out file exists and is empty for appending
    OUT_FN = out_path + 'data.' + cell_type + '.txt'
    OUT_Y_FN = out_path + 'y.' + cell_type + '.txt'
    prepare_OUT_FN(OUT_FN)
    prepare_OUT_FN(OUT_Y_FN)
    ADD_LABELS = True

    for chro in chros:
      print cell_type, chro
      num_f = add_features(cell_type, chro)
      num_y = add_y(cell_type, chro)
      ADD_LABELS = False
      if num_f != num_y:
        print 'Error: Number of rows in features and output do not agree.'
        print 'Num of rows in features:', num_f, '\t Num y:', num_y

  return

def add_y(cell_type, chro):
  curr_fn = INP_PATH + TYP + '.' + cell_type + '.' + chro + '.txt'
  ys = []
  with open(curr_fn) as f:
    for i, line in enumerate(f):
      if i == 0:
        ys.append('Labels')
        if line[0] != '>':
          print 'Error: Bad fgbg file format'
      if i > 0:
        ys.append(line.split()[2])

  if not ADD_LABELS:
    ys = ys[1:]

  with open(OUT_Y_FN, 'a') as f:
    for s in ys:
      f.write(s + '\n')

  return len(ys)

def prepare_OUT_FN(OUT_FN):
  with open(OUT_FN, 'w') as f:
    pass
  return

def add_features(cell_type, chro):
  # Handle features
  data = []   # List of columns
  for subdir in [MTT_FOLD + s + '/' for s in os.listdir(MTT_FOLD) if os.path.isdir(MTT_FOLD + s)]:
    for seek_fn in os.listdir(subdir):
      if fnmatch.fnmatch(seek_fn, cell_type + '.' + chro + '.*txt'):
        data = data_add(data, subdir + seek_fn)

  print 'Found', len(data), 'features/columns'
  data_rows = []
  for i in range(len(data[0])):
    # pdb.set_trace()
    data_rows.append([data[j][i] for j in range(len(data))])

  if not ADD_LABELS:
    data_rows = data_rows[1:]

  with open(OUT_FN, 'a') as f:
    for row in data_rows:
      f.write('\t'.join(row) + '\n')
  return len(data_rows)

def data_add(data, seek_fn):
  # print seek_fn
  rows = []
  with open(seek_fn) as f:
    for i, line in enumerate(f):
      rows.append(line.split())
  for i in range(len(rows[0])):
    data.append([rows[j][i] for j in range(len(rows))])

  return data


if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start