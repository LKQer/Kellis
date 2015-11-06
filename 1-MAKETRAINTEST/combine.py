# Combines features produced by maketraintest.py 
# into a single data matrix text file

import sys, string, datetime, random, copy, os, commands, fnmatch, re, argparse
import numpy as np
from collections import defaultdict
sys.path.append('/broad/compbio/maxwshen/Kellis/util')
from lib import *

def main():

  name = sys.argv[1]
  cell_types = ['IMR90']
  chros = ['1']

  common_fold = '/broad/compbio/maxwshen/data/1-MAKETRAINTEST/'
  mtt_fold = common_fold + 'traintest/'
  out_path = common_fold + 'combined/'
  ensure_dir_exists(out_path)

  # Ensure out file exists and is empty for appending
  out_fn = out_path + name + '.txt'
  with open(out_fn, 'w') as f:
    pass

  for cell_type in cell_types:
    for chro in chros:
      print cell_type, chro

      data = []   # List of columns
      for subdir in [mtt_fold + s + '/' for s in os.listdir(mtt_fold) if os.path.isdir(mtt_fold + s)]:
        for seek_fn in os.listdir(subdir):
          if fnmatch.fnmatch(seek_fn, cell_type + '.' + chro + '*.txt'):
            data = data_add(data, subdir + seek_fn)

      print 'Found', len(data), 'features/columns'
      data_rows = []
      for i in range(len(data[0])):
        data_rows.append([data[j][i] for j in range(len(data))])
      with open(out_fn, 'a') as f:
        for row in data_rows:
          f.write('\t'.join(row) + '\n')
  return


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