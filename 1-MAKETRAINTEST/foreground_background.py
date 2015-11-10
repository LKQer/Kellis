# Find foreground and background interactions
# Output is used by maketraintest.py
# Max Shen

import sys, string, datetime, random, copy, os, commands, fnmatch, re, argparse
import numpy as np
import operator
from collections import defaultdict
sys.path.append('/broad/compbio/maxwshen/Kellis/util')
from lib import *

OUT_PATH = '/broad/compbio/maxwshen/data/1-MAKETRAINTEST/fgbg/'

def main():
  global _NUM
  global LIMIT
  global OUT_PATH
  _NUM = 1000
  LIMIT = float('inf')
  
  # Using NATO phonetic alphabet
  name = sys.argv[1]
  OUT_PATH = OUT_PATH + name + '/'

  # celltypes = ['IMR90', 'GM12878', 'K562']
  celltypes = ['K562']
  chrs = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', \
          '12', '13', '14', '15', '16', '17', '18', '19', '20', \
          '21', '22', 'X']

  for ct in celltypes:
    for chro in chrs:
      print ct, chro
      datapath = DataPath(ct, chro)
      find_fgbg(datapath, name)

  return

def find_fgbg(datapath, name):
  intxns = dict()
  with open(datapath.path) as f:
    for i, line in enumerate(f):
      pass
      if i % 50000 == 0:
        print i, datetime.datetime.now()
      if i > LIMIT:
        break

      words = line.split()
      intxns[' '.join(words[:2])] = float(words[2])

  ordered_keys = sorted(intxns, key = intxns.get)
  high = ordered_keys[-_NUM:]
  low = ordered_keys[:_NUM]

  ensure_dir_exists(OUT_PATH)
  out_fn = OUT_PATH + datapath.celltype + '.' + datapath.chro + '.txt'
  with open(out_fn, 'w') as f:
    f.write('> ' + ' '.join([datapath.celltype, datapath.chro, str(LIMIT)]) + '\n')
    for key in high:
      f.write(key + '\t' + str(intxns[key]) + '\n')
    for key in low:
      f.write(key + '\t' + str(intxns[key]) + '\n')
  print 'Written to', out_fn
  return

if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start