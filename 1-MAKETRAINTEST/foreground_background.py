# Find foreground and background interactions
# Output is used by maketraintest.py
# Max Shen

import sys, string, datetime, random, copy, os, commands, fnmatch, re, argparse
import numpy as np
import operator
from collections import defaultdict
sys.path.append('/broad/compbio/maxwshen/Kellis/util')
from lib import *

out_path = '/broad/compbio/maxwshen/data/1-MAKETRAINTEST/fgbg/'

def main():
  global _NUM
  global limit
  _NUM = 500
  limit = 10000
  
  celltypes = ['IMR90']
  chrs = [1]
  for ct in celltypes:
    for chro in chrs:
      datapath = DataPath(ct, chro)
      find_fgbg(datapath)

  return

def find_fgbg(datapath):
  high = dict()
  low = dict()
  with open(datapath.path) as f:
    for i, line in enumerate(f):
      if i % 50000 == 0:
        print i, datetime.datetime.now()
      if i > limit:
        break

      intxn = ' '.join(line.split()[:2])
      num = float(line.split()[2])
      if len(high) < _NUM:
        high[intxn] = num
      else:
        if min(high.values()) < num:
          del high[min(high.iteritems(), key=operator.itemgetter(1))[0]]
          high[intxn] = num
      if len(low) < _NUM:
        low[intxn] = num
      else:
        if max(low.values()) > num:
          del low[max(low.iteritems(), key=operator.itemgetter(1))[0]]
          low[intxn] = num


  ensure_dir_exists(out_path)
  out_fn = out_path + datapath.celltype + '.' + datapath.chro + '.txt'
  with open(out_fn, 'w') as f:
    f.write('> ' + ' '.join([datapath.celltype, datapath.chro, str(limit)]) + '\n')
    for key in high:
      f.write(key + '\t' + str(high[key]) + '\n')
    for key in low:
      f.write(key + '\t' + str(low[key]) + '\n')
  print 'Written to', out_fn
  return

if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start