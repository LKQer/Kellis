# Generates training and testing data
# Includes lots of magic numbers reflecting internal organization of featurized data
# Max Shen

import sys, string, datetime, random, copy, os, commands, fnmatch, re, argparse
import numpy as np
from collections import defaultdict
sys.path.append('/broad/compbio/maxwshen/Kellis/util')
from lib import *

out_path = '/broad/compbio/maxwshen/data/1-MAKETRAINTEST/traintest/'

INTXN_LEN = 5000

def main():
  SIZE = 1000
  FRAC = .80

  inp_fold = '/broad/compbio/maxwshen/data/1-MAKETRAINTEST/fgbg/'
  inp_fn = inp_fold + 'sample.imr90.1.txt'

  INTXNS(inp_fn)
  print 'Processing', CELLTYPE, CHRO

  get_motifs()

  return

def INTXNS(inp_fn):
  global INTXNS
  global CELLTYPE
  global CHRO

  INTXNS = []
  with open(inp_fn) as f:
    for i, line in enumerate(f):
      if i == 0:
        [_, CELLTYPE, CHRO, limit] = line.split()
      if i > 0:
        INTXNS.append([int(s) for s in line.split()[:2]])
  return

def get_motifs():
  motif_path = '/broad/compbio/maxwshen/data/wouter/'
  step = 200
  for mfn in ['motifs_A_L/', 'motifs_M_Z/']:
    curr_fn = motif_path + mfn + 'GENOME_chr' + CHRO + '_binary.txt'
    data = dict()
    with open(curr_fn) as f:
      for i, line in enumerate(f):
        data[i * step] = line.strip()

    rows = []
    for itx in INTXNS:
      row = []
      for loci in itx:  
        bps = range(loci, loci + INTXN_LEN - step, step)
        for bp in bps:
          row.append(data[bp])
      rows.append('\t'.join(row))

    ensure_dir_exists(out_path + mfn)
    motif_out_fn = out_path + mfn + CELLTYPE + '.' + CHRO + '.txt'
    print 'Writing to', motif_out_fn
    with open(motif_out_fn, 'w') as f:
      for row in rows:
        # print len(row)
        f.write(row + '\n')

  return

if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start