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
WOUTER_STEP = 200

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
  # Writes a tab-delimited file where each row is an interaction
  # Each column is all the motifs
  # Each entry is the number of times a motif appears in the intxn
  #   -> Integer valued
  motif_path = '/broad/compbio/maxwshen/data/wouter/'
  for mfn in ['motifs_A_L/', 'motifs_M_Z/']:
    curr_fn = motif_path + mfn + 'GENOME_chr' + CHRO + '_binary.txt'
    data = dict()
    col_labels = ''
    with open(curr_fn) as f:
      for i, line in enumerate(f):
        if i == 1:
          col_labels = line.split()
          num_cols = len(col_labels)
        if i > 1:
          data[(i - 2) * WOUTER_STEP] = line.strip()

    rows = []
    for itx in INTXNS:
      row = []
      for loci in itx:  
        curr_loci = [0 for s in range(num_cols)]
        for bp in range(loci, loci + INTXN_LEN, WOUTER_STEP):
          cd = [int(s) for s in data[bp].split()]
          for i in range(len(cd)):
            curr_loci[i] += cd[i]
        row += [str(s) for s in curr_loci]
      rows.append(row)

    labels = []
    for i in range(2):
      for l in col_labels:
        print l
        labels.append(l + '_' + str(i))

    ensure_dir_exists(out_path + mfn)
    motif_out_fn = out_path + mfn + CELLTYPE + '.' + CHRO + '.txt'
    print 'Writing to', motif_out_fn
    with open(motif_out_fn, 'w') as f:
      f.write('\t'.join(labels) + '\n')
      for row in rows:
        # print len(row)
        f.write('\t'.join(row) + '\n')

  return

if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start