# Generates training and testing data
# Includes lots of magic numbers reflecting internal organization of featurized data
# Max Shen

import sys, string, datetime, random, copy, os, commands, fnmatch, re, argparse
import numpy as np
from collections import defaultdict
sys.path.append('/broad/compbio/maxwshen/Kellis/util')
from lib import *

OUT_PATH = '/broad/compbio/maxwshen/data/1-MAKETRAINTEST/traintest/'

INTXN_LEN = 5000
WOUTER_STEP = 200
NUM_CHROMHMM_STATES = 18            # VALID:  15, 18 or 25

_E = {'imr90': 'E017', 'gm12878': 'E116', 'k562': 'E123'}

def main():
  SIZE = 1000
  FRAC = .80

  inp_fold = '/broad/compbio/maxwshen/data/1-MAKETRAINTEST/fgbg/'
  inp_fn = inp_fold + 'sample.imr90.1.txt'

  INTXNS(inp_fn)
  print 'Processing', CELLTYPE, CHRO

  # get_motifs()
  get_chromHMM()

  return

def INTXNS(inp_fn):
  # Grabs interactions from the output of fgbg.py.
  # Format = Locus1 Locus2 IntxnFraction
  global INTXNS
  global CELLTYPE
  global CHRO
  global E_CELLTYPE

  INTXNS = []
  with open(inp_fn) as f:
    for i, line in enumerate(f):
      if i == 0:
        [_, CELLTYPE, CHRO, limit] = line.split()
        E_CELLTYPE = _E[CELLTYPE]
      if i > 0:
        INTXNS.append([int(s) for s in line.split()[:2]])
  return

def write_data_to_file(out_fn, labels, rows):
  # Expects labels to be a list of strings
  # Expects rows to be a list of strings
  with open(out_fn, 'w') as f:
    f.write('\t'.join(labels) + '\n')
    for row in rows:
      # print len(row)
      f.write('\t'.join(row) + '\n')
  return

def get_motifs():
  # Writes a tab-delimited file where each row is an interaction
  # Each column is all the motifs, times 2 (1 per locus)
  # Each row is total motif appearances in that 5kb locus
  #   -> Integer valued
  # First half of row corresponds to first loci in intxn, second to second
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

    ensure_dir_exists(OUT_PATH + mfn)
    motif_out_fn = OUT_PATH + mfn + CELLTYPE + '.' + CHRO + '.txt'
    print 'Writing to', motif_out_fn
    write_data_to_file(motif_out_fn)

  return

def get_chromHMM():
  # Writes a tab-delimited file where each row is an interaction
  # Each column is all the ChromHMM states, times 2 (1 per locus)
  # Each row is sum of ChromHMM states appearances in that 5kb locus
  #   -> Integer valued
  # First half of row corresponds to first loci in intxn, second to second
  ch_path = '/broad/compbio/maxwshen/data/wouter/ChromHMM_data/'
  if NUM_CHROMHMM_STATES not in [15, 18, 25]:
    print 'ERROR: Invalid NUM_CHROMHMM_STATES', NUM_CHROMHMM_STATES
  if NUM_CHROMHMM_STATES == 25:
    curr_fn = ch_path + 'imputed/' + E_CELLTYPE + '_25_imputed12marks_chr' + CHRO + '_statebyline.txt'
  if NUM_CHROMHMM_STATES == 15:
    curr_fn = ch_path + 'observed/' + E_CELLTYPE + '_15_coreMarks_chr' + CHRO + '_statebyline.txt'
  if NUM_CHROMHMM_STATES == 18:
    curr_fn = ch_path + 'observed_aux/' + E_CELLTYPE + '_18_core_K27ac_chr' + CHRO + '_statebyline.txt'

  data = dict()
  with open(curr_fn) as f:
    for i, line in enumerate(f):
      if i > 2:
        data[(i - 3) * WOUTER_STEP] = line.strip()

  rows = []
  for itx in INTXNS:
    row = [] 
    for loci in itx:
      curr = [0 for i in range(NUM_CHROMHMM_STATES)]
      for bp in range(loci, loci + INTXN_LEN, WOUTER_STEP):
        curr[int(data[bp]) - 1] += 1
      row += [str(s) for s in curr]
    rows.append(row)

  labels = ['ChromHMM ' + str(i + 1) + '_0' for i in range(NUM_CHROMHMM_STATES)] + ['ChromHMM ' + str(i + 1) + '_1' for i in range(NUM_CHROMHMM_STATES)]

  # Write to output
  out_fold = 'ChromHMM.imputed.' + str(NUM_CHROMHMM_STATES)
  ensure_dir_exists(OUT_PATH + out_fold)
  ch_out_fn = OUT_PATH + out_fold + '/' + CELLTYPE + '.' + CHRO + '.txt'
  print 'Writing to', ch_out_fn
  write_data_to_file(ch_out_fn)

  return

if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start