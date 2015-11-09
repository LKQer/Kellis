# Generates featurized data
# Common data format:
# 
# Row 1: Tab-separated columns
# ROw 2...n+1: n separate rows, with tab-separated value
#
# Output folder: /broad/compbio/maxwshen/data/1-MAKETRAINTEST/...
# 
# Combine into a single data matrix text file using combine.py
# 
# 
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
CURR_BED_FN = ''

_E = {'IMR90': 'E017', 'GM12878': 'E116', 'K562': 'E123'}

def main():
  SIZE = 1000
  FRAC = .80
  ensure_dir_exists(OUT_PATH + 'temp/')

  inp_fold = '/broad/compbio/maxwshen/data/1-MAKETRAINTEST/fgbg/'
  inp_fn = inp_fold + 'IMR90.1.txt'

  INTXNS(inp_fn)
  print 'Processing', CELLTYPE, CHRO

  # get_motifs()
  # get_chromHMM()
  # get_TFs()
  # imr90_narrowpeaks = ['DNase.macs2', 'H2A.Z', 'H2AK5ac', 'H2AK9ac', 'H2BK5ac', 'H2BK12ac', 'H2BK15ac', 'H2BK20ac', 'H2BK120ac', 'H3K4ac', 'H3K4me1', 'H3K4me2', 'H3K4me3', 'H3K9ac', 'H3K9me1', 'H3K9me3', 'H3K14ac', 'H3K18ac', 'H3K23ac', 'H3K27ac', 'H3K27me3', 'H3K36me3', 'H3K56ac', 'H3K79me1', 'H3K79me2', 'H4K5ac', 'H4K8ac', 'H4K20me1', 'H4K91ac']
  # for np in imr90_narrowpeaks:
    # get_narrowPeak(np)
  get_wgbs()

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
        CELLTYPE = CELLTYPE.upper()
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
  print '  DNA Motifs...'
  base_fold = '/broad/compbio/maxwshen/data/wouter/'
  for mfn in ['motifs_A_L/', 'motifs_M_Z/']:
    curr_fn = base_fold + mfn + 'GENOME_chr' + CHRO + '_binary.txt'
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
        labels.append(l + '_' + str(i))

    ensure_dir_exists(OUT_PATH + mfn)
    motif_out_fn = OUT_PATH + mfn + CELLTYPE + '.' + CHRO + '.txt'
    # print 'Writing to', motif_out_fn
    write_data_to_file(motif_out_fn, labels, rows)

  return

def get_chromHMM():
  # Writes a tab-delimited file where each row is an interaction
  # Each column is all the ChromHMM states, times 2 (1 per locus)
  # Each row is sum of ChromHMM states appearances in that 5kb locus
  #   -> Integer valued
  # First half of row corresponds to first loci in intxn, second to second
  print '  ChromHMM...'
  base_path = '/broad/compbio/maxwshen/data/wouter/ChromHMM_data/'
  if NUM_CHROMHMM_STATES not in [15, 18, 25]:
    print 'ERROR: Invalid NUM_CHROMHMM_STATES', NUM_CHROMHMM_STATES
  if NUM_CHROMHMM_STATES == 25:
    curr_fn = base_path + 'imputed/' + E_CELLTYPE + '_25_imputed12marks_chr' + CHRO + '_statebyline.txt'
  if NUM_CHROMHMM_STATES == 15:
    curr_fn = base_path + 'observed/' + E_CELLTYPE + '_15_coreMarks_chr' + CHRO + '_statebyline.txt'
  if NUM_CHROMHMM_STATES == 18:
    curr_fn = base_path + 'observed_aux/' + E_CELLTYPE + '_18_core_K27ac_chr' + CHRO + '_statebyline.txt'

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

  labels = ['ChromHMM' + str(i + 1) + '_0' for i in range(NUM_CHROMHMM_STATES)] + ['ChromHMM' + str(i + 1) + '_1' for i in range(NUM_CHROMHMM_STATES)]

  # Write to output
  out_fold = 'ChromHMM.imputed.' + str(NUM_CHROMHMM_STATES)
  ensure_dir_exists(OUT_PATH + out_fold)
  ch_out_fn = OUT_PATH + out_fold + '/' + CELLTYPE + '.' + CHRO + '.txt'
  # print 'Writing to', ch_out_fn
  write_data_to_file(ch_out_fn, labels, rows)

  return

def get_TFs():
  # Writes a tab-delimited file where each row is an interaction
  # Each column is all the TFs present in the celltype, times 2 (1 per locus)
  # Each row is sum of TF appearances in that 5kb locus
  #   -> Integer valued
  # First half of row corresponds to first loci in intxn, second to second
  print '  TFs...'
  base_path = '/broad/compbio/maxwshen/data/wouter/TFs/'
  curr_fn = base_path + CELLTYPE + '/' + CELLTYPE + '_chr' + CHRO + '_binary.txt'

  data = dict()
  labels = []
  with open(curr_fn) as f:
    for i, line in enumerate(f):
      if i == 1:
        labels = line.split()
      if i > 1:
        data[(i - 2) * WOUTER_STEP] = line.replace('2', '')

  example_line = line.split()
  new_labels = []
  for i in range(len(example_line)):
    if example_line[i] != '2':
      new_labels.append(labels[i])
  labels = new_labels

  rows = []
  for itx in INTXNS:
    row = []
    for loci in itx:
      curr = [0 for i in range(len(labels))]
      for bp in range(loci, loci + INTXN_LEN, WOUTER_STEP):
        cd = [int(s) for s in data[bp].split()]
        for i in range(len(cd)):
          curr[i] += cd[i]
      row += [str(s) for s in curr]
    rows.append(row)

  new_labels = []
  for i in range(2):
    for l in labels:
      new_labels.append(l + '_' + str(i))
  labels = new_labels

  out_fold = 'TFs'
  ensure_dir_exists(OUT_PATH + out_fold)
  ch_out_fn = OUT_PATH + out_fold + '/' + CELLTYPE + '.' + CHRO + '.txt'
  # print 'Writing to', ch_out_fn
  write_data_to_file(ch_out_fn, labels, rows)

  return

def get_narrowPeak(query):
  # For DNase, query should be 'DNase.macs2'
  # For each query, produces an Nx2 matrix where N = number of intxns
  # Values are floats, sum of signal values over 5kb region where at 
  # least 1bp of the peak overlaps
  print '  ' + query + '...'
  base_path = '/broad/compbio/maxwshen/data/narrowPeak/'
  curr_fn = base_path + E_CELLTYPE + '-' + query + '.narrowPeak'

  data = []
  with open(curr_fn) as f:
    for i, line in enumerate(f):
      if line.split()[0] == 'chr' + CHRO:
        data.append(NarrowPeak(line))

  rows = []
  for itx in INTXNS:
    row = []
    for loci in itx:
      curr = 0
      for peak in data:
        if peak.start <= loci <= peak.end or peak.start <= loci + INTXN_LEN <= peak.end:
          curr += peak.signal
      row.append(str(curr))
    rows.append(row)

  labels = [query + '_0', query + '_1']
  ensure_dir_exists(OUT_PATH + 'NarrowPeak')
  ch_out_fn = OUT_PATH + 'NarrowPeak/' + CELLTYPE + '.' + CHRO + '.' + query + '.txt'
  # print 'Writing to', ch_out_fn
  write_data_to_file(ch_out_fn, labels, rows)
  return

def write_bed():
  # Called from get_wgbs
  # Writes interactions into a bed file
  # BED format:
  # chrA startpos endpos uniquename
  global CURR_BED_FN
  CURR_BED_FN = OUT_PATH + 'temp/temp.bed'
  with open(CURR_BED_FN, 'w') as f:
    for i in range(len(INTXNS)):
      itx = INTXNS[i]
      for j in range(len(itx)):
        loci = itx[j]
        f.write('chr' + CHRO + '\t' + str(loci) + '\t' + str(loci + INTXN_LEN) + '\titxn' + str(i) + '-' + str(j) + '\n')
  return

def get_wgbs():
  # Calls bwtool extract on a bed file made from write_bed()
  # Appends data to the bed file. Format:
  # ChrA startpos endpos uniquename num_datapts data_pt1 data_pt2 ... 
  # This is written to the 'pre' file, which is then averaged into 'avg' file
  print '  wgbs...'
  write_bed()
  base_path = '/broad/compbio/maxwshen/data/wgbs/'
  curr_fn = base_path + E_CELLTYPE + '_WGBS_FractionalMethylation.bigwig'
  ensure_dir_exists(OUT_PATH + 'wgbs')
  ch_out_fn = OUT_PATH + 'wgbs/' + CELLTYPE + '.' + CHRO + '.pre.txt'

  status = commands.getstatusoutput('bwtool extract bed ' + CURR_BED_FN + ' ' + curr_fn + ' ' + ch_out_fn)
  bigwig_avg(ch_out_fn)

  # bwtool extract output files can be very large, so remove them 
  status = commands.getstatusoutput('rm -rf ' + ch_out_fn)
  return

def bigwig_avg(pre_fn):
  # Called from get_wgbs
  # Takes a 'pre' file from bwtool extract and averages the values, setting NA = 0
  # Output value is average fractional methylation over the intxn loci
  out_fn = '.'.join(pre_fn.split('.')[:-2]) + '.txt'
  with open(pre_fn) as f:
    with open(out_fn, 'w') as g:
      g.write('wgbs_0\twgbs_1')
      row = []
      for i, line in enumerate(f):
        if i % 2 == 0:
          g.write('\t'.join(row) + '\n')
          row = []
        data = line.split()[5].replace('NA', '0')
        data = [float(s) for s in data.split(',')]
        row.append(str(np.mean(data))) 
      g.write('\t'.join(row) + '\n')
  return

if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start