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
  global OUT_PATH
  name = sys.argv[1]
  inp_fold = '/broad/compbio/maxwshen/data/1-MAKETRAINTEST/fgbg/' + name + '/'
  OUT_PATH = OUT_PATH + name + '/'

  ensure_dir_exists(OUT_PATH + 'temp/')

  # EDIT 
  for inp_fn in os.listdir(inp_fold):
    INTXNS(inp_fold + inp_fn)
    print 'Processing', CELLTYPE, CHRO

    get_motifs()
    get_chromHMM()
    get_TFs()
    get_all_narrowpeak()
    get_wgbs()
    get_gc()
    get_entropy()
    get_basic()

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

def get_basic():
  # Stores basic information into the matrix
  # Genomic distance of interaction, cell type, chromosome
  print '  Basic...'
  rows = []
  for itx in INTXNS:
    rows.append([abs(itx[0] - itx[1])])

  for row in rows:
    row.append(CELLTYPE)
    row.append(CHRO)
  rows = ['\t'.join(row) for row in rows]

  labels = ['basic_genomic_dist', 'basic_celltype', 'basic_chromosome']

  ensure_dir_exists(OUT_PATH + 'basic')
  out_fn = OUT_PATH + 'basic/' + CELLTYPE + '.' + CHRO + '.txt'
  write_data_to_file(out_fn, labels, rows)
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

def get_all_narrowpeak():
  # Called from main
  for fn in os.listdir('/broad/compbio/maxwshen/data/narrowPeak/'):
    if fnmatch.fnmatch(fn, E_CELLTYPE + '*'):
      get_narrowPeak(fn)
  return

def get_narrowPeak(query):
  # For DNase, query should be 'DNase.macs2'
  # For each query, produces an Nx2 matrix where N = number of intxns
  # Values are floats, sum of signal values over 5kb region where at 
  # least 1bp of the peak overlaps
  print '  ' + query + '...'
  base_path = '/broad/compbio/maxwshen/data/narrowPeak/'
  curr_fn = base_path + query

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


def get_gc():
    print '  GC Content...'
    base_path = '/broad/compbio/maxwshen/data/hg19/'        
    #base_path = './'
    filename = base_path + 'chr' + CHRO + '.fa'
    seq = readSeq(filename)

    rows = []    
    for itx in INTXNS:
        row = []
        for loci in itx:
            gcnum = 0
        
            start  = loci
            stop = loci + INTXN_LEN
        
            for i in xrange(start,stop):
                letter = seq[i]
            
                if ((letter.upper() == 'G') | (letter.upper() == 'C')):
                    gcnum = gcnum + 1                
        
            gc_content = float(gcnum)/(stop-start)        
            row.append(str(gc_content))
        rows.append(row)
    
    # write to file    
    out_fn = OUT_PATH + 'gc/' + CELLTYPE + '.' + CHRO + '.txt'
    query = 'gc'
    labels = [query + '_0', query + '_1']
    write_data_to_file(out_fn, labels, rows)
    
    return

def get_entropy():
    print '  Entropy...'
    base_path = '/broad/compbio/maxwshen/data/hg19/'        
    #base_path = './'
    filename = base_path + 'chr' + CHRO + '.fa'
    seq = readSeq(filename)
    
    alphabet = ['A','C','G','T']
    # to avoid log 0
    eps = sys.float_info.epsilon
    
    rows = []    
    
    for itx in INTXNS:
        row = []
        for loci in itx:
            gcnum = 0
        
            start  = loci
            stop = loci + INTXN_LEN
        
            st = seq[start:stop].upper()
            stList = list(st)
        
            # make the alphabet - not sure if I should just have constant alphabet
            #alphabet = list(Set(stList))
            freqList = []
            for symbol in alphabet:
                ctr = 0
                for sym in stList:
                    if sym == symbol:
                        ctr += 1
                freqList.append(float(ctr) / len(stList))

            ent = 0.0
        
            for freq in freqList:
                freq = max(freq,eps)
                ent = ent + freq * math.log(freq, 2)
            ent = -ent
        
            # normalize
            ent_max = math.log(len(freqList), 2)
            ent = ent/ent_max
        
            row.append(str(ent))  
        rows.append(row)
    
    # write to file    
    out_fn = OUT_PATH + 'repeat/' + CELLTYPE + '.' + CHRO + '.txt'
    query = 'entropy'
    labels = [query + '_0', query + '_1']
    write_data_to_file(out_fn, labels, rows)
    
    return    
    

if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start