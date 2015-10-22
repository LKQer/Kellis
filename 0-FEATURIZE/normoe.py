# Normalizes and finds O/E for Hi-C Data from Rao
# Max Shen

import sys, string, datetime, random, copy, os, commands, fnmatch, re, argparse
import numpy as np
from collections import defaultdict

def main():
  # inp_fn = sys.args[0]

  resolution = 5000
  alg = 'KR'
  quality = 'MAPQG0'
  fold = '/broad/compbio/maxwshen/data/hic//5kb_resolution_intrachromosomal/'
  # chrs = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8']
  # chrs = ['chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16']
  chrs = ['chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX']

  for chro in chrs:
    print chro, datetime.datetime.now()
    fnf = fold + chro + '/' + quality + '/' + chro + '_5kb.' 
    raw_fn = fnf + 'RAWobserved'
    norm_fn = fnf + 'KRnorm'
    exp_fn = fnf + 'KRexpected'
    norm_oe(raw_fn, norm_fn, exp_fn, resolution)

  return

def norm_oe(raw_fn, norm_fn, exp_fn, resolution):
  norm = dict()
  with open(norm_fn) as f:
    for i, line in enumerate(f):
      word = line.strip()
      if word != 'NaN':
        norm[i] = float(word)

  exp = dict()
  with open(exp_fn) as f:
    for i, line in enumerate(f):
      word = float(line.strip())
      exp[i] = word

  out_path = '/broad/compbio/maxwshen/data/processed_hic/k562.5kb.mapqg0/'
  out_fn = out_path + raw_fn.split('/')[-1] + 'norm.txt'
  print out_fn
  with open(raw_fn) as f:
    with open(out_fn, 'w') as g:
      for i, line in enumerate(f):
        if i % 1000000 == 0:
          print i, datetime.datetime.now()
        words = line.split()
        start = int(words[0])
        end = int(words[1])
        val = float(words[2])

        ind1 = start / resolution
        ind2 = end / resolution
        if ind1 in norm and ind2 in norm:
          new_val = val / (norm[ind1] * norm[ind2])
          ind_exp = (end - start) / resolution
          new_val /= exp[ind_exp]
          g.write(str(start) + '\t' + str(end) + '\t' + str(new_val) + '\n')

  print 'Done'
  return

if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start