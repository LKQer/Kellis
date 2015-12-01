# Find foreground and background interactions
# Output is used by maketraintest.py
# Max Shen

import sys, string, datetime, random, copy, os, commands, fnmatch, re, argparse
import pdb, pickle
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
  global MAX_DIST
  global BG_MIN
  global BG_MAX
  global FG_MIN
  global CHRO
  global CELLTYPE
  _NUM = float('inf')     # Find this many of each set (foreground, background)
  MAX_DIST = 1000000      # maximum intxn distance to pull
  LIMIT = float('inf')    # No early stopping
  FG_MIN = 15
  BG_MIN = -1
  BG_MAX = 0.6

  # Using NATO phonetic alphabet
  name = sys.argv[1]
  celltype = sys.argv[2]
  cnum = sys.argv[3]
  OUT_PATH = OUT_PATH + name + '/'

  celltypes = [celltype]
  # celltypes = ['IMR90', 'GM12878', 'K562']
  # No chr9 (missing data), chr X
  # chrs = ['22', '21', '20', '19', '18', '17', '16', '15', '14', \
          # '13', '12', '11', '10', '8', '7', '6', '5', \
          # '4', '3', '2', '1']
  # chrs = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11']
  chros = [['1', '2', '3', '4', '5', '6', '7'], ['8', '10', '11', '12', '13', '14', '15'], ['16', '17', '18', '19', '20', '21', '22']]
  # celltypes = ['IMR90', 'K562']
  chrs = chros[int(cnum)]

  for ct in celltypes:
    for chro in chrs:
      CHRO = chro
      CELLTYPE = ct
      print ct, chro
      datapath = DataPath(ct, chro)
      find_fgbg(datapath, name)

  return


def find_fgbg(datapath, name):
  print '  Reading Hi-C intxn file...'
  intxns = dict()
  with open(datapath.path) as f:
    for i, line in enumerate(f):
      pass
      if i % 6000000 == 0:
        print '    ', i, datetime.datetime.now()
      if i > LIMIT:
        break
      words = line.split()
      intxns[' '.join(words[:2])] = float(words[2])

  print '  Sorting keys..'
  ordered_keys = sorted(intxns, key = intxns.get)
  # high = ordered_keys[-_NUM:]
  # low = ordered_keys[:_NUM]

  # Used to find independent intxns within each chromosomes
  keys = filter_match_distances(intxns, ordered_keys)

  ensure_dir_exists(OUT_PATH)
  write_file(keys, intxns)
  return
  

def write_file(keys, intxns):
  # type is either 'train' or 'test'
  # Also mirrors the interactions
  out_fn = OUT_PATH + CELLTYPE + '.' + CHRO + '.txt'
  with open(out_fn, 'w') as f:
    f.write('>  ' + ' '.join([CELLTYPE, CHRO, str(LIMIT)]) + '\n')
    for key in keys:
      f.write(key + '\t' + str(intxns[key]) + '\n')
      f.write(key.split()[1] + ' ' + key.split()[0] + '\t' + str(intxns[key]) + '\n')
  print 'Written to', out_fn
  return


def convert_to_loci(keys):
  # keys is a list of space-delimited strings with 2 components
  # loci is a set of strings
  loci = set()
  for itx in keys:
    loci.add(itx.split()[0])
    loci.add(itx.split()[1])
  return loci


def loci_match(itx, loci):
  # itx is a space-delimited string with 2 components
  # loci is a list of strings
  return itx.split()[0] in loci or itx.split()[1] in loci


def count_loci_enrich(typ, keys):
  # Want to ensure some loci are not overrepresented compared to others
  # Takes in a list of space-delimited strings (keys)
  # Builds a dictionary counting occurences of loci
  counter = dict()
  for k in keys:
    loci = k.split()
    for locus in loci:
      if locus not in counter:
        counter[locus] = 0
      counter[locus] += 1

  with open('ex' + typ + '.txt', 'w') as f:
    for key in counter:
      f.write(key + ' ' + str(counter[key]) + '\n')
  return



def get_dist(key):
  return abs(int(key.split()[1]) - int(key.split()[0]))


def filter_match_distances(intxns, ordered_keys):
  # Finds the top _NUM intxns within MAX_DIST and matches (if possible) 1-to-1
  # between foreground and background
  # Also ensures foreground and background sets are within the set range
  # Store keys in high, and store their dists
  print '  Building fg/bg by matching distances...'
  high = []
  dists = []
  for i in range(len(ordered_keys) - 1, 0, -1):
    ok = ordered_keys[i]
    dist = get_dist(ok)
    # print dist, ok, intxns[ok]
    # pdb.set_trace()
    if dist <= MAX_DIST:
      high.append(ok)
      dists.append(dist)
    if len(high) == _NUM:
      break
    if intxns[ok] < FG_MIN:
      print 'ERROR: Insufficient number of FOREGROUND', CELLTYPE, CHRO, '_NUM:', _NUM, 'found:', len(high), 'processed:', len(ordered_keys) - i
      break
      # sys.exit(0)

  low = []
  for i in range(len(ordered_keys)):
    ok = ordered_keys[i]
    dist = get_dist(ok)
    if dist in dists:
      # Find an intxn that is unique to all fg/bg intxns so far
      # matches a distance in foreground 1-to-1
      dists.remove(dist)
      low.append(ok)
    if len(dists) == 0:
      # Ensures we meet _NUM intxns      
      break
    if intxns[ok] > BG_MAX:
      print 'Insufficient number of BACKGROUND', CELLTYPE, CHRO, '_NUM:', _NUM, 'found:', len(low), 'processed:', i
      break
      # sys.exit(0)

  # If there are unmatched background/foreground...
  hc = copy.copy(high)
  for h in hc:
    if get_dist(h) in dists:
      dists.remove(get_dist(h))
      high.remove(h)
    if len(dists) == 0:
      break

  if len(high) != len(low):
    print 'ERROR: Mismatching sizes of high and low'
  print '  Found:', len(high), len(low)

  return high + low



if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start