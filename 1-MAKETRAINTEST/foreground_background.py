# Find foreground and background interactions
# Output is used by maketraintest.py
# Max Shen

import sys, string, datetime, random, copy, os, commands, fnmatch, re, argparse
import pdb
# import numpy as np
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
  global TEST_FRAC
  global BG_MIN
  global BG_MAX
  global FG_MIN
  global CHRO
  global CELLTYPE
  _NUM = 500             # Find this many of each set (foreground, background)
  MAX_DIST = 1500000      # maximum intxn distance to pull
  LIMIT = float('inf')    # No early stopping
  FG_MIN = 10
  BG_MIN, BGMAX = (0.5, 1.5)
  TEST_FRAC = 0.1         # Percent of test split

  # Using NATO phonetic alphabet
  name = sys.argv[1]
  OUT_PATH = OUT_PATH + name + '/'

  celltypes = ['IMR90', 'GM12878', 'K562']
  # celltypes = ['K562']
  chrs = ['1', '2', '3', '4', '5', '6', '7', '8', '10', '11', \
          '12', '13', '14', '15', '16', '17', '18', '19', '20', \
          '21', '22', 'X']
  # chrs = ['22']

  for ct in celltypes:
    for chro in chrs:
      CHRO = chro
      CELLTYPE = ct
      print ct, chro
      datapath = DataPath(ct, chro)
      find_fgbg(datapath, name)

  return

def find_fgbg(datapath, name):
  intxns = dict()
  with open(datapath.path) as f:
    for i, line in enumerate(f):
      pass
      if i % 6000000 == 0:
        print i, datetime.datetime.now()
      if i > LIMIT:
        break
      words = line.split()
      intxns[' '.join(words[:2])] = float(words[2])

  ordered_keys = sorted(intxns, key = intxns.get)
  # high = ordered_keys[-_NUM:]
  # low = ordered_keys[:_NUM]
  high, low = filter_match_distances(intxns, ordered_keys)
  train, test = training_test_split(high, low)
  train, test = augment(train, test, intxns)

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
  
def training_test_split(high, low):
  random.shuffle(high)
  random.shuffle(low)
  test = []
  test += high[ : len(high) * TEST_FRAC]
  test += low[ : len(low) * TEST_FRAC]
  train = []
  train += high[len(high) * TEST_FRAC : ]
  train += low[ len(low) * TEST_FRAC : ]
  return train, test

def convert_to_loci(keys):
  loci = []
  for itx in keys:
    loci += itx.split()
  return loci

def loci_match(itx, loci):
  # itx is a space-delimited string with 2 components
  # loci is a list of strings
  return itx.split()[0] in loci or itx.split()[1] in loci

def augment(train, test, intxns, ordered_keys):
  train_loci = convert_to_loci(train)
  test_loci = convert_to_loci(test)

  add_to_train = []
  add_to_test = []

  for i in range(len(ordered_keys)):
    ok = ordered_keys[i]
    if intxns[ok] < FG_MIN:
      break
    if loci_match(ok, train_loci) and not loci_match(ok, test_loci):
      add_to_train.append(ok)
    if loci_match(ok, test_loci) and not loci_match(ok, train_loci):
      add_to_test.append(ok)

  for i in range(len(ordered_keys) - 1, 0, -1):
    ok = ordered_keys[i]
    if intxns[ok] < BG_MIN:
      continue
    if intxns[ok] > BG_MAX:
      break
    if loci_match(ok, train_loci) and not loci_match(ok, test_loci):
      add_to_train.append(ok)
    if loci_match(ok, test_loci) and not loci_match(ok, train_loci):
      add_to_test.append(ok)

  # It is possible that we added two intxns that link fg/bg
  # Example: Train is 1, test is 2. We have a new intxn 1-3 and 2-3.
  # To avoid this, we remove overlap among the augmented sets
  cleanup_loci = convert_to_loci(add_to_train)
  add_to_test = [s for s in add_to_test if not loci_match(s, cleanup_loci)]

  train += add_to_train
  test += add_to_test

  return train, test

def get_dist(key):
  return abs(int(key.split()[1]) - int(key.split()[0]))

def filter_match_distances(intxns, ordered_keys, datapath):
  # Finds the top _NUM intxns within MAX_DIST and matches (if possible) 1-to-1
  # between foreground and background

  # Store keys in high, and store their dists
  high = []
  dists = []
  high_loci = []
  for i in range(len(ordered_keys) - 1, 0, -1):
    ok = ordered_keys[i]
    dist = get_dist(ok)
    if dist <= MAX_DIST and not loci_match(ok, high_loci):
      high.append(ok)
      dists.append(dist)
      high_loci += ok.split()
    if len(high) == _NUM:
      break
    if intxns[ok] < FG_MIN:
      print 'ERROR: Insufficient number of foreground', CELLTYPE, CHRO, '_NUM:', _NUM, 'found:', len(high), 'processed:', len(ordered_keys) - i
      print 'Suggested to restart fgbg with lower _NUM parameter'
      break
      # sys.exit(0)

  low = []
  low_loci = []
  for i in range(len(ordered_keys)):
    ok = ordered_keys[i]
    dist = get_dist(ok)
    if intxns[ok] < BG_MIN:
      continue
    if dist in dists and not loci_match(ok, high_loci) and not loci_match(ok, low_loci):
      # Find an intxn that is unique to all fg/bg intxns so far
      # matches a distance in foreground 1-to-1
      dists.remove(dist)
      low.append(ok)
      low_loci += ok.split()
    if len(dists) == 0:
      # Ensures we meet _NUM intxns      
      break
    if intxns[ok] > BG_MAX:
      print 'ERROR: Insufficient number of background', CELLTYPE, CHRO, '_NUM:', len(high), 'found:', len(low), 'processed:', i
      print 'Suggested to restart fgbg with lower _NUM parameter'
      break
      # sys.exit(0)
  return high, low


if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start