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
  global TEST_FRAC
  global BG_MIN
  global BG_MAX
  global FG_MIN
  global CHRO
  global CELLTYPE
  global AUGMENT_LIMIT
  _NUM = 500             # Find this many of each set (foreground, background)
  MAX_DIST = 1500000      # maximum intxn distance to pull
  LIMIT = float('inf')    # No early stopping
  FG_MIN = 10
  BG_MIN = 0.5
  BG_MAX = 1.5
  TEST_FRAC = 0.1         # Percent of test split
  AUGMENT_LIMIT = 40       # Limit on augmentation-fold for intxns

  # Using NATO phonetic alphabet
  name = sys.argv[1]
  OUT_PATH = OUT_PATH + name + '/'

  # celltypes = ['IMR90', 'GM12878', 'K562']
  # No chr9 (missing data), chr X
  chrs = ['22', '21', '20', '19', '18', '17', '16', '15', '14', \
          '13', '12', '11', '10', '8', '7', '6', '5', \
          '4', '3', '2', '1']
  celltypes = ['K562']
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

  ordered_keys = sorted(intxns, key = intxns.get)
  # high = ordered_keys[-_NUM:]
  # low = ordered_keys[:_NUM]
  high, low = filter_match_distances(intxns, ordered_keys)
  train, test = training_test_split(high, low)
  train, test = augment(train, test, intxns, ordered_keys)

  write_file('train', train, intxns)
  write_file('test', test, intxns)
  return
  

def write_file(typ, keys, intxns):
  # type is either 'train' or 'test'
  # Also mirrors the interactions
  ensure_dir_exists(OUT_PATH)
  out_fn = OUT_PATH + typ + '.' + CELLTYPE + '.' + CHRO + '.txt'
  with open(out_fn, 'w') as f:
    f.write('>  ' + ' '.join([CELLTYPE, CHRO, str(LIMIT)]) + '\n')
    for key in keys:
      f.write(key + '\t' + str(intxns[key]) + '\n')
      f.write(key.split()[1] + ' ' + key.split()[0] + '\t' + str(intxns[key]) + '\n')
  print 'Written to', out_fn
  return


def training_test_split(high, low):
  # Randomly grab TEST_FRAC from both high and low
  print '  Splitting into training/testing sets...'
  random.shuffle(high)
  random.shuffle(low)
  test = []
  test += high[ : int(len(high) * TEST_FRAC)]
  test += low[ : int(len(low) * TEST_FRAC)]
  train = []
  train += high[int(len(high) * TEST_FRAC) : ]
  train += low[ int(len(low) * TEST_FRAC) : ]
  return train, test


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

def test_augment(ok, typ, train_loci, test_loci, train_counter, test_counter):
  # Tests to see if intxn 'ok' can enter the augment set for 'typ'
  # In particular, ensures that we don't pass AUGMENT_LIMIT
  passes = False
  if typ == 'train':
    if loci_match(ok, train_loci) and not loci_match(ok, test_loci):
      for locus in ok.split():
        if locus in train_counter:
          if train_counter[locus] > AUGMENT_LIMIT - 1:
            return False
      passes = True
  if typ == 'test':
    if loci_match(ok, test_loci) and not loci_match(ok, train_loci):
      for locus in ok.split():
        if locus in test_counter:
          if test_counter[locus] > AUGMENT_LIMIT - 1:
            return False
      passes = True

  if not passes:
    return False
  if passes:
    if typ == 'train':
      for locus in ok.split():
        if locus in train_counter:
          train_counter[locus] += 1
        else:
          train_counter[locus] = 1
    if typ == 'test':
      for locus in ok.split():
        if locus in test_counter:
          test_counter[locus] += 1
        else:
          test_counter[locus] = 1
  return True 


def augment(train, test, intxns, ordered_keys):
  # Find interactions within foreground and background sets that
  # overlap with interactions in training and test. 
  # Expands training and test sets.
  # Note: There are many loci that are involved in both foreground and background intxns
  #   As a result, for example, foreground only will saturate all 
  #   loci found in both fg/bg.
  print '  Augmenting...'
  count_loci_enrich('trainpre', train)
  count_loci_enrich('testpre', test)

  train_loci = convert_to_loci(train)
  test_loci = convert_to_loci(test)

  # pdb.set_trace()

  add_to_train = []
  add_to_test = []

  print '    Background...'
  train_counter = {key: 0 for key in train_loci}
  test_counter = {key: 0 for key in test_loci}
  for i in range(len(ordered_keys)):
    ok = ordered_keys[i]
    if intxns[ok] < BG_MIN:
      continue
    if intxns[ok] > BG_MAX:
      break
    if test_augment(ok, 'train', train_loci, test_loci, train_counter, test_counter):
      add_to_train.append(ok)
    if test_augment(ok, 'test', train_loci, test_loci, train_counter, test_counter):
      add_to_test.append(ok)
  print '      -BG-  Train:', len(add_to_train), '  Test:', len(add_to_test)
  print '      Train Enrichment (avg, min, max):', np.mean(train_counter.values()), min(train_counter.values()), max(train_counter.values())

  print '    Foreground...'
  train_counter = {key: 0 for key in train_loci}
  test_counter = {key: 0 for key in test_loci}
  for i in range(len(ordered_keys) - 1, 0, -1):
    ok = ordered_keys[i]
    if intxns[ok] < FG_MIN:
      break
    if test_augment(ok, 'train', train_loci, test_loci, train_counter, test_counter):
      add_to_train.append(ok)
    if test_augment(ok, 'test', train_loci, test_loci, train_counter, test_counter):
      add_to_test.append(ok)
  print '      -FG+BG-  Train:', len(add_to_train), '  Test:', len(add_to_test)
  print '      Train Enrichment (avg, min, max):', np.mean(train_counter.values()), min(train_counter.values()), max(train_counter.values())

  # It is possible that we added two intxns that link fg/bg
  # Example: Train is 1, test is 2. We have a new intxn 1-3 and 2-3.
  # To avoid this, we remove overlap among the augmented sets
  # This also guarantees that training and test are completely separate.
  # We remove from test set b/c we prefer training set to be larger.  
  cleanup_loci = convert_to_loci(add_to_train)
  print '    Filtering test augment...'
  add_to_test = [s for s in add_to_test if not loci_match(s, cleanup_loci)]

  print '    Distance matching training...'
  add_to_train = distance_match_augment(add_to_train, intxns)

  print '    Final Train Augmentation:', len(add_to_train)
  print '    Final Test Augmentation:', len(add_to_test)

  count_loci_enrich('train', add_to_train)
  count_loci_enrich('test', add_to_test)
  train += add_to_train
  test += add_to_test

  return train, test


def distance_match_augment(add_to_train, intxns):
  # Distance matches fg/bg in training augment
  fg = []
  bg = []
  keep_fg = []
  keep_bg = []


  for itx in add_to_train:
    if intxns[itx] > BG_MAX:
      fg.append(itx)
    else:
      bg.append(itx)

  fg_dists = defaultdict(list)
  for key in fg:
    fg_dists[get_dist(key)].append(key)

  for key in bg:
    dist = get_dist(key)
    if dist in fg_dists and len(fg_dists[dist]) > 0:
      keep_fg.append(fg_dists[dist][0])
      fg_dists[dist] = fg_dists[dist][1:]
      keep_bg.append(key)

  return keep_fg + keep_bg


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