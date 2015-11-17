import sys, string, datetime, random, copy, os, commands, fnmatch, re, csv
from collections import defaultdict
import numpy as np

# CLASSES: COMPUTATIONAL BIOLOGY
pass

# CLASSES: MISCELLANEOUS UTILITY
class DataPath():
  # Object for processed Hi-C data paths
  def __init__(self, celltype, chro):
    self.celltype = celltype
    self.chro = str(chro)
    self.path = '/broad/compbio/maxwshen/data/0-PROCESS-HIC/' + self.celltype + '.5kb.mapqg0/chr' + self.chro + '_5kb.RAWobservednorm.txt'

  def __repr__(self):
    return 'Celltype: ' + self.celltype, '\nChr. ' + str(chro) + '\nPath: ' + path

class NarrowPeak():
  def __init__(self, input_str):
    self.start = int(input_str.split()[1])
    self.end = int(input_str.split()[2])
    self.length = self.end - self.start
    self.signal = float(input_str.split()[6])


# METHODS: COMPUTATIONAL BIOLOGY
def read_fasta(fasta_fn):
  # Builds header and reads lists from a fasta file
  headers = []
  reads = []
  get_read = False
  with open(fasta_fn) as f:
    curr_read = ''
    for i, line in enumerate(f):
      if line[0] == '>' or line[0] == '@':
        if curr_read != '':
          reads.append(curr_read)
          curr_read = ''
        headers.append(line.strip())
        get_read = True
      elif get_read:
        curr_read += line.strip().upper()
  if curr_read != '':
    reads.append(curr_read)
  if len(headers) == 0 or len(reads) == 0:
    print 'ERROR: Empty fasta file', fasta_fn
  return headers, reads


def read_fastq(fastq_fn):
  # There is no official single fastq format, so this may not always work
  headers = []
  reads = []
  quals = []
  get_read = False
  get_qual = False
  with open(fastq_fn) as f:
    curr_read = ''
    curr_qual = ''
    for i, line in enumerate(f):
      if line[0] == '@':
        if curr_qual != '':
          quals.append(curr_qual)
          curr_qual = ''
        headers.append(line.strip())
        get_qual = False
        get_read = True
      elif line == '+':
        if curr_read != '':
          reads.append(curr_read)
          curr_read = ''
        get_read = False
        get_qual = True
      else:
        if get_read:
          curr_read += line.strip().upper()
        if get_qual:
          curr_qual += line.strip()
  if curr_read != '':
    reads.append(curr_read)
  if curr_qual != '':
    quals.append(curr_qual)
  if len(headers) == 0 or len(reads) == 0 or len(quals) == 0:
    print 'ERROR: Empty fastq file', fastq_fn
  return headers, reads, qual


def reverse_complement(inp):
  if type(inp) != str:
    print 'Input not a string'
    return ''
  cor = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
  return ''.join([cor[c] for c in inp[::-1]])

# METHODS: MISCELLANEOUS UTILITY
def ensure_dir_exists(directory):
  if not os.path.exists(directory):
    os.makedirs(directory)
  return

def read_delimited_txt(inp_fn, dlm):
  # Reads in a text file with the given delimiter, like '\t'
  print 'Reading in', inp_fn, '...'
  with open(inp_fn) as f:
    reader = csv.reader(f, delimiter = dlm)
    d = list(reader)
  return d

def normalize_0_1(data):
  # Normalizes column-wise so that all values are between 0 and 1
  # Necessary preprocessing for a deep belief network
  data = np.transpose(data)
  print len(data), len(data[0])
  print '  0-1 normalizing data...'
  for i in range(len(data)):
    row = [float(s) for s in data[i]]
    if i % 100 == 0:
      print '  ', i, datetime.datetime.now()
    new_row = []
    (minr, maxr) = (min(row), max(row))
    for j in range(len(row)):
      new_row.append((row[j] - minr)/(maxr - minr))
    data[i] = new_row
  data = np.transpose(data)
  return data