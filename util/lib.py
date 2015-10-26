class DataPath():
  # Object for processed Hi-C data paths
  def __init__(self, celltype, chro):
    self.celltype = celltype
    self.chro = str(chro)
    self.path = '/broad/compbio/maxwshen/data/processed_hic/' + self.celltype + '.5kb.mapqg0/chr' + self.chro + '_5kb.RAWobservednorm.txt'

  def __repr__(self):
    return 'Celltype: ' + self.celltype, '\nChr. ' + str(chro) + '\nPath: ' + path