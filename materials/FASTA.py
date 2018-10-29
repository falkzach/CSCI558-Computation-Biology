class FASTA:
  def __init__(self, fname):
    self.accession_to_sequence = {}

    accession = ''
    seq = ''
    for line in open(fname).readlines():
      line = line.replace('\n', '')
      if line[0] == '>':
        if seq != '':
          self.accession_to_sequence[accession] = seq

        accession = line.split()[0][1:]
        seq = ''

      else:
        seq += line.upper()

    if seq != '':
      self.accession_to_sequence[accession] = seq

if __name__=='__main__':
  FASTA('test.fasta')
