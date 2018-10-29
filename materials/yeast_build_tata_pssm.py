from FASTA import *
import numpy
import pylab as P
P.ion()

nucleotides = ['G', 'A', 'T', 'C']

nucleotide_to_index = {}
for i,nuc in enumerate(nucleotides):
  nucleotide_to_index[nuc] = i

# build PSSM on yeast genome:
yeast = FASTA('s_cerevisiae.fasta')

# motif is TATAwxyzuv
motif_start = 'TATA'
motif_length=10

pseudo_count = 1
count_pssm = numpy.zeros( (motif_length, 4) ) + 1

num_matches = 0
for chromosome_name, chromosome_sequence in yeast.accession_to_sequence.items():
  print 'processing', chromosome_name
  for i in xrange(len(chromosome_sequence)-motif_length):
    sl = chromosome_sequence[i:i+motif_length]
    if sl.startswith(motif_start):
      num_matches += 1
      for i,nuc in enumerate(sl):
        nuc_index = nucleotide_to_index[nuc]
        count_pssm[i][nuc_index] += 1

print
for nuc in nucleotides:
  print '   {:8}'.format(nuc),
print
pssm = numpy.array(count_pssm, float) / (num_matches + pseudo_count)
print pssm

print '[press enter]'
raw_input()

# score using this PSSM
fig_num = 1
for chromosome_name, chromosome_sequence in yeast.accession_to_sequence.items():
  print 'processing', chromosome_name
  window_scores = []
  for i in xrange(len(chromosome_sequence)-motif_length):
    sl = chromosome_sequence[i:i+motif_length]

    # score sl against the pssm:
    motif_prob = 1.0
    for j in xrange(motif_length):
      nuc = sl[j]
      nuc_prob = pssm[j][ nucleotide_to_index[nuc] ]
      motif_prob *= nuc_prob

    if motif_prob > 0.002:
      print '\t', i, sl
    
    window_scores.append(motif_prob)

  P.figure(fig_num)
  P.xlabel('bp')
  P.ylabel('TATA... probability')
  fig_num += 1
  P.plot(window_scores)
  P.title(chromosome_name)

  print '[press enter]'
  raw_input()


