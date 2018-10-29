from nussinov import *

inf = 1e20

@Memoized
def zuker_energy(seq,i,j):
  # ignore short structures:
  if i<j and i>j-4:
    return inf

  if j<=i:
    return 0

  case_1 = zuker_energy(seq,i+1,j)
  case_2 = zuker_energy(seq,i,j-1)
  case_3 = zuker_energy_i_and_j_maybe_bonded(seq,i,j)

  if j-i>1:
    case_4 = min( [ zuker_energy(seq,i,k) + zuker_energy(seq,k+1,j) for k in xrange(i+1, j) ] )
    return min(case_1, case_2, case_3, case_4)
  return min(case_1, case_2, case_3)

# energies from https://ab.inf.uni-tuebingen.de/teaching/ws06/albi1/script/rna_structure_06Dec2006.pdf
def hairpin_energy(seq,i,j):
  # hairpin
  n = j-i+1
  if n<3:
    return inf

  if n == 3:
    return 4.1
  if n == 4:
    return 4.9
  if n == 5:
    return 4.4
  if n <= 10:
    return 5.3
  if n <= 20:
    return 6.1
  if n <= 25:
    return 6.3
  if n <= 30:
    return 6.5

  return 7

def single_half_of_bulge_or_interior(n):
  if n <= 1:
    return 3.9
  if n <= 2:
    return 3.1
  if n <= 3:
    return 3.5
  if n <= 4:
    return 4.2
  if n <= 5:
    return 4.8
  if n <= 10:
    return 5.
  if n <= 15:
    return 6.0
  if n <= 20:
    return 6.3
  if n <= 25:
    return 6.5
  if n <= 30:
    return 6.7

  return 7

def bulge_or_interior_loop_energy(seq,i,j,i2,j2):
  # bulge loop or interior loop

  # Note: we could customize this energy for these two cases:
  if i2-i>1 and j-j2>1:
    # interior loop
    pass
  else:
    # bulge loop
    pass

  return single_half_of_bulge_or_interior(i2-i) + single_half_of_bulge_or_interior(j-j2)

# stacking
CG=('C','G')
GC=('G','C')
GU=('G','U')
UG=('U','G')
AU=('A','U')
UA=('U','A')
  
pairs = [CG,GC,GU,UG,AU,UA]

stacking_energy_matrix = [[-2.4,-3.3,-2.1,-1.4,-2.1,-2.1],[-3.3,-3.4,-2.5,-1.5,-2.2,-2.4],[-2.1,-2.5,1.3,-0.5,-1.4,-1.3],[-1.4,-1.5,-0.5,0.3,-0.6,-1.0],[-2.1,-2.2,-1.4,-0.6,-1.1,-0.9],[-2.1,-2.4,-1.3,-1.0,-0.9,-1.3]]

paired_stacking_energies = {}
for i in xrange(len(stacking_energy_matrix)):
  row_dict = {}
  for j in xrange(len(stacking_energy_matrix[0])):
    row_dict[pairs[j]] = stacking_energy_matrix[i][j]
  paired_stacking_energies[pairs[i]] = row_dict

single_stacking_energies = {}
for row,p in zip(stacking_energy_matrix,pairs):
  single_stacking_energies[p] = sum(row) / float(len(row))

def stacking_energy(seq,i,j):
  nuc_i,nuc_j = seq[i],seq[j]

  # paired stacking energies:
  if i<j-1:
    nuc_i_p1,nuc_j_m1 = seq[i+1],seq[j-1]
    if (nuc_i,nuc_j) in paired_stacking_energies:
      row_dict = paired_stacking_energies[(nuc_i,nuc_j)]
      if (nuc_i_p1,nuc_j_m1) in row_dict:
        return row_dict[(nuc_i_p1,nuc_j_m1)]

  # single stacking energies:
  if (nuc_i,nuc_j) in single_stacking_energies:
    return single_stacking_energies[(nuc_i,nuc_j)]

  return inf

@Memoized
def zuker_energy_i_and_j_maybe_bonded(seq,i,j):
  # ignore short structures:
  if i<j and i>j-4:
    return inf

  if j<=i:
    return 0
    
  case_1 = hairpin_energy(seq,i,j)
  case_2 = zuker_energy_i_and_j_maybe_bonded(seq,i+1,j-1) + stacking_energy(seq,i,j)
  case_3 = best_bulge_or_interior_loop_energy(seq,i,j)
  if j-i>2:
    case_4 = best_multiloop_energy(seq,i,j)
    return min(case_1, case_2, case_3, case_4)
  return min(case_1, case_2, case_3)

@Memoized
def best_bulge_or_interior_loop_energy(seq,i,j):
  return min([ bulge_or_interior_loop_energy(seq,i,j,i2,j2)+zuker_energy_i_and_j_maybe_bonded(seq,i2,j2) for i2 in xrange(i+1,j-1) for j2 in xrange(i2+1,j) if i2-i + j-j2>2 ])

multiloop_penalty = 5
@Memoized
def best_multiloop_energy(seq,i,j):
  return min([ zuker_energy(seq,i+1,k)+zuker_energy(seq,k+1,j-1) for k in xrange(i+1,j-1) ]) + multiloop_penalty

@Memoized
def zuker_backtrack(seq,i,j):
  # ignore short structures:
  if i<j and i>j-4:
    return ''

  if j<i:
    return ''

  if j==i:
    return ' '

  w_score = zuker_energy(seq,i,j)

  case_1 = zuker_energy(seq,i+1,j)
  if w_score == case_1:
    return ' ' + zuker_backtrack(seq,i+1,j)

  case_2 = zuker_energy(seq,i,j-1)
  if w_score == case_2:
    return zuker_backtrack(seq,i,j-1) + ' '

  case_3 = zuker_energy_i_and_j_maybe_bonded(seq,i,j)
  if w_score == case_3:
    return zuker_backtrack_maybe_bonded(seq,i,j)

  if j-i>1:
    case_4,k = min( [ (zuker_energy(seq,i,k)+zuker_energy(seq,k+1,j), k) for k in xrange(i+1, j) ] )
    if w_score == case_4:
      return zuker_backtrack(seq,i,k) + zuker_backtrack(seq,k+1,j)

@Memoized
def zuker_backtrack_best_bulge_or_interior_loop(seq,i,j):
  bi_score = best_bulge_or_interior_loop_energy(seq,i,j)
  for i2 in xrange(i+1,j-1):
    for j2 in xrange(i2+1,j):
      if i2-i + j-j2>2:
        if bi_score == bulge_or_interior_loop_energy(seq,i,j,i2,j2)+zuker_energy_i_and_j_maybe_bonded(seq,i2,j2):
          return (' '*(i2-i)) + zuker_backtrack_v(seq,i2,j2) + (' '*(j-j2))

@Memoized
def zuker_backtrack_multiloop(seq,i,j):
  vm_score = best_multiloop_energy(seq,i,j)
  for k in xrange(i+1,j-1):
    if vm_score == zuker_energy(seq,i+1,k)+zuker_energy(seq,k+1,j-1)+multiloop_penalty:
      return ' ' + zuker_backtrack(seq,i+1,k) + zuker_backtrack(seq,k+1,j-1) + ' '

@Memoized
def zuker_backtrack_maybe_bonded(seq,i,j):
  # ignore short structures:
  if i<j and i>j-4:
    return ''

  if j<i:
    return ''

  if j==i:
    return ' '

  v_score = zuker_energy_i_and_j_maybe_bonded(seq,i,j)
    
  case_1 = hairpin_energy(seq,i,j)
  if v_score == case_1:
    return ' '*(j-i+1)
  case_2 = zuker_energy_i_and_j_maybe_bonded(seq,i+1,j-1) + stacking_energy(seq,i,j)
  if v_score == case_2:
    return '(' + zuker_backtrack_maybe_bonded(seq,i+1,j-1) + ')'
  case_3 = best_bulge_or_interior_loop_energy(seq,i,j)
  if v_score == case_3:
    return zuker_backtrack_best_bulge_or_interior_loop(seq,i,j)
  if j-i>2:
    case_4 = best_multiloop_energy(seq,i,j)
    if v_score == case_4:
      return zuker_backtrack_multiloop(seq,i,j)

def main(argv):
  if len(argv) == 1:
    seq=argv[0]
    seq = seq.upper()
    seq = seq.replace('T','U')

    print 'FREE ENERGY:', zuker_energy(seq,0,len(seq)-1)

    print seq
    print zuker_backtrack(seq,0,len(seq)-1)

    print 
    print 'scores of random sequences with same length:'
    for i in xrange(16):
      print zuker_energy(random_seq(len(seq)),0,len(seq)-1)

  else:
    print 'usage: <dna_or_rna_seq>'
    print 'e.g., try running python zuker.py GCCCGGAUGAUCCUCAGUGGUCUGGGGUGCAGGCUUCAAACCUGUAGCUGUCUAGCGACAGAGUGGUUCAAUUCCACCUUUCGGGCG'
    print '(tRNA-SeC-TCA-1-1 from http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi19/genes/tRNA-SeC-TCA-1-1.html)'

if __name__=='__main__':
  main(sys.argv[1:])
