import random
import sys
sys.setrecursionlimit(20000)

# class used as function "decorator":
class Memoized:
  def __init__(self, function):
    self._function = function
    self._cache = {}
  def __call__(self, *args):
    if args not in self._cache:
      # not in the cache: call the function and store the result in
      # the cache
      self._cache[args] = self._function(*args)
    # the result must now be in the cache:
    return self._cache[args]

comp = {'G':'C', 'A':'U', 'U':'A', 'C':'G'}  

def score_base_pair(nuc_a, nuc_b):
  if nuc_a == comp[nuc_b]:
    return 1
  return 0

# i and j are the start and end indices (inclusive)
@Memoized
def score_rna_structure(seq,i,j):
  if j < i:
    return 0

  if j == i:
    return 0

  # case 1: peel off and don't use left base
  case_1 = score_rna_structure(seq,i+1,j)

  # case 2: peel off and don't use right base
  case_2 = score_rna_structure(seq,i,j-1)

  # case 3: peel off left and right bases, form base pair
  case_3 = score_rna_structure(seq,i+1,j-1) + score_base_pair(seq[i],seq[j])

  if j-i>1:
    # case 4: for any k, form two sub-structures:
    # one using {i,i+1,...k} and another using {k+1,k+2,...j}
    case_4 = max( [ score_rna_structure(seq,i,k) + score_rna_structure(seq,k+1,j) for k in xrange(i+1, j) ] )
    return max(case_1, case_2, case_3, case_4)
  return max(case_1, case_2, case_3)

@Memoized
def backtrack(seq,i,j):
  result = ''

  if j < i:
    return ''
  if j == i:
    # this character has no one left to pair with; peel it off:
    return ' '

  # case 1: peel off and don't use left base
  case_1 = score_rna_structure(seq,i+1,j)
  if score_rna_structure(seq,i,j) == case_1:
    return ' ' + backtrack(seq,i+1,j)

  # case 2: peel off and don't use right base
  case_2 = score_rna_structure(seq,i,j-1)
  if score_rna_structure(seq,i,j) == case_2:
    return backtrack(seq,i,j-1) + ' '

  # case 3: peel off left and right bases, form base pair
  case_3 = score_rna_structure(seq,i+1,j-1) + score_base_pair(seq[i],seq[j])
  if score_rna_structure(seq,i,j) == case_3:
    return '(' + backtrack(seq,i+1,j-1) + ')'

  # case 1-3 were not used, so case 4 must be used; therefore, j-i>1:

  # case 4: for any k, form two sub-structures:
  # one using {i,i+1,...k} and another using {k+1,k+2,...j}
  for k in range(i+1,j):
    case_4_k = score_rna_structure(seq,i,k) + score_rna_structure(seq,k+1,j)
    if score_rna_structure(seq,i,j) == case_4_k:
      return backtrack(seq,i,k) + backtrack(seq,k+1,j)

  # it is impossible to run code down here; we must have returned

def get_structure(seq):
  n = len(seq)
  best_score,index = max( [ (score_rna_structure(seq,i,j),(i,j)) for i in xrange(n) for j in xrange(n) ] )
  print 'MATCHES:', best_score
  print seq
  print backtrack(seq,0,n-1)

def random_seq(n):
  GATC=['G','A','U','C']

  ints = [ random.randint(0,3) for i in xrange(n) ]
  a = [ GATC[i] for i in ints ]

  res = ''
  for c in a:
    res += c
  return res

def main(argv):
  if len(argv) == 1:
    seq=argv[0]
    seq = seq.upper()
    seq = seq.replace('T','U')

    get_structure(seq)

    print
    print 'scores of random sequences with same length:'
    for i in xrange(16):
      print score_rna_structure(random_seq(len(seq)),0,len(seq)-1)
  else:
    print 'usage: <dna_or_rna_seq>'
    print 'e.g., try running python nussinov.py GCCCGGAUGAUCCUCAGUGGUCUGGGGUGCAGGCUUCAAACCUGUAGCUGUCUAGCGACAGAGUGGUUCAAUUCCACCUUUCGGGCG'
    print '(tRNA-SeC-TCA-1-1 from http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi19/genes/tRNA-SeC-TCA-1-1.html)'

if __name__=='__main__':
  main(sys.argv[1:])
