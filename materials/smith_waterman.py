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

gap_penalty = -4
mismatch_penalty = -2
match_bonus = 1

nucleotides = ['G', 'A', 'T', 'C']

match_scores = {}
for nuc in nucleotides:
  match_scores[nuc] = match_bonus

def score_pair(nuc_a, nuc_b):
  if nuc_a == '-' or nuc_b == '-':
    return gap_penalty
  if nuc_a != nuc_b:
    return mismatch_penalty
  return match_bonus

@Memoized
def get_best_score_ending_here(A, B, index_a, index_b):
  if index_a == 0 or index_b == 0:
    return 0

  nuc_a, nuc_b = A[index_a-1], B[index_b-1]

  # always include case to start at this cell
  best_scores_ending_before_here = [0]

  best_scores_ending_before_here.append(get_best_score_ending_here(A, B, index_a-1, index_b) + score_pair(nuc_a, '-'))
  best_scores_ending_before_here.append(get_best_score_ending_here(A, B, index_a, index_b-1) + score_pair('-', nuc_b))
  best_scores_ending_before_here.append(get_best_score_ending_here(A, B, index_a-1, index_b-1) + score_pair(nuc_a, nuc_b))

  return max(best_scores_ending_before_here)

@Memoized
def get_best_score_starting_here(A, B, index_a, index_b):
  return get_best_score_ending_here(A[::-1], B[::-1], len(A)-index_a, len(B)-index_b)

@Memoized
def get_best_score_pasing_through(A, B, index_a, index_b):
  return get_best_score_starting_here(A, B, index_a, index_b) + get_best_score_ending_here(A, B, index_a, index_b)

def print_matrix(A, B, score):
  print '   ',
  for j in xrange(len(B)+1):
    if j == 0:
      print ' '*(5-1),
    else:
      print '{:5}'.format(B[j-1]),
  print
  
  for i in xrange(len(A)+1):
    if i != 0:
      print A[i-1],
    else:
      print ' ',
    for j in xrange(len(B)+1):
      print '{:5.2f}'.format(score(A, B, i,j)),
    print

# backtrack:
def backtrack(A, B, index_a, index_b):
  score_ending_here = get_best_score_ending_here(A, B, index_a, index_b)

  my_path = [(index_a, index_b)]

  if index_a == 0 or index_b == 0:
    return my_path

  nuc_a, nuc_b = A[index_a-1], B[index_b-1]
   
  # Note: beware of floating point arithmetic ==:
  if get_best_score_ending_here(A, B, index_a-1, index_b-1) + score_pair(nuc_a, nuc_b) == score_ending_here:
    return my_path + backtrack(A, B, index_a-1, index_b-1)
  elif get_best_score_ending_here(A, B, index_a-1, index_b) + score_pair(nuc_a, '-') == score_ending_here:
    return my_path + backtrack(A, B, index_a-1, index_b)
  elif get_best_score_ending_here(A, B, index_a, index_b-1) + score_pair('-', nuc_b) == score_ending_here:
    return my_path + backtrack(A, B, index_a, index_b-1)
  else: # best path began here (upper, upper-left, and left are not good):
    return my_path

def print_alignment_from_path(A, B, path):
  alignment_a = ''
  alignment_b = ''
  for i in xrange(len(path)-1):
    index_a, index_b = path[i]
    next_index_a, next_index_b = path[i+1]

    if (index_a + 1, index_b + 1) == (next_index_a, next_index_b):
      # diagonal move:
      alignment_a += A[next_index_a-1]
      alignment_b += B[next_index_b-1]
    elif index_a == next_index_a:
      # vertical move
      alignment_a += '-'
      alignment_b += B[next_index_b-1]
    else:
      # horizontal move
      alignment_a += A[next_index_a-1]
      alignment_b += '-'

  print 'alignment'
  print alignment_a
  print alignment_b

if __name__ == '__main__':
  #seq_a = 'TGATCGATCGATCGATCGATCA'
  #seq_b = 'CGATCATCGACGATCCATCT'
  #seq_a = 'TGATGATCGATCGATCGATCA'
  #seq_b = 'CGATCATCGACGATCCATCT'
  #seq_a = 'TGATCGA'
  #seq_b = 'CGATCAT'
  #seq_a = 'GATTACA'
  #seq_b = 'GATACCA'
  seq_a = 'GGGACCC'
  seq_b = 'TTTAGGG'

  print 'Starting here:'
  print_matrix(seq_a, seq_b, get_best_score_starting_here)
  print

  print 'Ending here:'
  print_matrix(seq_a, seq_b, get_best_score_ending_here)
  print

  print 'Passing through:'
  print_matrix(seq_a, seq_b, get_best_score_pasing_through)
  print

  best_score, best_index = max( [ (get_best_score_pasing_through(seq_a, seq_b, i,j),(i,j)) for i in xrange(len(seq_a)+1) for j in xrange(len(seq_b)+1) ] )

  print 'best score', best_score, 'from best index', best_index

  best_path = backtrack(seq_a, seq_b, *best_index)[::-1]
  print 'best path', best_path

  print_alignment_from_path(seq_a, seq_b, best_path)

