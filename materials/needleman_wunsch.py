from smith_waterman import *

@Memoized
def nw_get_best_score_ending_here(A, B, index_a, index_b):
  if index_a == 0 and index_b == 0:
    return 0

  nuc_a, nuc_b = A[index_a-1], B[index_b-1]

  best_scores_ending_before_here = []
  if index_a > 0:
    best_scores_ending_before_here.append(nw_get_best_score_ending_here(A, B, index_a-1, index_b) + score_pair(nuc_a, '-'))
  if index_b > 0:
    best_scores_ending_before_here.append(nw_get_best_score_ending_here(A, B, index_a, index_b-1) + score_pair('-', nuc_b))
  if index_a > 0 and index_b > 0:
    best_scores_ending_before_here.append(nw_get_best_score_ending_here(A, B, index_a-1, index_b-1) + score_pair(nuc_a, nuc_b))

  return max(best_scores_ending_before_here)

@Memoized
def nw_get_best_score_starting_here(A, B, index_a, index_b):
  return nw_get_best_score_ending_here(A[::-1], B[::-1], len(A)-index_a, len(B)-index_b)

@Memoized
def nw_get_best_score_pasing_through(A, B, index_a, index_b):
  return nw_get_best_score_starting_here(A, B, index_a, index_b) + nw_get_best_score_ending_here(A, B, index_a, index_b)

# backtrack:
def nw_backtrack(A, B, index_a, index_b):
  score_ending_here = nw_get_best_score_ending_here(A, B, index_a, index_b)

  my_path = [(index_a, index_b)]

  if index_a == 0 and index_b == 0:
    return my_path

  nuc_a, nuc_b = A[index_a-1], B[index_b-1]
   
  # Note: beware of floating point arithmetic ==:
  if index_a > 0 and index_b > 0:
    if nw_get_best_score_ending_here(A, B, index_a-1, index_b-1) + score_pair(nuc_a, nuc_b) == score_ending_here:
      return my_path + nw_backtrack(A, B, index_a-1, index_b-1)
  if index_a > 0:
    if nw_get_best_score_ending_here(A, B, index_a-1, index_b) + score_pair(nuc_a, '-') == score_ending_here:
      return my_path + nw_backtrack(A, B, index_a-1, index_b)
  if index_b > 0:
    if nw_get_best_score_ending_here(A, B, index_a, index_b-1) + score_pair('-', nuc_b) == score_ending_here:
      return my_path + nw_backtrack(A, B, index_a, index_b-1)
  else:
    return my_path

if __name__ == '__main__':
  seq_a = 'GATACAA'
  seq_b = 'GATACCAA'

  print 'Starting here:'
  print_matrix(seq_a, seq_b, nw_get_best_score_starting_here)
  print

  print 'Ending here:'
  print_matrix(seq_a, seq_b, nw_get_best_score_ending_here)
  print

  print 'Passing through:'
  print_matrix(seq_a, seq_b, nw_get_best_score_pasing_through)
  print

  best_index = (len(seq_a), len(seq_b))
  best_score = nw_get_best_score_ending_here(seq_a, seq_b, *best_index)

  print 'best score', best_score, 'from best index', best_index

  best_path = nw_backtrack(seq_a, seq_b, *best_index)[::-1]
  print 'best path', best_path

  print_alignment_from_path(seq_a, seq_b, best_path)

