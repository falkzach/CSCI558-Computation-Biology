from smith_waterman import *
import numpy

#seq_a = 'TGATCGATCGATCGATCGATCA'
#seq_b = 'CGATCATCGACGATCCATCT'
seq_a = 'TGATGATCGATCGATCGATCA'
seq_b = 'CGATCATCGACGATCCATCT'
#seq_a = 'TGATCGA'
#seq_b = 'CGATCAT'
#seq_a = 'GATTACA'
#seq_b = 'GATACCA'
#seq_a = 'GGGACCC'
#seq_b = 'TTTAGGG'

def safe_slice(list_or_none, start=None, end=None, gap=None):
  if list_or_none is None:
    return None
  else:
    return list_or_none[start:end:gap]

def advance_column_ending_here_rightward(A, B, prev_column, col_index, upper_seed_row=None):
  col_size = len(A)+1

  next_column = numpy.zeros(col_size)
  if upper_seed_row is not None:
    next_column[0] = upper_seed_row[col_index+1]

  for i in xrange(col_size):
    next_column[i] = max(next_column[i], prev_column[i] + gap_penalty)
    if i > 0:
      next_column[i] = max(next_column[i], next_column[i-1] + gap_penalty)
      next_column[i] = max(next_column[i], prev_column[i-1] + score_pair(A[i-1], B[col_index]) )

  return next_column

def advance_column_starting_here_leftward(A, B, next_column, col_index, lower_seed_row=None):
  return advance_column_ending_here_rightward(A[::-1], B[::-1], next_column[::-1], len(B)-col_index-1, upper_seed_row=safe_slice(lower_seed_row,gap=-1))[::-1]

def get_column_ending_here(A, B, n, upper_seed_row=None, left_seed_col=None):
  assert(left_seed_col is None or len(A)+1 == len(left_seed_col) )
  assert(upper_seed_row is None or len(B)+1 == len(upper_seed_row) )

  assert(n >=0 and n<len(B)+1)

  col_size = len(A)+1

  column = left_seed_col
  if left_seed_col is None:
    column = numpy.zeros(col_size)

  if n == 0:
    return column

  for col in xrange(n):
    column = advance_column_ending_here_rightward(A, B, column, col, upper_seed_row)

  return column

def get_column_starting_here(A, B, n, lower_seed_row=None, right_seed_col=None):
  assert(right_seed_col is None or len(A)+1 == len(right_seed_col) )
  assert(lower_seed_row is None or len(B)+1 == len(lower_seed_row) )

  result = numpy.zeros(len(A)+1)
  return get_column_ending_here(A[::-1], B[::-1], len(B)-n, upper_seed_row=safe_slice(lower_seed_row, gap=-1), left_seed_col=safe_slice(right_seed_col, gap=-1))[::-1]

def get_column_passing_through_here(A, B, n, upper_seed_row=None, left_seed_col=None, lower_seed_row=None, right_seed_col=None):
  return get_column_ending_here(A, B, n, upper_seed_row=upper_seed_row, left_seed_col=left_seed_col) + get_column_starting_here(A, B, n, lower_seed_row=lower_seed_row, right_seed_col=right_seed_col)

def get_row_ending_here(A, B, n, upper_seed_row=None, left_seed_col=None):
  assert(left_seed_col is None or len(A)+1 == len(left_seed_col) )
  assert(upper_seed_row is None or len(B)+1 == len(upper_seed_row) )

  return get_column_ending_here(B, A, n, upper_seed_row=left_seed_col, left_seed_col=upper_seed_row)

def get_row_starting_here(A, B, n, lower_seed_row=None, right_seed_col=None):
  assert(right_seed_col is None or len(A)+1 == len(right_seed_col) )
  assert(lower_seed_row is None or len(B)+1 == len(lower_seed_row) )

  return get_column_starting_here(B, A, n, lower_seed_row=right_seed_col, right_seed_col=lower_seed_row)

def get_row_passing_through_here(A, B, n, upper_seed_row=None, left_seed_col=None, lower_seed_row=None, right_seed_col=None):
  return get_row_ending_here(A, B, n, upper_seed_row=upper_seed_row, left_seed_col=left_seed_col) + get_row_starting_here(A, B, n, lower_seed_row=lower_seed_row, right_seed_col=right_seed_col)

def smith_waterman_linspace_helper(A, B, upper_seed_row=None, left_seed_col=None, lower_seed_row=None, right_seed_col=None, top_row_index=0, left_column_index=0):
  assert(left_seed_col is None or len(A)+1 == len(left_seed_col) )
  assert(upper_seed_row is None or len(B)+1 == len(upper_seed_row) )
  assert(right_seed_col is None or len(A)+1 == len(right_seed_col) )
  assert(lower_seed_row is None or len(B)+1 == len(lower_seed_row) )

  if len(A) == 0:
    return []

  result = []
  if len(B) < 3:
    # <=3 columns in matrix; cannot bisect further

    # store as matrices (this will take linear space since |B| is in O(1))
    matrix_ending_here = numpy.zeros( (len(A)+1, len(B)+1) )
    column = get_column_ending_here(A, B, 0, upper_seed_row, left_seed_col)
    matrix_ending_here.T[0] = column
    for i in xrange(len(B)):
      column = advance_column_ending_here_rightward(A, B, column, i, upper_seed_row)
      matrix_ending_here.T[i+1] = column

    matrix_starting_here = numpy.zeros( (len(A)+1, len(B)+1) )
    column = get_column_starting_here(A, B, len(B), lower_seed_row, right_seed_col)
    matrix_starting_here.T[-1] = column
    for i in range(len(B))[::-1]:
      column = advance_column_starting_here_leftward(A, B, column, i, lower_seed_row)
      matrix_starting_here.T[i] = column

    best_score_passing_through_region = max( [ matrix_ending_here[i,j] + matrix_starting_here[i,j] for i in xrange(len(A)+1) for j in xrange(len(B)+1) ] )

    for i in xrange(len(A)+1):
      for j in xrange(len(B)+1):
        score_passing_through = matrix_ending_here[i,j] + matrix_starting_here[i,j]
        if score_passing_through == best_score_passing_through_region:
          result.append( (score_passing_through, (i+top_row_index, j+left_column_index)) )

    return result

  # in case both start or end both occur in Q1 or Q3:
  start_index = get_best_starting_index(A, B)
  end_index = get_best_ending_index(A, B)
  
  # get the middle column and the row where that middle column has maximal score passing through:
  #mid_column = (len(B)+1)/2
  mid_column = (start_index[1] + end_index[1]) / 2
  column_passing_through = get_column_passing_through_here(A, B, mid_column, upper_seed_row=upper_seed_row, left_seed_col=left_seed_col, lower_seed_row=lower_seed_row, right_seed_col=right_seed_col)
  best_val, best_row = max( [ (v,i) for i,v in enumerate(column_passing_through) ] )
  row_passing_through = get_row_passing_through_here(A, B, best_row, upper_seed_row=upper_seed_row, left_seed_col=left_seed_col, lower_seed_row=lower_seed_row, right_seed_col=right_seed_col)

  column_ending_here = get_column_ending_here(A, B, mid_column, upper_seed_row=upper_seed_row, left_seed_col=left_seed_col)
  column_starting_here = get_column_starting_here(A, B, mid_column, lower_seed_row=lower_seed_row, right_seed_col=right_seed_col)
  row_ending_here = get_row_ending_here(A, B, best_row, upper_seed_row=upper_seed_row, left_seed_col=left_seed_col)
  row_starting_here = get_row_starting_here(A, B, best_row, lower_seed_row=lower_seed_row, right_seed_col=right_seed_col)

  A_sub, B_sub = A[start_index[0]:end_index[0]], B[start_index[1]:end_index[1]]
  #print 'subproblem', A_sub, B_sub

  # process the upper-left sub matrix
  result.extend( smith_waterman_linspace_helper(A[:best_row], B[:mid_column], upper_seed_row=safe_slice(upper_seed_row, end=mid_column+1), left_seed_col=safe_slice(left_seed_col, end=best_row+1), lower_seed_row=row_starting_here[:mid_column+1], right_seed_col=column_starting_here[:best_row+1], top_row_index=top_row_index, left_column_index=left_column_index) )

  # add every cell matching the best score in this column
  best_found = False
  for i,s in list(enumerate(column_passing_through))[::-1]:
    if s == best_row:
      best_found = True
    if best_found and s != best_row:
      break
    if best_found:
      result.append( (s, (i+top_row_index, mid_column+left_column_index)) )

  # process the lower-right sub matrix
  result.extend( smith_waterman_linspace_helper(A[best_row:], B[mid_column:], upper_seed_row=row_ending_here[mid_column:], left_seed_col=column_ending_here[best_row:], lower_seed_row=safe_slice(lower_seed_row, start=mid_column), right_seed_col=safe_slice(right_seed_col, start=best_row), top_row_index=top_row_index+best_row, left_column_index=left_column_index+mid_column) )
  
  return result

def argmax(vec):
  if len(vec) == 0:
    return None

  return max( [ (v,i) for i,v in enumerate(vec) ] )[1]

def get_best_ending_index(A, B):
  best_score = 0
  best_index = (len(A),0)
  
  col_size = len(A)+1
  column = numpy.zeros(col_size)

  for col in range(len(B)):
    column = advance_column_ending_here_rightward(A, B, column, col)

    best_index_for_col = argmax(column)
    if column[best_index_for_col] > best_score:
      best_score = column[best_index_for_col]
      best_index = (best_index_for_col, col+1)

  return best_index

def get_best_starting_index(A, B):
  best_score = 0
  best_index = (len(A),0)
  
  col_size = len(A)+1
  column = numpy.zeros(col_size)

  for col in range(len(B))[::-1]:
    column = advance_column_starting_here_leftward(A, B, column, col)

    best_index_for_col = argmax(column)
    if column[best_index_for_col] > best_score:
      best_score = column[best_index_for_col]
      best_index = (best_index_for_col, col)

  return best_index

def align(A, B):
  nodes = smith_waterman_linspace_helper(A, B)

  print 'best score', max([ x[0] for x in nodes ])
  best_paths = take_maximal_scores_and_remove_duplicates(nodes)
  print 'subset of best paths use', best_paths
  index_set = frozenset(best_paths)
  best_path = linspace_backtrack(A, B, index_set, *best_paths[-1])[::-1]
  print 'single path', best_path
  print_alignment_from_path(A, B, best_path)

def take_maximal_scores_and_remove_duplicates(score_path):
  max_score = max(score_path)[0]
  maximal_scoring_path = [ x for x in score_path if x[0] == max_score ]
  sorted_maximal_path = sorted( [ x[1] for x in maximal_scoring_path ] )

  unique_sorted_maximal_path = []
  node_set = set()
  for n in sorted_maximal_path:
    if n not in node_set:
      unique_sorted_maximal_path.append(n)
    node_set.add(n)

  return unique_sorted_maximal_path

@Memoized
def linspace_get_best_score_ending_here(A, B, index_set, index_a, index_b):
  if (index_a, index_b) not in index_set:
    return 0
  if index_a == 0 or index_b == 0:
    return 0

  nuc_a, nuc_b = A[index_a-1], B[index_b-1]

  # always include case to start at this cell
  best_scores_ending_before_here = [0]

  best_scores_ending_before_here.append(linspace_get_best_score_ending_here(A, B, index_set, index_a-1, index_b) + score_pair(nuc_a, '-'))
  best_scores_ending_before_here.append(linspace_get_best_score_ending_here(A, B, index_set, index_a, index_b-1) + score_pair('-', nuc_b))
  best_scores_ending_before_here.append(linspace_get_best_score_ending_here(A, B, index_set, index_a-1, index_b-1) + score_pair(nuc_a, nuc_b))

  return max(best_scores_ending_before_here)

def linspace_backtrack(A, B, index_set, index_a, index_b):
  score_ending_here = linspace_get_best_score_ending_here(A, B, index_set, index_a, index_b)

  my_path = [(index_a, index_b)]

  if index_a == 0 or index_b == 0:
    return my_path

  nuc_a, nuc_b = A[index_a-1], B[index_b-1]
   
  # Note: beware of floating point arithmetic ==:
  if linspace_get_best_score_ending_here(A, B, index_set, index_a-1, index_b-1) + score_pair(nuc_a, nuc_b) == score_ending_here:
    return my_path + backtrack(A, B, index_a-1, index_b-1)
  elif linspace_get_best_score_ending_here(A, B, index_set, index_a-1, index_b) + score_pair(nuc_a, '-') == score_ending_here:
    return my_path + backtrack(A, B, index_a-1, index_b)
  elif linspace_get_best_score_ending_here(A, B, index_set, index_a, index_b-1) + score_pair('-', nuc_b) == score_ending_here:
    return my_path + backtrack(A, B, index_a, index_b-1)
  else: # best path began here (upper, upper-left, and left are not good):
    return my_path

align(seq_a, seq_b)
