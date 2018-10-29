# The MIT License (MIT)

# Copyright (c) 2018 Oliver Serang

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import sys
import numpy

priors = {'normal':0.996, 'gc_rich':0.004}
emissions = { 'normal':{'G':0.209,'A':0.291,'T':0.291,'C':0.209}, 'gc_rich':{'G':0.331,'A':0.169,'T':0.169,'C':0.331} }
transitions = { 'normal':{'normal':0.999, 'gc_rich':0.001}, 'gc_rich':{'normal':0.01, 'gc_rich':0.99}}

states = ['normal', 'gc_rich']
nucleotides = ['G', 'A', 'T', 'C']

def viterbi():
  best_log_from_left = []
    
  base = genome[0]
  best_log_chances_of_arriving_at_current_layer = {}
  for s in states:
    best_log_chances_of_arriving_at_current_layer[s] = numpy.log(priors[s]) + numpy.log(emissions[s][base])

  best_log_from_left.append(best_log_chances_of_arriving_at_current_layer)

  for i in xrange(n-1):
    base = genome[i+1]

    new_best_log_chances_of_arriving_at_current_layer = {'normal':-numpy.inf, 'gc_rich':-numpy.inf}

    for start_s in states:
      for end_s in states:
        new_best_log_chances_of_arriving_at_current_layer[end_s] = max( new_best_log_chances_of_arriving_at_current_layer[end_s], best_log_chances_of_arriving_at_current_layer[start_s]+numpy.log(transitions[start_s][end_s])+numpy.log(emissions[end_s][base]) )

    best_log_chances_of_arriving_at_current_layer = new_best_log_chances_of_arriving_at_current_layer
    best_log_from_left.append(best_log_chances_of_arriving_at_current_layer)

  viterbi_log_prob,viterbi_end_s = max([ (b,a) for a,b in best_log_chances_of_arriving_at_current_layer.items()])

  # perform backtracking:
  viterbi_path = [viterbi_end_s]
  for i in range(n-1)[::-1]:
    base = genome[i+1]
    log_scores_and_start_states = []
    viterbi_end_s = max([ (best_log_from_left[i][start_s]+numpy.log(transitions[start_s][viterbi_end_s]),start_s) for start_s in states])[1]
    viterbi_path.append(viterbi_end_s)

  viterbi_path = viterbi_path[::-1]
  return viterbi_path, viterbi_log_prob

def print_gc_rich():
  viterbi_path, viterbi_log_prob = viterbi()

  print 'Log probability', viterbi_log_prob
  print 

  for i in xrange(n-1):
    if viterbi_path[i] == 'normal' and viterbi_path[i+1] == 'gc_rich':
      print 'GC-rich region from', i+1,
    if viterbi_path[i] == 'gc_rich' and viterbi_path[i+1] == 'normal':
      print 'to', i
      print

fname = sys.argv[1]
genome = open(fname).read().replace('\n','').upper()
n = len(genome)

print_gc_rich()

