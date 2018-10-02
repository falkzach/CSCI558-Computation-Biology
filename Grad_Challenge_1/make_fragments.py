import sys
import numpy

def main(argv):
  if len(argv) != 4:
    print 'usage: <genome> <number_of_reads> <min_read_length> max_read_length>'
    print
  else:
    genome, reads, min_read_length, max_read_length = argv
    reads, min_read_length, max_read_length = [ int(x) for x in (reads, min_read_length, max_read_length) ]
    n = len(genome)

    for i in xrange(reads):
      start = numpy.random.randint(n - min_read_length)
      length = numpy.random.randint(min_read_length, max_read_length+1)
      print genome[start : start+length]

if __name__=='__main__':
  main(sys.argv[1:])
