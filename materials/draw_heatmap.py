import pylab
import numpy
import sys

def main(argv):
  if len(argv) != 1:
    print 'usage: matrix_filename'
    print
  else:
    fname = argv[0]
    mat = numpy.array([ [ float(val) for val in line.split() ] for line in open(fname).readlines() ])
    pylab.imshow(mat)
    pylab.ylabel('Sequence A')
    pylab.xlabel('Sequence B')
    pylab.colorbar()
    pylab.show()

if __name__ == '__main__':
  main(sys.argv[1:])
