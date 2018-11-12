CC=g++
CFLAGS=-I -O3 -std=c++11 -march=native -fopenmp -Wall
PYTHON=python3
TARGET=find_rna
TEST_INPUT=pyrococcus_furiosus.genome
TEST_OUTPUT=results.txt
TEST_COMPARE=ranges.txt

all: $(TARGET).o

$(TARGET).o: $(TARGET).cpp
	$(CC) $(CFLAGS) -o $(TARGET).o $(TARGET).cpp

clean:
	rm -f *.o 
	
check: all
	time ./$(TARGET).o $(TEST_INPUT) > $(TEST_OUTPUT)

score: all
	$(PYTHON) score.py $(TEST_COMPARE) $(TEST_OUTPUT)
