CC=g++
CFLAGS=-I -O3 -std=c++11 -march=native -fopenmp -Wall
TARGET=find_rna
TEST_INPUT=pyrococcus_furiosus.genome
TEST_INPUT=small.txt

all: $(TARGET).o

$(TARGET).o: $(TARGET).cpp
	$(CC) $(CFLAGS) -o $(TARGET).o $(TARGET).cpp

clean:
	rm -f *.o 
	
check: all
	./$(TARGET).o $(TEST_INPUT)
