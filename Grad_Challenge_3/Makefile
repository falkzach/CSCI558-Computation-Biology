CC=g++
CFLAGS=-I -O3 -std=c++11 -march=native -fopenmp -Wall
TARGET=cluster_ms
TEST_INPUT=120810_JF_HNU142_5.mgf
# TEST_INPUT=demo758Dacentroid.mgf
# TEST_INPUT=big.mgf 
TEST_OUTPUT=test.txt

all: $(TARGET).o

$(TARGET).o: $(TARGET).cpp
	$(CC) $(CFLAGS) -o $(TARGET).o $(TARGET).cpp

clean:
	rm -f *.o $(TEST_OUTPUT)
	
check: all
	time ./$(TARGET).o $(TEST_INPUT) > $(TEST_OUTPUT)
