CC=g++
CFLAGS=-I -O3 -std=c++11
TARGET=assemble
TEST_INPUT=fragments.txt

all: $(TARGET).o

$(TARGET).o: $(TARGET).cpp
	$(CC) $(CFLAGS) -o $(TARGET).o $(TARGET).cpp

clean:
	rm -f *.o 
	
check: all
	./$(TARGET).o $(TEST_INPUT)
