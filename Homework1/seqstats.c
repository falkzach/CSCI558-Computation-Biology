#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

static const char INPUT_FILE[] = "sequences.fa";

FILE * inputFile;
char * sequence;
int * sequenceIds;
char ** sequences;

char * line = NULL;
size_t len = 0;
ssize_t read;

int readFile() { //TODO: pass in pointers for sequences arrays
  int sequences;
  int sequenceId;

  inputFile = fopen(INPUT_FILE, "r");
  if (inputFile == NULL) {
    perror("file not found!\n");
    exit(EXIT_FAILURE);
  }
  while((read = getline(&line, &len, inputFile)) != -1) { //TODO: getline is part of POXIX, won't compile on windows, start working on the server i guess
    printf("line: %s\n", line);
  }
  fclose(inputFile);
  if (line) {
      free(line);
  }

  //TODO: mallock memory for sequences and sequenceIds

  inputFile = fopen(INPUT_FILE, "r");
  if (inputFile == NULL) {
    perror("file not found!\n");
    exit(EXIT_FAILURE);
  }
  while((read = getline(&line, &len, inputFile)) != -1) {
    //TODO: read in sequences
  }
  if (line) {
      free(line);
  }
  printf("hello?");
  return 0;
}

int main(int argc, char **argv) {
  printf("Hello C World\n");
  readFile();
  return 0;
}
