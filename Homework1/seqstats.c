#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <regex.h>

static const char INPUT_FILE[] = "sequences.fa"; //TODO: get from args

//file reading
FILE * inputFile;
char * line = NULL;
size_t len = 0;
ssize_t read;

//regex
char * fastaSequenceStartRegex =  ">(.*)\n";
regex_t sequenceRegex;
int reti;
char msgbuf[100];

//sequences
size_t * sequenceLengths;
size_t numberOfSequences = 0;

char * sequence;
char * sequenceIds;
char ** sequences;



int readFile() { //TODO: pass in pointers for sequences arrays
  int sequences;
  int sequenceId;

  /* Compile regular expression */
  reti = regcomp(&sequenceRegex, fastaSequenceStartRegex, 0);
  if (reti) {
      perror("could not compile regex\n");
      exit(1);
  }

  inputFile = fopen(INPUT_FILE, "r");
  if (inputFile == NULL) {
    perror("file not found!\n");
    exit(EXIT_FAILURE);
  }
  while((read = getline(&line, &len, inputFile)) != -1) {

    reti = regexec(&sequenceRegex, line, 0, NULL, 0);
    if (!reti) {
        printf("SEQUENCE START!\n");
        printf("line: %s\n", line);
    }
  }
  fclose(inputFile);

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
  /* Free memory allocated to the pattern buffer by regcomp() */
  regfree(&sequenceRegex);
  return 0;
}

int main(int argc, char **argv) {
  readFile();
  return 0;
}
