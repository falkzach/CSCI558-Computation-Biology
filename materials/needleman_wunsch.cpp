#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

float & index(float*matrix, long num_cols, int r, int c) {
  long flat_index = r*num_cols + c;
  return matrix[flat_index];
}

float match_bonus = 1.0;
float mismatch_penalty = -2.0;
float gap_penalty = -4.0;

void init_first_row(float*matrix, long num_cols) {
  for (long j=1; j<num_cols; ++j) {
    matrix[j] = j*gap_penalty;
  }
}

void init_first_column(float*matrix, long num_rows, long num_cols) {
  for (long j=1; j<num_rows; ++j) {
    matrix[j*num_cols] = j*gap_penalty;
  }
}

void stripe_row(const std::string & seq_a, const std::string & seq_b, float*matrix, long num_cols, long row) {
  // From above (gap):
  for (long i=0; i<num_cols; ++i)
    index(matrix, num_cols, row, i) = index(matrix, num_cols, row-1, i) + gap_penalty;

  char nuc_a = seq_a[row-1];

  // From upper-left (match or mismatch):
  for (long i=1; i<num_cols; ++i) {
    char nuc_b = seq_b[i-1];

    float match_or_mismatch;
    // Note: could use table instead of if statement to make much faster:
    if (nuc_a == nuc_b)
      match_or_mismatch = match_bonus;
    else
      match_or_mismatch = mismatch_penalty;

    index(matrix, num_cols, row, i) = std::max(index(matrix, num_cols, row, i), index(matrix, num_cols, row-1, i-1) + match_or_mismatch);
  }

  // From left (gap):
  for (long i=1; i<num_cols; ++i)
    index(matrix, num_cols, row, i) = std::max(index(matrix, num_cols, row, i), index(matrix, num_cols, row, i-1) + gap_penalty);
}

float*create_matrix(const std::string & seq_a, const std::string & seq_b) {
  long n_a = seq_a.size();
  long n_b = seq_b.size();

  long num_rows = n_a+1;
  long num_cols = n_b+1;

  float*matrix = new float[num_rows*num_cols];

  // top-left corner is 0.0
  matrix[0] = 0.0;
  
  // stripe first column going down:
  init_first_row(matrix, num_cols);
  init_first_column(matrix, num_rows, num_cols);

  for (int i=1; i<num_rows; ++i)
    stripe_row(seq_a, seq_b, matrix, num_cols, i);

  return matrix;
}

void needleman_wunsch(std::string seq_a, std::string seq_b) {
  float*ending_mat = create_matrix(seq_a, seq_b);

  std::reverse(seq_a.begin(), seq_a.end());
  std::reverse(seq_b.begin(), seq_b.end());
  float*starting_mat_reversed = create_matrix(seq_a, seq_b);

  for (unsigned long i=0; i<=seq_a.size(); ++i) {
    for (unsigned long j=0; j<=seq_b.size(); ++j)
      //std::cout << index(starting_mat_reversed, seq_b.size()+1, seq_a.size()-i, seq_b.size()-j) + index(ending_mat, seq_b.size()+1, i, j) << "\t";
      std::cout << index(ending_mat, seq_b.size()+1, i, j) << "\t";
    std::cout << std::endl;
  }

  std::cerr << "BEST_SCORE " << ending_mat[(seq_a.size()+1)*(seq_b.size()+1)-1] << std::endl;
}

int main(int argc, char**argv) {
  if (argc <= 2) {
    std::cerr << "usage: <seq_a_fname> <seq_b_fname>" << std::endl << std::endl;
  }
  else {
    std::string seq_a;
    std::ifstream fin_a(argv[1]);
    fin_a >> seq_a;

    std::string seq_b;
    std::ifstream fin_b(argv[2]);
    fin_b >> seq_b;

    needleman_wunsch(seq_a, seq_b);
  }
  return 0;
}
