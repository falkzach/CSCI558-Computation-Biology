// The MIT License (MIT)

// Copyright (c) 2018 Oliver Serang

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.



// Two-state GC-rich/non-GC-rich HMM analysis with multistart EM

#include <vector>
#include <iostream>
#include <fstream>
#include <array>
#include <random>
#include "LogDouble.hpp"

class HMMParameters {
private:
  std::array<double, 2> priorProbabilities;
  std::array<std::array<double, 2>, 2> transitionProbabilities;
  std::array<std::array<double, 256>, 2> emissionProbabilities;

public:
  static const std::array<char, 4> emissionCharacters;

  // Prior: Pr(H_0 = state)
  double getPriorProbability(int state) const {
    return priorProbabilities[state];
  }

  // Transition: Pr(H_{i+1} = toState | H_i = fromState)
  double getTransitionProbability(int fromState, int toState) const {
    return transitionProbabilities[fromState][toState];
  }

  // Emission: Pr(D_i = output | H_i = state)
  double getEmissionProbability(int state, char output) const {
    double result = emissionProbabilities[state][(int) output];
    
    // Given the state, the table of emission probabilities for
    // different characters is essentially a map<char, double>.
    // 
    // It's faster to use vector<double>, and use the char as an
    // integer key (i.e. the index).
    if (result >= 0.0)
      return result;

    // A negative probability indicates that an emission for that
    // character has never been initialized, and is thus not
    // supported:
    throw std::exception();
  }

  void enforceEmissionSymmetry() {
    // in every state, make the emission of G == the emission of C (same for A and T)
    for (int state=0; state<2; ++state) {
      double newGC = (emissionProbabilities[state]['G'] + emissionProbabilities[state]['C'])/2.0;
      double newAT = (emissionProbabilities[state]['A'] + emissionProbabilities[state]['T'])/2.0;
      emissionProbabilities[state]['G'] = emissionProbabilities[state]['C'] = newGC;
      emissionProbabilities[state]['A'] = emissionProbabilities[state]['T'] = newAT;
    }
  }

  int getGCRichState() const {
    double state_0 = emissionProbabilities[0]['G'] + emissionProbabilities[0]['C'];
    double state_1 = emissionProbabilities[1]['G'] + emissionProbabilities[1]['C'];
    
    if (state_0 > state_1)
      return 0;
    return 1;
  }

  void setParameters(const std::array<double, 2> & newPriors, const std::array<std::array<double, 2>, 2> & newTransitionProbabilities, const std::array<std::array<double, 256>, 2> & newEmissionProbabilities) {
    for (int state=0; state<2; ++state) {
      priorProbabilities[state] = newPriors[state];
      for (int toState=0; toState<2; ++toState)
	transitionProbabilities[state][toState] = newTransitionProbabilities[state][toState];
      for (int emissionChar=0; emissionChar<256; ++emissionChar)
	emissionProbabilities[state][emissionChar] = newEmissionProbabilities[state][emissionChar];
    }

    enforceEmissionSymmetry();
  }

  HMMParameters() {
    // Initialize all probabilities with -1 so that we can easily
    // detect a query for a disallowed state/emission:

    for (int state=0; state<2; ++state) {
      priorProbabilities[state] = -1.0;
      for (int toState=0; toState<2; ++toState)
	transitionProbabilities[state][toState] = -1.0;
      for (int emissionCharIndex=0; emissionCharIndex<256; ++emissionCharIndex)
	emissionProbabilities[state][emissionCharIndex] = -1.0;
    }
  }

  HMMParameters(const std::string & modelPath) {
    // Initialize all probabilities with -1 so that we can easily
    // detect a query for a disallowed state/emission:

    for (int state=0; state<2; ++state) {
      priorProbabilities[state] = -1.0;
      for (int toState=0; toState<2; ++toState)
	transitionProbabilities[state][toState] = -1.0;
      for (int emissionCharIndex=0; emissionCharIndex<256; ++emissionCharIndex)
	emissionProbabilities[state][emissionCharIndex] = -1.0;
    }
    std::ifstream fin(modelPath);

    double probability, remainingProbability;
    
    // read the prior probabilities (do not read the final one; it
    // will be 1 - the sum of the others)
    remainingProbability = 1.0;
    for (int state=0; state<2-1; ++state) {
      fin >> probability;
      remainingProbability -= probability;
      priorProbabilities[state] = probability;
    }
    priorProbabilities[2-1] = remainingProbability;

    // read the transition probabilities
    for (int fromState=0; fromState<2; ++fromState) {
      remainingProbability = 1.0;
      for (int toState=0; toState<2-1; ++toState) {
	fin >> probability;
	remainingProbability -= probability;
	transitionProbabilities[fromState][toState] = probability;
      }
      transitionProbabilities[fromState][2-1] = remainingProbability;
    }

    // read the emission characters and their probabilities
    for (int state=0; state<2; ++state) {
      remainingProbability = 1.0;
      for (int charIndex=0; charIndex<emissionCharacters.size()-1; ++charIndex) {
	fin >> probability;
	remainingProbability -= probability;
	emissionProbabilities[state][ emissionCharacters[charIndex] ] = probability;
      }
      emissionProbabilities[state][ emissionCharacters[emissionCharacters.size()-1] ] = remainingProbability;
    }
  }

  void randomize() {
    std::uniform_real_distribution<double> uniDist(0.0, 1.0);
    std::random_device rd;
    std::mt19937 gen(rd());

    double totPrior = 0.0;
    for (int state=0; state<2; ++state) {
      double prob = uniDist(gen);
      priorProbabilities[state] = prob;
      totPrior += prob;
      
      double totTrans = 0.0;
      for (int toState=0; toState<2; ++toState) {
	double prob = uniDist(gen);
	transitionProbabilities[state][toState] = prob;
	totTrans += prob;
      }
      for (int toState=0; toState<2; ++toState)
	transitionProbabilities[state][toState] /= totTrans;

      double totEmit = 0.0;
      for (char c : emissionCharacters) {
	double prob = uniDist(gen);
	emissionProbabilities[state][c] = prob;
	totEmit += prob;
      }
      for (char c : emissionCharacters)
	emissionProbabilities[state][c] /= totEmit;
    }

    for (int state=0; state<2; ++state)
      priorProbabilities[state] /= totPrior;

    enforceEmissionSymmetry();
  }

  void print() {
    std::cout << "Priors:" << std::endl;
    for (int state=0; state<2; ++state) {
      std::cout << state << " " << getPriorProbability(state) << "\t";
    }
    std::cout << std::endl;

    std::cout << "Transitions:" << std::endl;
    for (int state=0; state<2; ++state) {
      for (int toState=0; toState<2; ++toState)
	std::cout << state << "-->" << toState << " " << getTransitionProbability(state, toState) << "\t";
      std::cout << std::endl;
    }
    
    std::cout << "Emissions:" << std::endl;
    for (int state=0; state<2; ++state) {
      for (char emissionChar : HMMParameters::emissionCharacters)
	std::cout << state << " :-> " << emissionChar << " " << getEmissionProbability(state, emissionChar) << "\t";
      std::cout << std::endl;
    }
    std::cout << "GC-rich state is " << getGCRichState() << std::endl;

    std::cout << std::endl;
  }

  void save() {
    std::ofstream fout("resultParams.params");
    for (int state=0; state<2-1; ++state) {
      fout << getPriorProbability(state) << "\t";
    }
    fout << std::endl;

    for (int state=0; state<2; ++state) {
      for (int toState=0; toState<2-1; ++toState)
  	fout << getTransitionProbability(state, toState) << "\t";
      fout << std::endl;
    }
    
    for (int state=0; state<2; ++state) {
      for (int charInd=0; charInd < HMMParameters::emissionCharacters.size()-1; ++charInd) {
  	char emissionChar = HMMParameters::emissionCharacters[charInd];
  	fout << getEmissionProbability(state, emissionChar) << "\t";
      }
      fout << std::endl;
    }
  }
};

const std::array<char, 4> HMMParameters::emissionCharacters = {'G', 'A', 'T', 'C'};

class HMM {
private:
  HMMParameters hmmParams;
  std::vector<int> viterbiStates;
  LogDouble viterbiLikelihood;
public:
  const std::string * emissions;
  HMM(const std::string & genome, const HMMParameters & parameters) {
    emissions = & genome;
    hmmParams = parameters;

    std::cout << "constructed HMM with parameters:" << std::endl;
    hmmParams.print();

    viterbiStates = std::vector<int>(emissions->size(), -1); // initialize states to -1
  }

  const HMMParameters & parameters() const {
    return hmmParams;
  }

  int argmax(const std::array<LogDouble, 2> & arr) {
    int maxIndex = -1;
    LogDouble maxVal(-1.0);
    for (int k=0; k<2; ++k)
      if (arr[k] > maxVal) {
	maxVal = arr[k];
	maxIndex = k;
      }
    return maxIndex;
  }

  LogDouble max(const std::array<LogDouble, 2> & arr) {
    return arr[argmax(arr)];
  }

  void estimatePriorsTransitionsAndEmissionsFromViterbi() {
    std::array<double, 2> newPriors;
    std::array<std::array<double, 2>, 2> newTransitions;
    std::array<std::array<double, 256>, 2> newEmissions;

    for (int state=0; state<2; ++state) {
      newPriors[state] = 0.0;
      for (int toState=0; toState<2; ++toState)
	newTransitions[state][toState] = 0.0;
      for (int emissionChar=0; emissionChar<256; ++emissionChar)
	newEmissions[state][emissionChar] = 0.0;
    }

    // Count the events:
    for (int n=0; n<emissions->size(); ++n)
      ++newPriors[ viterbiStates[n] ];

    for (int n=0; n<emissions->size()-1; ++n)
      ++newTransitions[ viterbiStates[n] ][ viterbiStates[n+1] ];

    for (int n=0; n<emissions->size(); ++n)
      ++newEmissions[ viterbiStates[n] ][ (*emissions)[n] ];

    // Normalize the counts:
    for (int state=0; state<2; ++state) {
      newPriors[state] /= viterbiStates.size();
    }

    for (int state=0; state<2; ++state) {
      double totalGivenState = 0.0;
      for (int toState=0; toState<2; ++toState)
	totalGivenState += newTransitions[state][toState];
	
      if (totalGivenState > 0)
	for (int toState=0; toState<2; ++toState)
	  newTransitions[state][toState] /= totalGivenState;
    }

    for (int state=0; state<2; ++state) {
      double totalGivenState = 0.0;
      for (char emissionChar : HMMParameters::emissionCharacters)
	totalGivenState += newEmissions[state][ int(emissionChar) ];

      if (totalGivenState > 0)
	for (char emissionChar : HMMParameters::emissionCharacters)
	  newEmissions[state][int(emissionChar)] /= totalGivenState;
    }

    hmmParams.setParameters(newPriors, newTransitions, newEmissions);
  }

  void emTrain(int numberOfSteps) {
    LogDouble lastViterbiLikelihood(-1.0);
    for (int count=0; count<numberOfSteps; ++count) {
      computeViterbiStates();

      // Note: not completely safe (because of floating point comparison):
      if ( lastViterbiLikelihood == viterbiLikelihood )
	break;
      lastViterbiLikelihood = viterbiLikelihood;

      estimatePriorsTransitionsAndEmissionsFromViterbi();
    }

    // change for the final update of Viterbi states using the newest parameters:
    computeViterbiStates();
    std::cout << "Likelihood after EM " << viterbiLikelihood << std::endl;
  }

  void multiStartEMTrain(int numberStarts, int numberStepsPerEM) {
    std::cerr << "Multi train..." << std::endl;

    // Perform first step to get initial best likelihood:
    emTrain(numberStepsPerEM);
    LogDouble bestLikelihood = viterbiLikelihood;
    HMMParameters bestParams = hmmParams;

    std::cerr << "First likelihood: " << bestLikelihood << std::endl;
    hmmParams.print();

    for (int k=0; k<numberStarts-1; ++k) {
      std::cout << k+1 << "/" << numberStarts << std::endl;
      
      hmmParams.randomize();
      emTrain(numberStepsPerEM);
      if (viterbiLikelihood > bestLikelihood) {
	bestLikelihood = viterbiLikelihood;
	bestParams = hmmParams;
	std::cerr << "New max likelihood: " << bestLikelihood << std::endl;
	hmmParams.print();
      }
    }

    hmmParams = bestParams;
    computeViterbiStates();

    std::cerr << std::endl;
    std::cerr << "Best overall likelihood: " << bestLikelihood << std::endl;
    hmmParams.print();
    hmmParams.save();
  }
  
  const std::vector<int> & computeViterbiStates() {
    std::vector<std::array<LogDouble,2> > maximumProbabilitiesFromLeft(emissions->size());
    std::vector<std::array<int,2> > previousBestAtState(emissions->size());

    // Set prior and emission probailities for the first location:
    for (int k=0; k<2; ++k) {
      maximumProbabilitiesFromLeft[0][k] = LogDouble(hmmParams.getEmissionProbability(k, (*emissions)[0])) * LogDouble(hmmParams.getPriorProbability(k));
    }

    // Initialize maximum probabilities from left:
    for (int n=1; n<emissions->size(); ++n) {
      for (int k=0; k<2; ++k) {

	std::array<LogDouble, 2> newFromLeft;

	for (int fromState=0; fromState<2; ++fromState)
	  newFromLeft[fromState] = LogDouble(hmmParams.getTransitionProbability(fromState, k)) * LogDouble(maximumProbabilitiesFromLeft[n-1][fromState]);

	int bestPrevious = argmax(newFromLeft);
	previousBestAtState[n][k] = bestPrevious;
	maximumProbabilitiesFromLeft[n][k] = LogDouble(hmmParams.getEmissionProbability(k, (*emissions)[n])) * newFromLeft[bestPrevious];
      }
    }

    viterbiLikelihood = max(maximumProbabilitiesFromLeft[emissions->size()-1]);
    
    // Use backtracking to get the Viterbi path:
    int currentState = argmax(maximumProbabilitiesFromLeft[emissions->size()-1]);
    viterbiStates[emissions->size()-1] = currentState;
    for (int n=emissions->size()-1; n>=1; --n) {
      currentState = previousBestAtState[n][currentState];
      viterbiStates[n-1] = currentState;
    }

    return viterbiStates;
  }
};

std::vector<std::pair<unsigned int, unsigned int> > getGCRich(HMM & hmm, int numberStarts, int numberStepsPerEM) {
  hmm.multiStartEMTrain(numberStarts, numberStepsPerEM);

  std::vector<int> viterbiStates = hmm.computeViterbiStates();
  std::vector<std::pair<unsigned int, unsigned int> > contigs;

  int gcRichState = hmm.parameters().getGCRichState();

  bool inContig = false;
  unsigned int start, end;
  for (unsigned int k=0; k<viterbiStates.size(); ++k) {
    if (viterbiStates[k] == gcRichState && ! inContig) {
      start = k;
      end = k;
    }
    if (viterbiStates[k] == gcRichState && inContig)
      end = k;
    
    if (viterbiStates[k] == gcRichState)
      inContig = true;
    else {
      if (inContig) {
	std::cerr << start << " " << end << std::endl;
	contigs.push_back( std::pair<unsigned int, unsigned int>(start, end) );
      }
      inContig = false;
    }
  }

  if (inContig)
    contigs.push_back( std::pair<unsigned int, unsigned int>(start, end) );

  return contigs;
}

std::string loadFile(char*path) {
  std::ifstream fin(path);
  return std::string(std::istreambuf_iterator<char>(fin), std::istreambuf_iterator<char>());
}

int main(int argc, char**argv) {
  if (argc == 4 or argc == 5) {
    std::cerr << "loading..." << std::endl;
    std::string genome = loadFile(argv[1]);

    int numberStarts = atoi(argv[2]);
    int numberIterationsPerEM = atoi(argv[3]);

    std::cerr << "Genome contains " << genome.size() << " characters" << std::endl;
    std::cerr << "processing..." << std::endl;

    HMMParameters parameters;
    if (argc == 5)
      parameters = HMMParameters(argv[4]);
    else
      parameters.randomize();
    
    HMM hmm(genome, parameters);

    std::vector<std::pair<unsigned int, unsigned int> > rich = getGCRich(hmm, numberStarts, numberIterationsPerEM);
    for (const std::pair<unsigned int, unsigned int> & contig : rich) {
      std::cout << contig.first << " " << contig.second << std::endl;
      for (unsigned int k=contig.first; k<=contig.second; ++k)
	std::cout << genome[k];
      std::cout << std::endl << std::endl;
    }
  }
  else
    std::cerr << "usage: GCRichGenomeViterbi <genome> <number_starts> <number_iterations_per_start> [hmm params]" << std::endl;

  return 0;
}
