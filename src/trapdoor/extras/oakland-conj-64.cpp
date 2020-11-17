// @file
// @author TPOC: contact@palisade-crypto.org
//
// @copyright Copyright (c) 2019, New Jersey Institute of Technology (NJIT)
// All rights reserved.
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// 1. Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution. THIS SOFTWARE IS
// PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
// IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
// MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
// EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#define PROFILE  // define this to enable PROFILELOG and TIC/TOC
// Note must be before all headers

#include <getopt.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <random>
#include "obfuscation/lweconjunctionchcprf.h"
#include "utils/debug.h"
#include "utils/parallel.h"

using namespace lbcrypto;

bool CONJOBF(size_t n_evals, int n);  // defined later

string RandomBooleanString(usint length);

int main(int argc, char* argv[]) {
  PalisadeParallelControls.Enable();

  int opt;  // used in getting options
  int n_evals = 100;
  int n_bits = 10;

  while ((opt = getopt(argc, argv, "e:k:")) != -1) {
    switch (opt) {
      case 'e':
        n_evals = atoi(optarg);
        if (n_evals < 0) n_evals = 0;
        break;
      case 'k':
        n_bits = atoi(optarg);
        if (n_bits < 5) {
          n_bits = 5;
          std::cout << "setting n_bits to minimum size of 5" << std::endl;
        } else if (n_bits >= 13) {
          n_bits = 13;
          std::cout << "setting n_bits to maximum size of 13" << std::endl;
        }
        break;
      case 'h':
      default: /* '?' */
        std::cerr << "Usage: " << argv[0] << " <arguments> " << std::endl
                  << "arguments:" << std::endl
                  << "  -e  number of evaluations (0) {If >3, then all "
                     "evaluations will be random}"
                  << std::endl
                  << "  -k  bitsize of security parameter (ring dimension = "
                     "2^k) [8:13] (10)"
                  << std::endl
                  << "  -h  (false) prints this message" << std::endl;
        exit(EXIT_FAILURE);
    }
  }

  std::cerr << "Running " << argv[0] << " with security parameter "
            << std::to_string(1 << n_bits) << ". Pattern length is 40 bits."
            << std::endl;

  std::cerr << "Running " << argv[0] << " with "
            << ParallelControls::GetNumProcs() << " processors and "
            << ParallelControls().GetNumThreads() << " threads. " << std::endl;

  bool errorflag = CONJOBF(n_evals, 1 << n_bits);

  return (static_cast<int>(errorflag));
}

bool CONJOBF(size_t n_evals, int n) {
  bool errorflag = false;

  TimeVar t;

  TimeVar t1;  // for TIC TOC

  double processingTime(0.0), MSKeyGenTime(0.0), conKeyGenTime(0.0);

  std::string pattern =
      "1?1?10??????10111?1?10??????10111?1?10??????10111?1?10??????1011";
  std::string input1 =
      "1011101110111011101110111011101110111011101110111011101110111011";
  std::string input2 =
      "1011101110111010101110111011101010111011101110101011101110111010";

  size_t len = pattern.length();

  TIC(t);
  LWEConjunctionCHCPRFAlgorithm<DCRTPoly> algorithm(1 << 20, 8, 64, n);
  processingTime = TOC(t);
  std::cout << "Parameter Generation: " << processingTime << "ms" << std::endl;

  std::cout << "n = " << algorithm.GetRingDimension() << std::endl;
  std::cout << "log2 q = " << algorithm.GetLogModulus() << std::endl;

  TIC(t);
  auto key = algorithm.KeyGen();
  MSKeyGenTime = TOC(t);
  std::cout << "Master Secret (Unconstrained) Key Generation: " << MSKeyGenTime
            << "ms" << std::endl;

  TIC(t);
  auto constrainedKey = algorithm.Constrain(key, pattern);
  conKeyGenTime = TOC(t);
  std::cout << "Contstrained Key Generation: " << conKeyGenTime << "ms"
            << std::endl;

  TIC(t);
  const auto value1 = algorithm.Evaluate(key, input1);
  const auto value3 = algorithm.Evaluate(key, input2);
  processingTime = TOC(t);
  std::cout << "Evaluation (unconstrained): 2 * " << processingTime / 2 << "ms"
            << std::endl;
  TIC(t);
  const auto value2 = algorithm.Evaluate(constrainedKey, input1);
  const auto value4 = algorithm.Evaluate(constrainedKey, input2);
  processingTime = TOC(t);
  std::cout << "Evaluation (constrained): 2 * " << processingTime / 2 << "ms"
            << std::endl;
  // std::cout << value1 << std::endl;
  // std::cout << value2 << std::endl;
  std::cout << "pattern: " << pattern << std::endl;
  std::cout << "input 1: " << input1 << std::endl;
  std::cout << (algorithm.EqualTest(value2, value1)
                    ? "Matched (Correct)"
                    : "Did not match (Incorrect)")
            << std::endl;

  if (!algorithm.EqualTest(value2, value1)) errorflag = true;
  // std::cout << value3 << std::endl;
  // std::cout << value4 << std::endl;
  std::cout << "input 2: " << input2 << std::endl;
  std::cout << (algorithm.EqualTest(value4, value3) ? "Matched (Incorrect)"
                                                    : "Did not match (Correct)")
            << std::endl;

  if (algorithm.EqualTest(value4, value3)) errorflag = true;

  if (errorflag) return (errorflag);

  vector<string> inputStr(n_evals);

  for (usint i = 0; i < n_evals; i++) {
    inputStr[i] = RandomBooleanString(len);
  }

  // Variables for timing
  vector<double> timeTokenEval(n_evals);
  vector<double> timeEval(n_evals);
  vector<double> timeEqualTest(n_evals);

  ////////////////////////////////////////////////////////////
  // test the obfuscated pattern
  ////////////////////////////////////////////////////////////
  PROFILELOG("\nEvaluation started");
  for (usint i = 0; i < n_evals; i++) {
    TIC(t1);
    const auto uresult = algorithm.Evaluate(key, inputStr[i]);
    timeTokenEval[i] = TOC_US(t1);

    TIC(t1);
    const auto cresult = algorithm.Evaluate(constrainedKey, inputStr[i]);
    timeEval[i] = TOC_US(t1);

    TIC(t1);
    algorithm.EqualTest(uresult, cresult);
    timeEqualTest[i] = TOC_US(t1);
  }

  // print output timing results
  // note one could use PROFILELOG for these lines
  float sumTokenTime = 0.0;
  float sumTime = 0.0;
  float sumETtime = 0.0;
  for (usint i = 0; i < n_evals; i++) {
    sumTokenTime += timeTokenEval[i];
    std::cout << "T: Token Eval " << i << " execution time:  "
              << "\t" << timeTokenEval[i] / 1000 << " ms" << std::endl;
    sumTime += timeEval[i];
    std::cout << "T: Eval " << i << " execution time:  "
              << "\t" << timeEval[i] / 1000 << " ms" << std::endl;
    sumETtime += timeEqualTest[i];
    std::cout << "T: Equal test " << i << " execution time:  "
              << "\t" << timeEqualTest[i] / 1000 << " ms" << std::endl;
  }
  float aveTime = sumTime / static_cast<float>(n_evals);
  float aveTokenTime = sumTokenTime / static_cast<float>(n_evals);
  float aveETtime = sumETtime / static_cast<float>(n_evals);

  // print output timing results
  // note one could use PROFILELOG for these lines
  std::cout << "T: MSK generation time:        "
            << "\t" << MSKeyGenTime << " ms" << std::endl;
  std::cout << "T: Obfuscation execution time: "
            << "\t" << conKeyGenTime << " ms" << std::endl;
  std::cout << "T: Average token evaluation time:  "
            << "\t" << aveTokenTime / 1000 << " ms" << std::endl;
  std::cout << "T: Average evaluation time:  "
            << "\t" << aveTime / 1000 << " ms" << std::endl;
  std::cout << "T: Average EqualTest time:  "
            << "\t" << aveETtime / 1000 << " ms" << std::endl;
  std::cout << "T: Average total evaluation time:       "
            << "\t" << aveTokenTime / 1000 + aveTime / 1000 + aveETtime / 1000
            << " ms" << std::endl;

  return errorflag;
}

string RandomBooleanString(usint length) {
  static std::default_random_engine e{};
  static std::uniform_int_distribution<int> d{0, 1};

  string str("");
  for (usint i = 0; i < length; i++) {
    str += std::to_string(d(e));
  }

  return str;
}
