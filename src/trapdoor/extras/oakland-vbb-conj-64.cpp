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
// Note must must be before all headers

#include <getopt.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <random>
#include "obfuscation/lweconjunctionobfuscate.h"
#include "utils/debug.h"
#include "utils/parallel.h"

using namespace lbcrypto;

bool CONJOBF(bool dbg_flag, size_t n_evals, int n);  // defined later

string RandomBooleanString(usint length);

// main()   need this for Kurts makefile to ignore this.
int main(int argc, char *argv[]) {
  PalisadeParallelControls.Enable();

  int opt;                 // used in getting options
  bool dbg_flag = false;   // if true print debugging info
  usint pattern_size(64);  // size of the cleartext pattern
  int n_bits = 8;          // number of bits in underlying vector length
  int n_evals = 100;       // number of evaluations to run

  while ((opt = getopt(argc, argv, "e:k:d")) != -1) {
    switch (opt) {
      case 'd':
        dbg_flag = true;
        std::cout << "setting dbg_flag true" << std::endl;
        break;
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
        std::cerr
            << "Usage: " << argv[0] << " <arguments> " << std::endl
            << "arguments:" << std::endl
            << "  -e  number of evaluations (100) {evaluations will be random}"
            << std::endl
            << "  -k  bitsize of security parameter (ring dimension = 2^k) "
               "[8:13] (10)"
            << std::endl
            << "  -d  (false) sets debug flag true " << std::endl
            << "  -h  (false) prints this message" << std::endl;
        exit(EXIT_FAILURE);
    }
  }

  std::cerr << "\n\nStarting the demo\n" << std::endl;

  DEBUG("DEBUG IS TRUE");
  PROFILELOG("PROFILELOG IS TRUE");
#ifdef PROFILE
  std::cout << "PROFILE is defined" << std::endl;
#endif
#ifdef NDEBUG
  std::cout << "NDEBUG is defined" << std::endl;
#endif

  std::cerr << "Running " << argv[0] << " with security parameter "
            << std::to_string(1 << n_bits) << ". Pattern length "
            << pattern_size << "." << std::endl;

  std::cerr << "Running " << argv[0] << " with "
            << ParallelControls::GetNumProcs() << " processors and "
            << ParallelControls().GetNumThreads() << " threads. " << std::endl;

  // 32 bit test would run n_bits = 10 .. < 10+bit range no max, default 3
  // 48 bit test runs 10.. 12
  // 64 bit test ran from 1..13

  unsigned int n = 1 << n_bits;

  bool errorflag = CONJOBF(dbg_flag, n_evals, n);

  return static_cast<int>(errorflag);
}

//////////////////////////////////////////////////////////////////////
bool CONJOBF(bool dbg_flag, size_t n_evals, int n) {
  // if dbg_flag == true; print debug outputs
  // n_evals = 1,2,3 number of evaluations to perform
  // returns
  //  errorflag = # of bad evaluations

  DEBUG("DEBUG IS TRUE");
  PROFILELOG("PROFILELOG IS TRUE");
#ifdef PROFILE
  std::cout << "PROFILE is defined" << std::endl;
#endif
#ifdef NDEBUG
  std::cout << "NDEBUG is defined" << std::endl;
#endif

  TimeVar t1, t_total;  // for TIC TOC
  TIC(t_total);         // start timer for total time

  usint m = 2 * n;

  usint chunkSize = 8;
  usint base = 1 << 20;

  // if (n > 1<<11)
  //  base = 1<<18;

  // Generate the test pattern
  std::string inputPattern =
      "1?10?10?1?10?10?1?10?10?1?10??0?1?10?10?1?10?10?1?10?10?1?10??0?";
  ClearLWEConjunctionPattern<DCRTPoly> clearPattern(inputPattern);

  ObfuscatedLWEConjunctionPattern<DCRTPoly> obfuscatedPattern;
  obfuscatedPattern.SetChunkSize(chunkSize);
  obfuscatedPattern.SetBase(base);
  obfuscatedPattern.SetLength(clearPattern.GetLength());
  obfuscatedPattern.SetRootHermiteFactor(1.006);

  LWEConjunctionObfuscationAlgorithm<DCRTPoly> algorithm;

  // Variables for timing
  double timeKeyGen(0.0), timeObf(0.0);

  double stdDev = SIGMA;
  DCRTPoly::DggType dgg(stdDev);  // Create the noise generator

  // Finds q using the correctness constraint for the given value of n
  algorithm.ParamsGen(dgg, &obfuscatedPattern, m / 2);

  // this code finds the values of q and n corresponding to the root Hermite
  // factor in obfuscatedPattern algorithm.ParamsGen(dgg, &obfuscatedPattern);

  const shared_ptr<typename DCRTPoly::Params> ilParams =
      obfuscatedPattern.GetParameters();

  const BigInteger &modulus = ilParams->GetModulus();
  const BigInteger &rootOfUnity = ilParams->GetRootOfUnity();
  m = ilParams->GetCyclotomicOrder();

  PROFILELOG("\nq = " << modulus);
  PROFILELOG("rootOfUnity = " << rootOfUnity);
  PROFILELOG("n = " << m / 2);
  PROFILELOG(printf("delta=%lf", obfuscatedPattern.GetRootHermiteFactor()));
  PROFILELOG("\nbase = " << base);

  typename DCRTPoly::DugType dug;
  typename DCRTPoly::TugType tug;

  PROFILELOG("\nCryptosystem initialization: Performing precomputations...");

  // This code is run only when performing execution time measurements

  // Precomputations for FTT
  DiscreteFourierTransform::PreComputeTable(m);

  std::cout << " \nCleartext pattern: " << std::endl;
  std::cout << clearPattern.GetPatternString() << std::endl;

  PROFILELOG("Key generation started");
  TIC(t1);
  algorithm.KeyGen(dgg, &obfuscatedPattern);
  timeKeyGen = TOC(t1);
  PROFILELOG("Key generation time: "
             << "\t" << timeKeyGen << " ms");

  BinaryUniformGenerator dbg = BinaryUniformGenerator();

  DEBUG("Obfuscation Execution started");
  TIC(t1);
  algorithm.Obfuscate(clearPattern, dgg, tug, &obfuscatedPattern);
  timeObf = TOC(t1);
  PROFILELOG("Obfuscation time: "
             << "\t" << timeObf << " ms");

  vector<string> inputStr(n_evals);

  for (usint i = 0; i < n_evals; i++) {
    inputStr[i] = RandomBooleanString(clearPattern.GetLength());
  }

  vector<bool> out(n_evals);
  vector<bool> result(n_evals);
  bool errorflag = false;

  // Variables for timing
  vector<double> timeEval(n_evals);

  // 3 known patterns

  double timeTotal(0.0);

  std::string inputStr1 =
      "1110010011100100111001001110010011100100111001001110010011100100";
  std::string inputStr2 =
      "1100110111001101110011011100111111001101110011011100110111001111";
  std::string inputStr3 =
      "1010110110101101101011011010110110101101101011011010110110101101";

  bool out1 = algorithm.Evaluate(clearPattern, inputStr1);
  DEBUG(" \nCleartext pattern evaluation of: " << inputStr1 << " is " << out1);

  bool out2 = algorithm.Evaluate(clearPattern, inputStr2);
  DEBUG(" \nCleartext pattern evaluation of: " << inputStr2 << " is " << out2);

  bool out3 = algorithm.Evaluate(clearPattern, inputStr3);
  DEBUG(" \nCleartext pattern evaluation of: " << inputStr3 << " is " << out3);

  ////////////////////////////////////////////////////////////
  // Generate and test the obfuscated pattern
  ////////////////////////////////////////////////////////////
  double timeEval1(0.0), timeEval2(0.0), timeEval3(0.0);

  // todo make this a loop
  bool result1 = false;
  bool result2 = false;
  bool result3 = false;
  std::cout << " \nCleartext pattern: " << std::endl;
  std::cout << clearPattern.GetPatternString() << std::endl;

  PROFILELOG("Evaluation started");
  DEBUG("====== just before eval ");
  DEBUGEXP(*(obfuscatedPattern.GetParameters()));

  DEBUG("====== ");
  TIC(t1);
  result1 = algorithm.Evaluate(obfuscatedPattern, inputStr1);
  timeEval1 = TOC(t1);
  DEBUG(" \nCleartext pattern evaluation of: " << inputStr1 << " is " << result1
                                               << ".");
  PROFILELOG("Evaluation 1 execution time: "
             << "\t" << timeEval1 << " ms");

  errorflag = false;
  if (result1 != out1) {
    std::cout << "ERROR EVALUATING 1 "
              << " got " << result1 << " wanted " << out1 << std::endl;
    errorflag |= true;
  }
  if (n_evals > 1) {
    PROFILELOG("Evaluation 2 started");
    TIC(t1);
    result2 = algorithm.Evaluate(obfuscatedPattern, inputStr2);
    timeEval2 = TOC(t1);
    DEBUG(" \nCleartext pattern evaluation of: " << inputStr2 << " is "
                                                 << result2 << ".");
    PROFILELOG("Evaluation 2 execution time: "
               << "\t" << timeEval2 << " ms");

    if (result2 != out2) {
      std::cout << "ERROR EVALUATING 2"
                << " got " << result2 << " wanted " << out2 << std::endl;
      errorflag |= true;
    }
  }

  if (n_evals > 2) {
    PROFILELOG("Evaluation 3 started");
    TIC(t1);
    result3 = algorithm.Evaluate(obfuscatedPattern, inputStr3);
    timeEval3 = TOC(t1);
    DEBUG("\nCleartext pattern evaluation of: " << inputStr3 << " is "
                                                << result3 << ".");
    PROFILELOG("Evaluation 3 execution time: "
               << "\t" << timeEval3 << " ms");
    if (result3 != out3) {
      std::cout << "ERROR EVALUATING 3"
                << " got " << result3 << " wanted " << out3 << std::endl;
      errorflag |= true;
    }
  }

  size_t counter = 0;

  ////////////////////////////////////////////////////////////
  // main performance evaluation based on n_eval random inputs
  ////////////////////////////////////////////////////////////
  PROFILELOG("\nEvaluation started");
  for (usint i = 0; i < n_evals; i++) {
    out[i] = algorithm.Evaluate(clearPattern, inputStr[i]);
    DEBUG(" \nCleartext pattern evaluation of: " << inputStr[i] << " is "
                                                 << out[i]);

    if (out[i]) counter++;

    TIC(t1);
    result[i] = algorithm.Evaluate(obfuscatedPattern, inputStr[i]);
    timeEval[i] = TOC_US(t1);

    DEBUG(" \nObfuscated pattern evaluation of: " << inputStr[i] << " is "
                                                  << result[i] << ".");
    // PROFILELOG("Evaluation "<<i<<" execution time: " << "\t" << timeEval[i]
    // << " ms");

    if (result[i] != out[i]) {
      std::cout << "ERROR EVALUATING " << i << " got " << result[i]
                << " wanted " << out[i] << std::endl;
      errorflag |= true;
    }
  }  // end eval loop
  // get the total program run time.
  timeTotal = TOC(t_total);

  // print output timing results
  // note one could use PROFILELOG for these lines
  std::cout << "Timing Summary for n = " << m / 2 << std::endl;
  float aveTime = 0.0;
  for (usint i = 0; i < n_evals; i++) {
    aveTime += timeEval[i];
    std::cout << "T: Eval " << i << " execution time:  "
              << "\t" << timeEval[i] / 1000 << " ms" << std::endl;
  }
  aveTime /= static_cast<float>(n_evals);

  if (errorflag) {
    std::cout << "FAIL " << std::endl;
  } else {
    std::cout << "SUCCESS " << std::endl;
  }

  // print output timing results
  // note one could use PROFILELOG for these lines
  std::cout << "\nTiming Summary for n = " << m / 2 << std::endl;
  std::cout << "T: Key generation time:        "
            << "\t" << timeKeyGen << " ms" << std::endl;
  std::cout << "T: Obfuscation execution time: "
            << "\t" << timeObf << " ms" << std::endl;
  std::cout << "T: Average evaluation execution time:  "
            << "\t" << aveTime / 1000 << " ms" << std::endl;
  std::cout << "T: Total execution time:       "
            << "\t" << timeTotal << " ms" << std::endl;

  DiscreteFourierTransform::Reset();

  return (errorflag);
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
