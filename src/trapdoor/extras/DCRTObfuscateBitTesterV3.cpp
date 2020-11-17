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

#include <fstream>
#include <iostream>
#include "obfuscation/lweconjunctionobfuscate.h"
#include "utils/debug.h"
#include "utils/parallel.h"

using namespace lbcrypto;

bool CONJOBF(bool dbg_flag, int n_evals, int n);  // defined later

// main()   need this for Kurts makefile to ignore this.
int main(int argc, char *argv[]) {
  bool errorflag = false;

  if (argc < 2) {  // called with no arguments
    std::cout << "Usage is `" << argv[0]
              << " arg1 arg2 arg3' where: " << std::endl;
    std::cout << "  arg1 indicate verbosity of output. Possible values are 0 "
                 "or 1 with 1 being verbose.  Default is 0."
              << std::endl;
    std::cout << "  arg2 indicates number of evaluation operations to run.  "
                 "Possible values are 1, 2 or 3.  Default is 1."
              << std::endl;
    std::cout << "If no input is given, then this message is displayed, "
                 "defaults are assumed and user is prompted for ring dimension."
              << std::endl;
  }
  bool dbg_flag = false;

  if (argc >= 2) {
    if (atoi(argv[1]) != 0) {
#if !defined(NDEBUG)
      dbg_flag = true;
      // std::cout << "setting dbg_flag true" << std::endl;
#endif
    }
  }

  int n_evals = 1;

  if (argc >= 3) {
    if (atoi(argv[2]) < 0) {
      n_evals = 1;
    } else if (atoi(argv[2]) >= 3) {
      n_evals = 3;
    } else {
      n_evals = atoi(argv[2]);
    }
  }

  std::cerr << "Running " << argv[0] << " with "
            << ParallelControls::GetNumProcs() << " processors and " << n_evals
            << " evaluation[s]." << std::endl;

  for (auto i = 0; i < 3; i++) {
    std::cout << "Run " << i << std::endl;
    for (usint n = 1024; n < 2048; n = 2 * n) {
      errorflag = CONJOBF(dbg_flag, n_evals, n);
      if (errorflag) return static_cast<int>(errorflag);
    }
    xalloc_stats();
  }
  // system("PAUSE");

  return static_cast<int>(errorflag);
}

//////////////////////////////////////////////////////////////////////
bool CONJOBF(bool dbg_flag, int n_evals, int n) {
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
  // 54 bits
  // BigInteger modulus("9007199254741169");
  // BigInteger rootOfUnity("7629104920968175");

  usint chunkSize = 1;
  usint base = 1 << 20;

  // Generate the test pattern
  std::string inputPattern = "1?10";
  ClearLWEConjunctionPattern<DCRTPoly> clearPattern(inputPattern);

  ObfuscatedLWEConjunctionPattern<DCRTPoly> obfuscatedPattern;
  obfuscatedPattern.SetChunkSize(chunkSize);
  obfuscatedPattern.SetBase(base);
  obfuscatedPattern.SetLength(clearPattern.GetLength());
  obfuscatedPattern.SetRootHermiteFactor(1.006);

  LWEConjunctionObfuscationAlgorithm<DCRTPoly> algorithm;

  // Variables for timing
  double timeDGGSetup(0.0), timeKeyGen(0.0), timeObf(0.0), timeEval1(0.0),
      timeEval2(0.0), timeEval3(0.0), timeTotal(0.0);

  double stdDev = SIGMA;
  typename DCRTPoly::DggType dgg =
      DCRTPoly::DggType(stdDev);  // Create the noise generator

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
  PROFILELOG("base=" << base);

  typename DCRTPoly::DugType dug;
  typename DCRTPoly::TugType tug;

  PROFILELOG("\nCryptosystem initialization: Performing precomputations...");

  // This code is run only when performing execution time measurements

  // Precomputations for FTT
  DiscreteFourierTransform::PreComputeTable(m);

  ////////////////////////////////////////////////////////////
  // Test the cleartext pattern
  ////////////////////////////////////////////////////////////

  DEBUG(" \nCleartext pattern: ");
  DEBUG(clearPattern.GetPatternString());

  DEBUG(" \nCleartext pattern length: ");
  DEBUG(clearPattern.GetLength());

  std::string inputStr1 = "1010";
  bool out1 = algorithm.Evaluate(clearPattern, inputStr1);
  DEBUG(" \nCleartext pattern evaluation of: " << inputStr1 << " is " << out1);

  std::string inputStr2 = "1111";
  bool out2 = algorithm.Evaluate(clearPattern, inputStr2);
  DEBUG(" \nCleartext pattern evaluation of: " << inputStr2 << " is " << out2);

  std::string inputStr3 = "1110";
  bool out3 = algorithm.Evaluate(clearPattern, inputStr3);
  DEBUG(" \nCleartext pattern evaluation of: " << inputStr3 << " is " << out3);

  ////////////////////////////////////////////////////////////
  // Generate and test the obfuscated pattern
  ////////////////////////////////////////////////////////////

  bool result1 = false;
  bool result2 = false;
  bool result3 = false;

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

  PROFILELOG("Evaluation 1 started");
  TIC(t1);
  result1 = algorithm.Evaluate(obfuscatedPattern, inputStr1);
  timeEval1 = TOC(t1);
  DEBUG(" \nCleartext pattern evaluation of: " << inputStr1 << " is " << result1
                                               << ".");
  PROFILELOG("Evaluation 1 execution time: "
             << "\t" << timeEval1 << " ms");

  bool errorflag = false;
  if (result1 != out1) {
    std::cout << "ERROR EVALUATING 1" << std::endl;
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
      std::cout << "ERROR EVALUATING 2" << std::endl;
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
      std::cout << "ERROR EVALUATING 3" << std::endl;
      errorflag |= true;
    }
  }

  // get the total program run time.
  timeTotal = TOC(t_total);

  // print output timing results
  // note one could use PROFILELOG for these lines
  std::cout << "Timing Summary for n = " << m / 2 << std::endl;
  std::cout << "T: DGG setup time:        "
            << "\t" << timeDGGSetup << " ms" << std::endl;
  std::cout << "T: Key generation time:        "
            << "\t" << timeKeyGen << " ms" << std::endl;
  std::cout << "T: Obfuscation execution time: "
            << "\t" << timeObf << " ms" << std::endl;
  std::cout << "T: Eval 1 execution time:  "
            << "\t" << timeEval1 << " ms" << std::endl;
  std::cout << "T: Eval 2 execution time:  "
            << "\t" << timeEval2 << " ms" << std::endl;
  std::cout << "T: Eval 3 execution time:  "
            << "\t" << timeEval3 << " ms" << std::endl;
  std::cout << "T: Total execution time:       "
            << "\t" << timeTotal << " ms" << std::endl;

  if (errorflag) {
    std::cout << "FAIL " << std::endl;
  } else {
    std::cout << "SUCCESS " << std::endl;
  }

  DiscreteFourierTransform::Reset();

  return (errorflag);
}
