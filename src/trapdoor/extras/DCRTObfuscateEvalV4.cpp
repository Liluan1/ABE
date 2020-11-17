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
#include "utils/serialize-binary.h"

using namespace lbcrypto;

string RandomBooleanString(usint length);

bool EvaluateConjObfs(bool dbg_flag, int n, usint pattern_size, usint num_evals,
                      bool rand_evals, bool single_flag);  // defined later

// main()   need this for Kurts makefile to ignore this.
int main(int argc, char* argv[]) {
  int opt;                // used in getting options
  bool dbg_flag = false;  // if true print debugging info
  usint pattern_size(8);  // size of the cleartext pattern
  int n_bits = 8;         // number of bits in underlying vector length
  usint n_evals = 3;      // number of evaluations to run
  bool rand_evals = false;
  bool single_flag = false;
  while ((opt = getopt(argc, argv, "de:k:b:hs")) != -1) {
    switch (opt) {
      case 'd':
        dbg_flag = true;
        std::cout << "setting dbg_flag true" << std::endl;
        break;
      case 'e':
        n_evals = atoi(optarg);
        if (n_evals > 3) {
          rand_evals = true;
        }
        break;
      case 'k':
        n_bits = atoi(optarg);
        if (n_bits < 8) {
          n_bits = 8;
          std::cout << "setting n_bits to minimum size of 8" << std::endl;
        } else if (n_bits >= 13) {
          n_bits = 13;
          std::cout << "setting n_bits to maximum size of 13" << std::endl;
        }
        break;
      case 'b':
        pattern_size = atoi(optarg);
        if ((pattern_size != 8) && (pattern_size != 16) &&
            (pattern_size != 32) && (pattern_size != 40) &&
            (pattern_size != 64)) {
          std::cout << "bad pattern size: " << pattern_size << std::endl;
          exit(EXIT_FAILURE);
        }
        break;
      case 's':
        single_flag = true;
        std::cout << "reading from single file" << std::endl;
        break;
      case 'h':
      default: /* '?' */
        std::cerr
            << "Usage: " << argv[0] << " <arguments> " << std::endl
            << "required arguments:" << std::endl
            << "  -k  bitsize of security parameter (ring dimension = 2^k) "
               "[8:13] (10)"
            << std::endl
            << "  -b  pattern_bitsize [8 16 32 40 64] (8)" << std::endl
            << "  -e num_evals (3) {If >3, then all evaluations will be random}"
            << std::endl
            << "optional arguments:" << std::endl
            << "  -d  (false) sets debug flag true " << std::endl
            << "  -s  (false) use single input file" << std::endl
            << "  -h  prints this message" << std::endl;
        exit(EXIT_FAILURE);
    }
  }

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

  unsigned int n = 1 << n_bits;
  bool errorflag = EvaluateConjObfs(dbg_flag, n, pattern_size, n_evals,
                                    rand_evals, single_flag);
  return (static_cast<int>(errorflag));
}

//////////////////////////////////////////////////////////////////////
bool EvaluateConjObfs(bool dbg_flag, int n, usint pattern_size, usint n_evals,
                      bool rand_evals, bool single_flag) {
  // if dbg_flag == true; print debug outputs
  // n = size of vectors to use

  // returns
  //  errorflag = 1 if fail

  TimeVar t1, t_total;  // for TIC TOC
  float timeRead(0.0);

  usint m = 2 * n;

  // Read the test pattern from the file
  ClearLWEConjunctionPattern<DCRTPoly> clearPattern("");
  string clearFileName =
      "cp" + std::to_string(n) + "_" + std::to_string(pattern_size) + ".serial";

  DEBUG("reading clearPattern from file: " << clearFileName << ".serial");
  TIC(t1);
  Serial::DeserializeFromFile(clearFileName, clearPattern, SerType::BINARY);
  timeRead = TOC(t1);
  PROFILELOG("Read time: "
             << "\t" << timeRead << " ms");

  string obfFileName =
      "op" + std::to_string(n) + "_" + std::to_string(pattern_size) + ".serial";
  // note this is for debug -- will move to evaluate program once it all works
  ObfuscatedLWEConjunctionPattern<DCRTPoly> obfuscatedPattern;

  if (single_flag) {
    std::cout << "Deserializing Obfuscated Pattern from file " << obfFileName
              << std::endl;
  } else {
    std::cout << "Deserializing Obfuscated Pattern from fileset " << obfFileName
              << std::endl;
  }
  TIC(t1);
  Serial::DeserializeFromFile(obfFileName, obfuscatedPattern, SerType::BINARY);
  timeRead = TOC(t1);
  PROFILELOG("Done, Read time: "
             << "\t" << timeRead << " ms");

  TIC(t_total);  // start timer for total time

  LWEConjunctionObfuscationAlgorithm<DCRTPoly> algorithm;

  // Variables for timing
  // todo make eval an array
  vector<double> timeEval(n_evals);

  double timeTotal(0.0);

  const shared_ptr<typename DCRTPoly::Params> ilParams =
      obfuscatedPattern.GetParameters();

  m = ilParams->GetCyclotomicOrder();

  double stdDev = SIGMA;
  DCRTPoly::DggType dgg(stdDev);  // Create the noise generator

  // Precomputations for FTT
  DiscreteFourierTransform::PreComputeTable(m);

  ////////////////////////////////////////////////////////////
  // Test the cleartext pattern
  ////////////////////////////////////////////////////////////

  DEBUG(" \nCleartext pattern: ");
  DEBUG(clearPattern.GetPatternString());

  DEBUG(" \nCleartext pattern length: ");
  DEBUG(clearPattern.GetLength());

  vector<string> inputStr(std::max<usint>(
      3, n_evals));  // this must be at least size three
                     // to support nonrandom histoircal test patterns

  if (!rand_evals) {  // use historic evaluation patterns
    switch (pattern_size) {
      case 8:
        inputStr[0] = "11100100";
        inputStr[1] = "11001101";
        inputStr[2] = "10101101";
        break;

      case 16:
        inputStr[0] = "1110010011100100";
        inputStr[1] = "1100110111001101";
        inputStr[2] = "1010110110101101";
        break;

      case 32:  // 32 bit test
        inputStr[0] = "11100100111001001110010011100100";
        inputStr[1] = "11001101110011011100110111001111";
        inputStr[2] = "10101101101011011010110110101101";
        break;

      case 40:  // 40 bit test
        inputStr[0] = "1110010011100100111001001110010011100100";
        inputStr[1] = "1100110111001101110011011100111111001101";
        inputStr[2] = "1010110110101101101011011010110110101101";
        break;

      case 64:  // 64 bit test
        inputStr[0] =
            "1110010011100100111001001110010011100100111001001110010011100100";
        inputStr[1] =
            "1100110111001101110011011100111111001101110011011100110111001111";
        inputStr[2] =
            "1010110110101101101011011010110110101101101011011010110110101101";
        break;

      default:
        std::cout << "bad input pattern length selected (must be 8, 16, 32, 40 "
                     "or 64). "
                  << std::endl;
        exit(-1);
    }

  } else {  // use random evaluation patterns
    for (usint i = 0; i < n_evals; i++) {
      inputStr[i] = RandomBooleanString(pattern_size);
    }
  }
  vector<bool> out(n_evals);
  vector<bool> result(n_evals);
  bool errorflag = false;

  ////////////////////////////////////////////////////////////
  // test the obfuscated pattern
  ////////////////////////////////////////////////////////////
  PROFILELOG("Evaluation started");
  for (usint i = 0; i < n_evals; i++) {
    out[i] = algorithm.Evaluate(clearPattern, inputStr[i]);
    DEBUG(" \nCleartext pattern evaluation of: " << inputStr[i] << " is "
                                                 << out[i]);

    TIC(t1);
    result[i] = algorithm.Evaluate(obfuscatedPattern, inputStr[i]);
    timeEval[i] = TOC(t1);

    DEBUG(" \nObfuscated pattern evaluation of: " << inputStr[i] << " is "
                                                  << result[i] << ".");
    PROFILELOG("Evaluation " << i << " execution time: "
                             << "\t" << timeEval[i] << " ms");

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
              << "\t" << timeEval[i] << " ms" << std::endl;
  }
  aveTime /= static_cast<float>(n_evals);

  std::cout << "T: Average evaluation execution time:  "
            << "\t" << aveTime << " ms" << std::endl;
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

string RandomBooleanString(usint length) {
  static std::default_random_engine e{};
  static std::uniform_int_distribution<int> d{0, 1};

  string str("");
  for (usint i = 0; i < length; i++) {
    str += std::to_string(d(e));
  }

  return str;
}
