// @file palisade-vbb-conj.cpp - Example of distributed Virtual Black Box
// Obfuscation for conjunctions with wild cards
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

// forward definitions to be defined later
bool GenerateConjObfs(bool dbg_flag, int n, usint pattern_size,
                      bool pretty_flag, bool verify_flag, bool single_flag,
                      string inputFileName);

bool EvaluateConjObfs(bool dbg_flag, int n, usint pattern_size, usint num_evals,
                      bool rand_evals, bool single_flag);  // defined later

string RandomBooleanString(usint length);

// main()   need this for Kurts makefile to ignore this.
int main(int argc, char *argv[]) {
  // todo: allow clear pattern as input

  int opt;                // used in getting options
  bool dbg_flag = false;  // if true print debugging info
  usint pattern_size(8);  // size of the cleartext pattern
  bool eval_flag =
      false;       // if true, also run evaluation to verify correct generation.
  int n_bits = 8;  // number of bits in underlying vector length
  int n_evals = 0;  // number of evaluations to run
  bool rand_evals = false;
  bool pretty_flag = false;
  bool verify_flag = false;
  bool single_flag = false;
  bool gen_flag = false;
  string inputFileName("");
  while ((opt = getopt(argc, argv, "d:p:v:g:e:i:k:b:s:h")) != -1) {
    switch (opt) {
      case 'd':
        dbg_flag = true;
        std::cout << "setting dbg_flag true" << std::endl;
        break;
      case 'p':
        pretty_flag = true;
        std::cout << "setting prettyprint_flag true" << std::endl;
        break;
      case 'v':
        verify_flag = true;
        std::cout << "setting verify_flag true" << std::endl;
        break;
      case 'g':
        gen_flag = true;
        std::cout << "setting gen_flag true" << std::endl;
        break;
      case 'e':
        n_evals = atoi(optarg);
        if (n_evals < 0) n_evals = 0;
        if (n_evals > 3) {
          rand_evals = true;
        }
        eval_flag = true;
        break;
      case 'i':
        inputFileName = optarg;
        std::cout << "reading input from file " << inputFileName << std::endl;
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
        std::cout << "writing to single file" << std::endl;
        break;
      case 'h':
      default: /* '?' */
        std::cerr << "Usage: " << argv[0] << " <arguments> " << std::endl
                  << "arguments:" << std::endl
                  << "  -g  (false) generate the keys and obfuscated program "
                  << std::endl
                  << "  -e  number of evaluations (0) {If >3, then all "
                     "evaluations will be random}"
                  << std::endl
                  << "  -k  bitsize of security parameter (ring dimension = "
                     "2^k) [8:13] (10)"
                  << std::endl
                  << "  -b  pattern_bitsize [8 16 32 40 64] (8)" << std::endl
                  << "  -i  input_file  reads  pattern from input file "
                     "(concatenated or zeropadded to -b}"
                  << std::endl
                  << "  -d  (false) sets debug flag true " << std::endl
                  << "  -p  (false) enable prettyprint json output" << std::endl
                  << "  -v  (false) enable verification of serialization"
                  << std::endl
                  << "  -h  (false) prints this message" << std::endl;
        exit(EXIT_FAILURE);
    }
  }

  if (!((gen_flag) || (eval_flag))) {
    std::cerr
        << "At least one of -g (obfuscation generation) or -e (evaluation) "
           "flags has to be set. Run the command with the -h flag for help."
        << std::endl;
    exit(EXIT_FAILURE);
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

  // 32 bit test would run n_bits = 10 .. < 10+bit range no max, default 3
  // 48 bit test runs 10.. 12
  // 64 bit test ran from 1..13

  bool errorgenflag = false;
  bool errorevalflag = false;
  unsigned int n = 1 << n_bits;

  if (gen_flag) {
    errorgenflag = GenerateConjObfs(dbg_flag, n, pattern_size, pretty_flag,
                                    verify_flag, single_flag, inputFileName);
  }
  if (eval_flag) {
    errorevalflag = EvaluateConjObfs(dbg_flag, n, pattern_size, n_evals,
                                     rand_evals, single_flag);
  }

  return static_cast<int>(errorgenflag || errorevalflag);
}

//////////////////////////////////////////////////////////////////////
bool GenerateConjObfs(bool dbg_flag, int n, usint pattern_size,
                      bool pretty_flag, bool verify_flag, bool single_flag,
                      string inputFileName) {
  // if dbg_flag == true; print debug outputs
  // n = size of vectors to use (power of 2)
  // pattern_size = size of patterns (8, 32, 40, 64)
  // n_evals number of evals to run if 0..3 use preloaded, if > 3 use random.

  // returns
  //  errorflag = 1 if fail

  TimeVar t1, t_total;  // for TIC TOC
  TIC(t_total);         // start timer for total time

  usint m = 2 * n;

  usint chunkSize = 8;
  usint base = 1 << 20;

  // set inputPattern and adjust base for input pattern size
  std::string inputPattern("");
  if (inputFileName == "") {
    // set to historical input pattern
    switch (pattern_size) {
      case 8:
        inputPattern = "1?10?10?";  // 8 bit test
        break;

      case 16:
        inputPattern = "1?10?10?1?10?10?";  // 16 bit test
        break;

      case 32:
        inputPattern = "1?10?10?1?10?10?1?10?10?1?10??0?";  // 32 bit test
        break;

      case 40:
        inputPattern =
            "1?10?10?1?10?10?1?10?10?1?10??0?1?10?10?";  // 40 bit test
        break;

      case 64:
        // 64 bit test
        inputPattern =
            "1?10?10?1?10?10?1?10?10?1?10??0?1?10?10?1?10?10?1?10?10?1?10??0?";
        break;
      default:
        std::cout << "bad input pattern length selected (must be 8, 16, 32, 40 "
                     "or 64). "
                  << std::endl;
        exit(-1);
    }

  } else {
    // read the input pattern from a file

    std::ifstream inFile;
    inFile.open(inputFileName);
    if (inFile.is_open()) getline(inFile, inputPattern);

    inFile.close();
    std::cout << "Read input pattern: " << inputPattern << std::endl;

    // and adjust to pattern_size
    if (inputPattern.length() > pattern_size) {
      inputPattern =
          inputPattern.substr(0, pattern_size);  // lop off last characters
      std::cout << "Shortening input pattern to " << pattern_size
                << " bits: " << inputPattern << std::endl;
    }
    size_t l = inputPattern.length();
    if (l < pattern_size) {
      size_t n_needed = pattern_size - l;
      inputPattern.append(n_needed, '0');  // zero fill lop off last characters
      std::cout << "zero filling input pattern to " << pattern_size
                << " bits: " << inputPattern << std::endl;
    }
  }
  // now adjust base based in input size.
  switch (pattern_size) {
    case 8:
      if (n > 1 << 10)  // adjust for 8 bit test (use 32 bit test values)
        base = 1 << 15;
      break;

    case 16:
      if (n > 1 << 10)  // adjust for 16 bit test( use 32 bit test values)
        base = 1 << 15;
      break;

    case 32:
      if (n > 1 << 10)  // adjust for 32 bit test
        base = 1 << 15;
      break;

    case 40:
      break;

    case 64:
      if (n > 1 << 11)  // adjust for 64 bit test
        base = 1 << 18;
      break;
    default:
      std::cout
          << "bad input pattern length selected (must be 8, 16, 32, 40 or 64). "
          << std::endl;
      exit(-1);
  }
  DEBUG("clear pattern " << inputPattern);
  ClearLWEConjunctionPattern<DCRTPoly> clearPattern(inputPattern);

  string clearFileName =
      "cp" + std::to_string(n) + "_" + std::to_string(pattern_size) + ".serial";
  Serial::SerializeToFile(clearFileName, clearPattern, SerType::BINARY);

  if (verify_flag) {  // verify the serialization
    ClearLWEConjunctionPattern<DCRTPoly> testClearPattern("");

    Serial::DeserializeFromFile(clearFileName, testClearPattern,
                                SerType::BINARY);

    if (clearPattern.GetPatternString() ==
        testClearPattern.GetPatternString()) {
      std::cout << "Clear Pattern Serialization succeed" << std::endl;
    } else {
      std::cout << "Clear Pattern Serialization FAILED" << std::endl;
      std::cout << "    clear pattern:           "
                << clearPattern.GetPatternString() << std::endl;
      std::cout << "    recovered clear pattern: "
                << testClearPattern.GetPatternString() << std::endl;
    }
  }

  ObfuscatedLWEConjunctionPattern<DCRTPoly> obfuscatedPattern;
  obfuscatedPattern.SetChunkSize(chunkSize);
  obfuscatedPattern.SetBase(base);
  obfuscatedPattern.SetLength(clearPattern.GetLength());
  obfuscatedPattern.SetRootHermiteFactor(1.006);

  LWEConjunctionObfuscationAlgorithm<DCRTPoly> algorithm;

  // Variables for timing
  double timeKeyGen(0.0), timeObf(0.0), timeSerial(0.0), timeTotal(0.0);

  double stdDev = SIGMA;
  DCRTPoly::DggType dgg(stdDev);  // Create the noise generator

  // Finds q using the correctness constraint for the given value of n
  algorithm.ParamsGen(dgg, &obfuscatedPattern, m / 2);

  // this code finds the values of q and n corresponding to the root
  // Hermite factor in obfuscatedPattern

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

  ////////////////////////////////////////////////////////////
  // Generate and save the obfuscated pattern
  ////////////////////////////////////////////////////////////
  bool errorflag = false;

  PROFILELOG("Key generation started");
  TIC(t1);
  algorithm.KeyGen(dgg, &obfuscatedPattern);
  timeKeyGen = TOC(t1);
  PROFILELOG("Key generation time: "
             << "\t" << timeKeyGen << " ms");

  BinaryUniformGenerator dbg = BinaryUniformGenerator();

  DEBUG("Obfuscation Generation started");
  TIC(t1);
  algorithm.Obfuscate(clearPattern, dgg, tug, &obfuscatedPattern);
  timeObf = TOC(t1);
  PROFILELOG("Obfuscation time: "
             << "\t" << timeObf << " ms");
  // get the total program run time.
  timeTotal = TOC(t_total);

  DEBUG("Serializing Obfuscation");
  string obfFileName =
      "op" + std::to_string(n) + "_" + std::to_string(pattern_size) + ".serial";
  TIC(t1);
  Serial::SerializeToFile(obfFileName, obfuscatedPattern, SerType::BINARY);
  timeSerial = TOC(t1);
  PROFILELOG("Serialization  time: "
             << "\t" << timeSerial << " ms");

  if (verify_flag) {  // verify the serialization
    std::cout << "Verifying Serialization" << std::endl;
    ObfuscatedLWEConjunctionPattern<DCRTPoly> testObfuscatedPattern;

    Serial::DeserializeFromFile(obfFileName, testObfuscatedPattern,
                                SerType::BINARY);

    if (!obfuscatedPattern.Compare(testObfuscatedPattern)) {
      std::cout << "Serialization did verify" << std::endl;
    } else {
      std::cout << "Serialization verified" << std::endl;
    }

    DEBUG("Done");
  }

  // print output timing results
  // note one could use PROFILELOG for these lines
  std::cout << "Timing Summary for n = " << m / 2 << std::endl;
  std::cout << "T: Key generation time:        "
            << "\t" << timeKeyGen << " ms" << std::endl;
  std::cout << "T: Obfuscation execution time: "
            << "\t" << timeObf << " ms" << std::endl;
  std::cout << "T: Total execution time:       "
            << "\t" << timeTotal << " ms" << std::endl;
  std::cout << "T: Serialization execution time: "
            << "\t" << timeSerial << " ms" << std::endl;

  DiscreteFourierTransform::Reset();

  return (errorflag);
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
      "op" + std::to_string(n) + "_" + std::to_string(pattern_size);
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

  std::cout << "T: Average evaluation execution time:  "
            << "\t" << aveTime / 1000 << " ms" << std::endl;
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
