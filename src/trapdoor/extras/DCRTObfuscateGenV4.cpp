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
#include <fstream>
#include <iostream>
#include "obfuscation/lweconjunctionobfuscate.h"
#include "utils/debug.h"
#include "utils/parallel.h"
#include "utils/serialize-binary.h"

using namespace lbcrypto;

// forward definitions to be defined later
bool GenerateConjObfs(bool dbg_flag, int n, usint pattern_size, bool eval_flag,
                      usint n_evals, bool rand_evals, bool pretty_flag,
                      bool verify_flag, bool single_flag, string inputFileName);

// main()   need this for Kurts makefile to ignore this.
int main(int argc, char *argv[]) {
  // todo: allow clear pattern as input

  int opt;                // used in getting options
  bool dbg_flag = false;  // if true print debugging info
  usint pattern_size(8);  // size of the cleartext pattern
  // if true, also run evaluation to verify correct generation.
  bool eval_flag = false;
  int n_bits = 8;     // number of bits in underlying vector length
  usint n_evals = 0;  // number of evaluations to run
  bool rand_evals = false;
  bool pretty_flag = false;
  bool verify_flag = false;
  bool single_flag = false;
  string inputFileName("");
  while ((opt = getopt(argc, argv, "dpve:i:k:b:sh")) != -1) {
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
      case 'e':
        n_evals = atoi(optarg);
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
        std::cerr
            << "Usage: " << argv[0] << " <arguments> " << std::endl
            << "required arguments:" << std::endl
            << "  -k  bitsize of security parameter (ring dimension = 2^k) "
               "[8:13] (10)"
            << std::endl
            << "  -b  pattern_bitsize [8 16 32 40 64] (8)" << std::endl
            << "optional arguments:" << std::endl
            << "  -d  (false) sets debug flag true " << std::endl
            << "  -e num_evals (0) {If >3, then all evaluations will be random}"
            << std::endl
            << "  -i input_file  reads  pattern from inut file (concatenated "
               "or zeropadded to -b}"
            << std::endl
            << "  -p  (false) enable prettyprint json output" << std::endl
            << "  -v  (false) enable verification of serialization" << std::endl
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

  // 32 bit test would run n_bits = 10 .. < 10+bit range no max, default 3
  // 48 bit test runs 10.. 12
  // 64 bit test ran from 1..13

  unsigned int n = 1 << n_bits;

  bool errorflag = GenerateConjObfs(dbg_flag, n, pattern_size, eval_flag,
                                    n_evals, rand_evals, pretty_flag,
                                    verify_flag, single_flag, inputFileName);

  return static_cast<int>(errorflag);
}

//////////////////////////////////////////////////////////////////////
bool GenerateConjObfs(bool dbg_flag, int n, usint pattern_size, bool eval_flag,
                      usint n_evals, bool rand_evals, bool pretty_flag,
                      bool verify_flag, bool single_flag,
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
  double timeDGGSetup(0.0), timeKeyGen(0.0), timeObf(0.0), timeSerial(0.0),
      timeTotal(0.0);

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
  std::cout << "T: DGG setup time:        "
            << "\t" << timeDGGSetup << " ms" << std::endl;
  std::cout << "T: Key generation time:        "
            << "\t" << timeKeyGen << " ms" << std::endl;
  std::cout << "T: Obfuscation execution time: "
            << "\t" << timeObf << " ms" << std::endl;
  std::cout << "T: Total execution time:       "
            << "\t" << timeTotal << " ms" << std::endl;
  std::cout << "T: Serialization execution time: "
            << "\t" << timeSerial << " ms" << std::endl;

  if (eval_flag) {
    ////////////////////////////////////////////////////////////
    // Test the cleartext and obfuscated pattern
    ////////////////////////////////////////////////////////////

    std::cout << "Random evaluation not implemented yet";

    if (n_evals > 3) {
      n_evals = 3;
      std::cout << " limiting num_evals to " << n_evals << std::endl;
    }
    DEBUG(" \nCleartext pattern: ");
    DEBUG(clearPattern.GetPatternString());

    DEBUG(" \nCleartext pattern length: ");
    DEBUG(clearPattern.GetLength());

    std::string inputStr1("");
    std::string inputStr2("");
    std::string inputStr3("");
    switch (pattern_size) {
      case 8:  // 8 bit test
        inputStr1 = "11100100";
        inputStr2 = "11001101";
        inputStr3 = "10101101";
        break;

      case 16:
        inputStr1 = "1110010011100100";
        inputStr2 = "1100110111001101";
        inputStr3 = "1010110110101101";
        break;

      case 32:  // 32 bit test
        inputStr1 = "11100100111001001110010011100100";
        inputStr2 = "11001101110011011100110111001111";
        inputStr3 = "10101101101011011010110110101101";
        break;

      case 40:  // 40 bit test
        inputStr1 = "1110010011100100111001001110010011100100";
        inputStr2 = "1100110111001101110011011100111111001101";
        inputStr3 = "1010110110101101101011011010110110101101";
        break;

      case 64:  // 64 bit test
        inputStr1 =
            "1110010011100100111001001110010011100100111001001110010011100100";
        inputStr2 =
            "1100110111001101110011011100111111001101110011011100110111001111";
        inputStr3 =
            "1010110110101101101011011010110110101101101011011010110110101101";
        break;

      default:
        std::cout
            << "bad input pattern length selected (must be 32, 40 or 64). "
            << std::endl;
        exit(-1);
    }

    // std::string inputStr1 = "11100100"; //8 bit test
    // std::string inputStr2 = "11001101"; //8 bit test
    // std::string inputStr3 = "10101101"; //8 bit test

    bool out1 = algorithm.Evaluate(clearPattern, inputStr1);
    DEBUG(" \nCleartext pattern evaluation of: " << inputStr1 << " is "
                                                 << out1);

    bool out2 = algorithm.Evaluate(clearPattern, inputStr2);
    DEBUG(" \nCleartext pattern evaluation of: " << inputStr2 << " is "
                                                 << out2);

    bool out3 = algorithm.Evaluate(clearPattern, inputStr3);
    DEBUG(" \nCleartext pattern evaluation of: " << inputStr3 << " is "
                                                 << out3);

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
    DEBUG(" \nCleartext pattern evaluation of: " << inputStr1 << " is "
                                                 << result1 << ".");
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

    // get the total program run time.
    timeTotal = TOC(t_total);

    // print output timing results
    // note one could use PROFILELOG for these lines
    std::cout << "Timing Summary for n = " << m / 2 << std::endl;
    std::cout << "T: Eval 1 execution time:  "
              << "\t" << timeEval1 << " ms" << std::endl;
    std::cout << "T: Eval 2 execution time:  "
              << "\t" << timeEval2 << " ms" << std::endl;
    std::cout << "T: Eval 3 execution time:  "
              << "\t" << timeEval3 << " ms" << std::endl;
    std::cout << "T: Average evaluation execution time:  "
              << "\t" << (timeEval1 + timeEval2 + timeEval3) / 3 << " ms"
              << std::endl;
    std::cout << "T: Total execution time:       "
              << "\t" << timeTotal << " ms" << std::endl;

    if (errorflag) {
      std::cout << "FAIL " << std::endl;
    } else {
      std::cout << "SUCCESS " << std::endl;
    }
  }  // if eval_flag

  DiscreteFourierTransform::Reset();

  return (errorflag);
}
