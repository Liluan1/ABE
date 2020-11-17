// @file palisade-tbo-linear.cpp - Example of Token-Based Obfuscation for linear
// functions
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

#include <iostream>
#include "obfuscation/lwetbolinearsecret.h"

#include "utils/debug.h"

using namespace lbcrypto;

shared_ptr<vector<NativeInteger>> BuildWeightVector(
    const vector<uint32_t> &thresholds, PlaintextModulus p, uint32_t N,
    uint32_t wordSize) {
  DiscreteUniformGeneratorImpl<NativeVector> dug;
  dug.SetModulus(p);

  auto weights = std::make_shared<vector<NativeInteger>>(N);

  for (size_t k = 0; k < thresholds.size(); k++)
    for (size_t i = 0; i < wordSize - 1; i++) {
      if (i > thresholds[k])
        (*weights)[i + k * wordSize] = 0;
      else
        (*weights)[i + k * wordSize] = dug.GenerateInteger();
    }

  return weights;
}

shared_ptr<vector<uint32_t>> BuildDataVector(const vector<uint32_t> &input,
                                             uint32_t wordSize) {
  auto indices = std::make_shared<vector<uint32_t>>(input.size());
  for (size_t k = 0; k < input.size(); k++)
    (*indices)[k] = input[k] + k * wordSize;

  return indices;
}

int main(int argc, char *argv[]) {
  TimeVar t;

  double processingTime(0.0);

  uint32_t N = 1280;
  uint32_t n = 2048;
  uint32_t numAtt = 5;
  PlaintextModulus p = 1099511627776;  // 2^40
  uint32_t wordSize = 256;

  // vector of thresholds
  vector<uint32_t> thresholds = {134, 90, 56, 89, 200};
  shared_ptr<vector<NativeInteger>> weights =
      BuildWeightVector(thresholds, p, N, wordSize);

  vector<vector<uint32_t>> inputs = {
      {135, 100, 60, 95, 210}, {34, 100, 60, 95, 210},
      {135, 80, 60, 95, 210},  {135, 100, 40, 95, 210},
      {135, 100, 60, 85, 210}, {255, 254, 200, 150, 215},
  };

  std::cout << "\n======================================================"
            << std::endl;
  std::cout << "OBFUSCATION PROTOTYPE USING PRECOMPUTED SECRET KEYS"
            << std::endl;
  std::cout << "========================================================"
            << std::endl;

  TIC(t);
  LWETBOLinearSecret algorithm(N, n, p, numAtt);
  processingTime = TOC_US(t);
  std::cout << "Parameter Generation: " << processingTime / 1000 << "ms"
            << std::endl;

  std::cout << "\nn = " << algorithm.GetSecurityParameter() << std::endl;
  std::cout << "log2 q = " << algorithm.GetLogModulus() << std::endl;
  std::cout << "Number of attributes = " << algorithm.GetNumAtt() << std::endl;
  std::cout << "plaintext modulus = " << algorithm.GetPlaintextModulus()
            << std::endl;
  std::cout << "Dimension of weight/data vectors = " << algorithm.GetDimension()
            << std::endl;

  TIC(t);
  shared_ptr<LWETBOKeys> keys = algorithm.KeyGen();
  processingTime = TOC_US(t);
  std::cout << "\nKey generation time: " << processingTime / 1000 << "ms"
            << std::endl;

  std::cout << "MODE = " << keys->GetMode() << std::endl;

  TIC(t);
  shared_ptr<NativeVector> ciphertext = algorithm.Obfuscate(keys, *weights);
  processingTime = TOC_US(t);
  std::cout << "\nObfuscation time: " << processingTime / 1000 << "ms"
            << std::endl;

  bool errorFlag = false;

  double evalTokenTime(0.0);
  double evalTime(0.0);

  std::cout << "\nThresholds vector: " << thresholds << std::endl;

  for (size_t i = 0; i < inputs.size(); i++) {
    shared_ptr<vector<uint32_t>> indices = BuildDataVector(inputs[i], wordSize);

    std::cout << "\nInput #" << i + 1 << ": " << inputs[i] << std::endl;

    TIC(t);
    shared_ptr<NativeVector> token = algorithm.TokenGen(keys, *indices);
    processingTime = TOC_US(t);
    evalTokenTime += processingTime;
    std::cout << "Token generation time: " << processingTime / 1000 << "ms"
              << std::endl;

    TIC(t);
    NativeInteger result = algorithm.EvaluateClassifier(
        *indices, ciphertext, keys->GetPublicRandomVector(),
        keys->GetPublicRandomVectorPrecon(), token);
    processingTime = TOC_US(t);
    evalTime += processingTime;
    std::cout << "Evaluation time: " << processingTime / 1000 << "ms"
              << std::endl;

    std::cout << "result (encrypted computation) = " << result << std::endl;

    TIC(t);
    NativeInteger resultClear =
        algorithm.EvaluateClearClassifier(*indices, *weights);
    processingTime = TOC_US(t);
    std::cout << "Evaluation time (in clear): " << processingTime / 1000 << "ms"
              << std::endl;

    std::cout << "result (plaintext computation) = " << resultClear
              << std::endl;

    string flagResult;

    if (resultClear == 0)
      flagResult = "MATCH: ";
    else
      flagResult = "NO MATCH: ";

    if (result == resultClear) {
      flagResult = flagResult + " CORRECT";
    } else {
      flagResult = flagResult + " INCORRECT";
      errorFlag = true;
    }

    std::cout << flagResult << std::endl;
  }

  std::cout << "\nT: Average token evaluation time:  "
            << "\t" << evalTokenTime / static_cast<double>(1000 * inputs.size())
            << " ms" << std::endl;
  std::cout << "T: Average evaluation time:  "
            << "\t" << evalTime / static_cast<double>(1000 * inputs.size())
            << " ms" << std::endl;
  std::cout << "T: Average total evaluation time:       "
            << "\t"
            << (evalTokenTime + evalTime) /
                   static_cast<double>(1000 * inputs.size())
            << " ms" << std::endl;

  std::cout << "\n============================================================="
               "========================"
            << std::endl;
  std::cout << "OBFUSCATION PROTOTYPE WITH ON-DEMAND GENERATION OF SECRET KEYS "
               "USING AES-CTR"
            << std::endl;
  std::cout << "==============================================================="
               "========================"
            << std::endl;

  unsigned char key[32] = {1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11,
                           12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22,
                           23, 24, 25, 26, 27, 28, 29, 30, 31, 32};

  TIC(t);
  LWETBOLinearSecret algorithmAES(N, n, p, numAtt);
  processingTime = TOC_US(t);
  std::cout << "Parameter Generation: " << processingTime / 1000 << "ms"
            << std::endl;

  std::cout << "\nn = " << algorithmAES.GetSecurityParameter() << std::endl;
  std::cout << "log2 q = " << algorithmAES.GetLogModulus() << std::endl;
  std::cout << "Number of attributes = " << algorithmAES.GetNumAtt()
            << std::endl;
  std::cout << "plaintext modulus = " << algorithmAES.GetPlaintextModulus()
            << std::endl;
  std::cout << "Dimension of weight/data vectors = "
            << algorithmAES.GetDimension() << std::endl;

  TIC(t);
  shared_ptr<LWETBOKeys> keysAES = algorithmAES.KeyGen(key, 1);
  processingTime = TOC_US(t);
  std::cout << "\nKey generation time: " << processingTime / 1000 << "ms"
            << std::endl;

  std::cout << "MODE = " << keysAES->GetMode() << std::endl;

  TIC(t);
  shared_ptr<NativeVector> ciphertextAES =
      algorithmAES.Obfuscate(keysAES, *weights);
  processingTime = TOC_US(t);
  std::cout << "\nObfuscation time: " << processingTime / 1000 << "ms"
            << std::endl;

  evalTokenTime = 0.0;
  evalTime = 0.0;

  std::cout << "\nThresholds vector: " << thresholds << std::endl;

  for (size_t i = 0; i < inputs.size(); i++) {
    shared_ptr<vector<uint32_t>> indices = BuildDataVector(inputs[i], wordSize);

    std::cout << "\nInput #" << i + 1 << ": " << inputs[i] << std::endl;

    TIC(t);
    shared_ptr<NativeVector> tokenAES =
        algorithmAES.TokenGen(keysAES, *indices);
    processingTime = TOC_US(t);
    evalTokenTime += processingTime;
    std::cout << "Token generation time: " << processingTime / 1000 << "ms"
              << std::endl;

    TIC(t);
    NativeInteger result = algorithmAES.EvaluateClassifier(
        *indices, ciphertextAES, keysAES->GetPublicRandomVector(),
        keysAES->GetPublicRandomVectorPrecon(), tokenAES);
    processingTime = TOC_US(t);
    evalTime += processingTime;
    std::cout << "Evaluation time: " << processingTime / 1000 << "ms"
              << std::endl;

    std::cout << "result (encrypted computation) = " << result << std::endl;

    TIC(t);
    NativeInteger resultClear =
        algorithmAES.EvaluateClearClassifier(*indices, *weights);
    processingTime = TOC_US(t);
    std::cout << "Evaluation time (in clear): " << processingTime / 1000 << "ms"
              << std::endl;

    std::cout << "result (plaintext computation) = " << resultClear
              << std::endl;

    string flagResult;

    if (resultClear == 0)
      flagResult = "MATCH: ";
    else
      flagResult = "NO MATCH: ";

    if (result == resultClear) {
      flagResult = flagResult + " CORRECT";
    } else {
      flagResult = flagResult + " INCORRECT";
      errorFlag = true;
    }

    std::cout << flagResult << std::endl;
  }

  std::cout << "\nT: Average token evaluation time:  "
            << "\t" << evalTokenTime / static_cast<double>(1000 * inputs.size())
            << " ms" << std::endl;
  std::cout << "T: Average evaluation time:  "
            << "\t" << evalTime / static_cast<double>(1000 * inputs.size())
            << " ms" << std::endl;
  std::cout << "T: Average total evaluation time:       "
            << "\t"
            << (evalTokenTime + evalTime) /
                   static_cast<double>(1000 * inputs.size())
            << " ms" << std::endl;

  return errorFlag;
}
