// @file palisade-tbo-gbp-32.cpp - Example of Token-Based Obfuscation for
// general branching programs
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
#include "obfuscation/lwebpchcprf.h"
#include "utils/parallel.h"

#include "utils/debug.h"

using namespace lbcrypto;

void test(const BPCHCPRF<DCRTPoly>& algorithm,
          const vector<vector<Matrix<int>>>& M,
          const vector<std::pair<string, bool>> cases) {
  TimeVar t;
  double processingTime;

  std::cout << "n = " << algorithm.GetRingDimension() << std::endl;
  std::cout << "log2 q = " << algorithm.GetLogModulus() << std::endl;

  TIC(t);
  auto key = algorithm.KeyGen();
  processingTime = TOC(t);
  std::cout << "Master Secret (Unconstrained) Key Generation: "
            << processingTime << "ms" << std::endl;

  TIC(t);
  auto constrainedKey = algorithm.Constrain(key, M);
  processingTime = TOC(t);
  std::cout << "Constrained Key Generation: " << processingTime << "ms"
            << std::endl;

  for (const auto& value : cases) {
    std::cout << "input: " << value.first << std::endl;
    TIC(t);
    const auto value1 = algorithm.Evaluate(key, value.first);
    processingTime = TOC(t);
    std::cout << "EvalToken: " << processingTime << "ms" << std::endl;
    TIC(t);
    const auto value2 = algorithm.Evaluate(constrainedKey, value.first);
    processingTime = TOC(t);
    std::cout << "Evaluation: " << processingTime << "ms" << std::endl;

    // cout << value1 << std::endl;
    // cout << value2 << std::endl;
    bool match = *value1 == *value2;
    std::cout << (match ? "Matched " : "Did not match ")
              << (match == value.second ? "(Correct)" : "(Incorrect)")
              << std::endl;
  }
}

void CC17Manual() {
  // M accepts 010? and 100?
  auto zero_alloc = []() { return 0; };
  Matrix<int> M_00(zero_alloc, 3, 3);
  M_00(0, 2) = 1;
  M_00(1, 1) = 1;
  M_00(2, 0) = 1;
  Matrix<int> M_01(zero_alloc, 3, 3);
  M_01(0, 2) = 1;
  M_01(1, 0) = 1;
  M_01(2, 1) = 1;
  Matrix<int> M_10(zero_alloc, 3, 3);
  M_10(0, 2) = 1;
  M_10(1, 1) = 1;
  M_10(2, 0) = 1;
  Matrix<int> M_11(zero_alloc, 3, 3);
  M_11(0, 1) = 1;
  M_11(1, 2) = 1;
  M_11(2, 0) = 1;
  Matrix<int> M_20(zero_alloc, 3, 3);
  M_20(0, 0) = 1;
  M_20(1, 2) = 1;
  M_20(2, 1) = 1;
  Matrix<int> M_21(zero_alloc, 3, 3);
  M_21(0, 1) = 1;
  M_21(1, 0) = 1;
  M_21(2, 2) = 1;
  Matrix<int> M_30(zero_alloc, 3, 3);
  M_30(0, 0) = 1;
  M_30(1, 1) = 1;
  M_30(2, 2) = 1;
  Matrix<int> M_31(zero_alloc, 3, 3);
  M_31(0, 0) = 1;
  M_31(1, 1) = 1;
  M_31(2, 2) = 1;
  vector<vector<Matrix<int>>> M = {
      {M_00, M_01}, {M_10, M_11}, {M_20, M_21}, {M_30, M_31}};

  CC17Algorithm<DCRTPoly> algorithm(1 << 15, 2, 4, 1024, 3);
  test(algorithm, M, {{"1001", true}, {"0110", false}});
}

void CVW18Disjunction(const string& pattern,
                      const vector<std::pair<string, bool>>& cases) {
  // "10*0" -> x0 V -x1 V -x3
  auto zero_alloc = []() { return 0; };
  Matrix<int> I(zero_alloc, 1, 1);
  I(0, 0) = 1;
  Matrix<int> N(zero_alloc, 1, 1);
  Matrix<int> v(zero_alloc, 1, 1);
  v(0, 0) = 1;

  vector<vector<Matrix<int>>> M;
  for (const char& value : pattern) {
    if (value == '*') {
      M.push_back({I, I});
    } else if (value == '1') {
      M.push_back({I, N});
    } else {
      M.push_back({N, I});
    }
  }

  CVW18Algorithm<DCRTPoly> algorithm(1 << 15, 2, pattern.length(), 1024, v);
  test(algorithm, M, cases);
}

void CVW18HammingCloseness(const string& pattern, usint threshold,
                           const vector<std::pair<string, bool>>& cases) {
  // true if distance is smaller than threshold
  auto zero_alloc = []() { return 0; };
  Matrix<int> I(zero_alloc, threshold + 1, threshold + 1);
  for (usint i = 0; i <= threshold; i++) {
    I(i, i) = 1;
  }
  Matrix<int> N(zero_alloc, threshold + 1, threshold + 1);
  for (usint i = 0; i < threshold; i++) {
    N(i, i + 1) = 1;
  }
  N(threshold, threshold) = 1;
  Matrix<int> R(zero_alloc, threshold + 1, threshold + 1);
  R(threshold, threshold) = 1;
  vector<vector<Matrix<int>>> M;
  for (const char& value : pattern) {
    if (value == '0') {
      M.push_back({I, N});
    } else if (value == '1') {
      M.push_back({N, I});
    } else {
      M.push_back({I, I});
    }
  }
  M.back()[0] = M.back()[0] * R;
  M.back()[1] = M.back()[1] * R;
  Matrix<int> v(zero_alloc, 1, threshold + 1);
  v(0, 0) = 1;
  CVW18Algorithm<DCRTPoly> algorithm(1 << 20, 4, M.size(), 4096, v);
  test(algorithm, M, cases);
}

void CVW18WitnessEncryption() {
  WitnessEncryption<DCRTPoly> algorithm(1 << 15, 2, 1024, 4, 6);

  TimeVar t;
  double processingTime;

  std::cout << "n = " << algorithm.GetRingDimension() << std::endl;
  std::cout << "log2 q = " << algorithm.GetLogModulus() << std::endl;

  TIC(t);
  // clauses
  // -x0 V -x1 V -x2 V -x3
  // -x0 V -x1 V  x2
  // -x0 V  x1 V -x2 V  x3
  //  x0 V -x1 V -x2 V -x3
  //  x0 V -x1 V x2
  //  x0 V  x1
  auto ciphertext =
      algorithm.Encrypt({"0000", "001*", "0101", "1000", "101*", "11**"}, 0);
  processingTime = TOC(t);
  std::cout << "Encrypt: " << processingTime << "ms" << std::endl;

  TIC(t);
  const string input0 = "1001";
  const string input1 = "0010";
  usint value0 = algorithm.Decrypt(ciphertext, input0);
  usint value1 = algorithm.Decrypt(ciphertext, input1);
  processingTime = TOC(t);
  std::cout << "Decrypt: 2 * " << processingTime / 2 << "ms" << std::endl;
  std::cout << "input: " << input0 << std::endl;
  std::cout << value0 << (value0 == 0 ? " (Correct)" : " (Incorrect)")
            << std::endl;
  std::cout << "input: " << input1 << std::endl;
  std::cout << value1 << (value1 == 1 ? " (Correct)" : " (Incorrect)")
            << std::endl;
}

void CVW18CNF(const vector<vector<int>>& cnf,
              const vector<std::pair<string, bool>>& cases) {
  auto zero_alloc = []() { return 0; };
  Matrix<int> Good(zero_alloc, 3, 3);
  Good(0, 0) = 1;
  Good(1, 2) = 1;
  Good(2, 2) = 1;
  Matrix<int> Bad(zero_alloc, 3, 3);
  Bad(0, 0) = 1;
  Bad(1, 1) = 1;
  Bad(2, 2) = 1;
  Matrix<int> R(zero_alloc, 3, 3);
  R(0, 0) = 1;
  R(1, 0) = 1;
  R(2, 1) = 1;
  Matrix<int> T(zero_alloc, 3, 3);
  T(0, 0) = 1;
  vector<vector<Matrix<int>>> M;
  for (const auto& clause : cnf) {
    for (int variable : clause) {
      if (variable < 0) {
        M.push_back({Good, Bad});
      } else if (variable > 0) {
        M.push_back({Bad, Good});
      }
    }
    M.back()[0] = M.back()[0] * R;
    M.back()[1] = M.back()[1] * R;
  }
  M.back()[0] = M.back()[0] * T;
  M.back()[1] = M.back()[1] * T;
  Matrix<int> v(zero_alloc, 1, 3);
  v(0, 1) = 1;
  vector<std::pair<string, bool>> adjustedCases;
  for (const auto& c : cases) {
    string input = "";
    for (const auto& clause : cnf) {
      for (int variable : clause) {
        input.push_back(c.first[(variable > 0 ? variable : -variable) - 1]);
      }
    }
    adjustedCases.push_back({input, c.second});
  }
  CVW18Algorithm<DCRTPoly> algorithm(1 << 15, 2,
                                     adjustedCases[0].first.length(), 1024, v);
  test(algorithm, M, adjustedCases);
}

int main(int argc, char* argv[]) {
  PalisadeParallelControls.Enable();

  // CVW18Disjunction("10*000*1", {{"00111110", true}, {"01011100", false}});
  CVW18HammingCloseness(
      "0*100*10", 3,
      {{"10100110", true}, {"11100110", true}, {"11110111", false}});
  // CVW18WitnessEncryption();
  // CVW18CNF({{1, -2, 3}, {-1, 4, 5}}, {{"00000", true}, {"11000", false},
  // {"11001", true}});

  // CC17Manual();
  // CVW18Disjunction("10*000*110*000*110*000*110*000*1",
  // {{"00111110001111100011111000111110", true},
  // {"01011100010111000101110001011100", false}});
}
