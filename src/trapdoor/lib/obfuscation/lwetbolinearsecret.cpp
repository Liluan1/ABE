// @file lwetbolinearsecret.cpp Implementation of token-based obfuscation of
// linear functions (secret-key version)
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

#ifndef LBCRYPTO_OBFUSCATE_LWETBOLINEARSECRET_CPP
#define LBCRYPTO_OBFUSCATE_LWETBOLINEARSECRET_CPP

#include "obfuscation/lwetbolinearsecret.h"

namespace lbcrypto {

uint32_t LWETBOLinearSecret::GetLogModulus() const {
  double q = m_modulus.ConvertToDouble();
  uint32_t logModulus = floor(log2(q - 1.0) + 1.0);
  return logModulus;
}

LWETBOLinearSecret::LWETBOLinearSecret(uint32_t N, uint32_t n, uint32_t wmax,
                                       uint32_t pmax, uint32_t numAtt)
    : m_N(N), m_n(n), m_wmax(wmax), m_pmax(pmax), m_numAtt(numAtt) {
  double q = EstimateModulus();

  m_modulus = FirstPrime<NativeInteger>(ceil(log2(q)), 2 * m_n);

  m_dgg.SetStd(3.2);
}

LWETBOLinearSecret::LWETBOLinearSecret(uint32_t N, uint32_t n,
                                       PlaintextModulus p, uint32_t numAtt)
    : m_N(N), m_n(n), m_p(p), m_numAtt(numAtt) {
  double q = EstimateModulusClassifier();

  m_modulus = FirstPrime<NativeInteger>(ceil(log2(q)), 2 * m_n);

  m_dgg.SetStd(3.2);
}

double LWETBOLinearSecret::EstimateModulus() {
  // distribution parameter
  double sigma = 3.2;

  // assurance measure
  double alpha = 144;

  // Bound of the Gaussian error
  double Berr = sigma * sqrt(alpha);

  m_p = m_N * m_wmax * m_pmax;

  // Correctness constraint
  return 4 * m_N * m_pmax * m_p * Berr;
}

double LWETBOLinearSecret::EstimateModulusClassifier() {
  // distribution parameter
  double sigma = 3.2;

  // assurance measure
  double alpha = 144;

  // Bound of the Gaussian error
  double Berr = sigma * sqrt(alpha);

  // Correctness constraint
  return 4 * m_numAtt * m_p * Berr;
}

shared_ptr<LWETBOKeys> LWETBOLinearSecret::KeyGen() const {
  DiscreteUniformGeneratorImpl<NativeVector> dug;
  dug.SetModulus(m_modulus);

  vector<shared_ptr<NativeVector>> secretKey(m_N);

  for (size_t i = 0; i < m_N; i++)
    secretKey[i] = std::make_shared<NativeVector>(dug.GenerateVector(m_n));

  NativeVector publicRandomVector = dug.GenerateVector(m_n);

  NativeVector publicRandomVectorPrecon(m_n, m_modulus);

  for (size_t i = 0; i < m_n; i++)
    publicRandomVectorPrecon[i] =
        publicRandomVector[i].PrepModMulConst(m_modulus);

  auto keys = std::make_shared<LWETBOKeys>(secretKey, publicRandomVector,
                                           publicRandomVectorPrecon);

  return keys;
}

shared_ptr<LWETBOKeys> LWETBOLinearSecret::KeyGen(unsigned char *aes_key,
                                                  uint32_t seed) const {
  DiscreteUniformGeneratorImpl<NativeVector> dug;
  dug.SetModulus(m_modulus);

  NativeVector publicRandomVector = dug.GenerateVector(m_n);

  NativeVector publicRandomVectorPrecon(m_n, m_modulus);

  for (size_t i = 0; i < m_n; i++)
    publicRandomVectorPrecon[i] =
        publicRandomVector[i].PrepModMulConst(m_modulus);

  auto keys =
      std::make_shared<LWETBOKeys>(publicRandomVector, publicRandomVectorPrecon,
                                   aes_key, seed, m_n, m_modulus);

  return keys;
}

shared_ptr<NativeVector> LWETBOLinearSecret::TokenGen(
    const vector<shared_ptr<NativeVector>> &keys,
    const vector<NativeInteger> &input) const {
  auto token = std::make_shared<NativeVector>(m_n, m_modulus);

  for (size_t Ni = 0; Ni < input.size(); Ni++)
    token->ModAddEq(keys[Ni]->ModMul(input[Ni]));

  return token;
}

shared_ptr<NativeVector> LWETBOLinearSecret::TokenGen(
    shared_ptr<LWETBOKeys> &keys, const vector<uint32_t> &inputIndices) const {
  auto token = std::make_shared<NativeVector>(m_n, m_modulus);

  vector<NativeVector> secretKey = vector<NativeVector>(inputIndices.size());

#pragma omp parallel for schedule(dynamic)
  for (size_t i = 0; i < inputIndices.size(); i++) {
    secretKey[i] = *keys->GetSecretKey(inputIndices[i]);
  }

  for (size_t i = 0; i < inputIndices.size(); i++) {
    token->ModAddEq(secretKey[i]);
  }

  return token;
}

shared_ptr<NativeVector> LWETBOLinearSecret::Obfuscate(
    const shared_ptr<LWETBOKeys> keyPair,
    const vector<NativeInteger> &weights) const {
  auto ciphertext = std::make_shared<NativeVector>(m_N, m_modulus);

  NativeInteger mu = m_modulus.ComputeMu();

#pragma omp parallel for schedule(dynamic)
  for (size_t Ni = 0; Ni < m_N; Ni++) {
    shared_ptr<NativeVector> secretKey = keyPair->GetSecretKey(Ni);

    for (size_t ni = 0; ni < m_n; ni++) {
      (*ciphertext)[Ni].ModAddEq(
          (*secretKey)[ni].ModMulFast(keyPair->GetPublicRandomVector()[ni],
                                      m_modulus, mu),
          m_modulus);
    }

    (*ciphertext)[Ni].ModAddFastEq(
        NativeInteger(m_p).ModMulFast(m_dgg.GenerateInteger(m_modulus),
                                      m_modulus, mu),
        m_modulus);

    (*ciphertext)[Ni].ModAddFastEq(weights[Ni], m_modulus);
  }

  return ciphertext;
}

NativeInteger LWETBOLinearSecret::EvaluateClassifier(
    const vector<uint32_t> &inputIndices,
    const shared_ptr<NativeVector> ciphertext,
    const NativeVector &publicRandomVector,
    const NativeVector &publicRandomVectorPrecon,
    const shared_ptr<NativeVector> token) const {
  NativeInteger result;

  for (size_t Ni = 0; Ni < inputIndices.size(); Ni++)
    result.ModAddFastEq((*ciphertext)[inputIndices[Ni]], m_modulus);

  for (size_t ni = 0; ni < m_n; ni++)
    result.ModSubEq(
        (*token)[ni].ModMulFastConst(publicRandomVector[ni], m_modulus,
                                     publicRandomVectorPrecon[ni]),
        m_modulus);

  NativeInteger halfQ(m_modulus >> 1);

  if (result > halfQ)
    result.ModSubEq(m_modulus, m_p);
  else
    result.ModEq(m_p);

  return result;
}

NativeInteger LWETBOLinearSecret::EvaluateClearClassifier(
    const vector<uint32_t> &inputIndices,
    const vector<NativeInteger> weights) const {
  NativeInteger result;

  for (size_t Ni = 0; Ni < inputIndices.size(); Ni++)
    result += weights[inputIndices[Ni]];

  result.ModEq(m_p);

  return result;
}

}  // namespace lbcrypto

#endif
