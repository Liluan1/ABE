// @file subgaussian.h Provides implementation of subgaussian sampling
// algorithms
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

#include <memory>
#include <random>
#include <vector>

#include "math/matrix.h"

#ifndef LBCRYPTO_LATTICE_SUBGAUSSIAN_H
#define LBCRYPTO_LATTICE_SUBGAUSSIAN_H

namespace lbcrypto {

template <class Integer>
class LatticeSubgaussianUtility {
 public:
  LatticeSubgaussianUtility()
      : m_base(2), m_modulus(1), m_k(1), m_baseDigits{1} {};

  LatticeSubgaussianUtility(const uint32_t &base, const Integer &modulus,
                            const uint32_t &k)
      : m_base(base), m_modulus(modulus), m_k(k) {
    Precompute();
  }

  void InverseG(const Integer &u, PRNG &prng, vector<int64_t> *output) const;

  void BcBD(const vector<float> &target, PRNG &prng, vector<int64_t> *x) const;

  const uint32_t GetK() const { return m_k; }

 private:
  void Precompute();

  // input parameters
  uint32_t m_base;
  Integer m_modulus;
  uint32_t m_k;

  uint32_t m_baseDigits;

  // precomputed tables
  vector<int64_t> m_qvec;
  vector<float> m_d;
};

template <class Element>
void InverseRingVector(
    const LatticeSubgaussianUtility<typename Element::Integer> &util,
    const shared_ptr<typename Element::Params> ilParams,
    const Matrix<Element> &pubElemB, uint32_t seed, Matrix<Element> *psi) {
  PALISADE_THROW(
      not_implemented_error,
      "InverseRingVector is not implemented for this polynomial type.");
}

template <class Element>
void InverseRingVectorSpecial(
    const LatticeSubgaussianUtility<typename Element::Integer> &util,
    const shared_ptr<typename Element::Params> ilParams,
    const Matrix<Element> &pubElemB, uint32_t seed, Matrix<DCRTPoly> *psi) {
  PALISADE_THROW(
      not_implemented_error,
      "InverseRingVector is not implemented for this polynomial type.");
}

shared_ptr<Matrix<DCRTPoly>> InverseRingVectorDCRT(
    const std::vector<LatticeSubgaussianUtility<NativeInteger>> &util,
    const Matrix<DCRTPoly> &pubElemB, uint32_t seed);

template <>
void InverseRingVector<Poly>(
    const LatticeSubgaussianUtility<typename Poly::Integer> &util,
    const shared_ptr<typename Poly::Params> ilParams,
    const Matrix<Poly> &pubElemB, uint32_t seed, Matrix<Poly> *psi);

template <>
void InverseRingVectorSpecial<Poly>(
    const LatticeSubgaussianUtility<typename Poly::Integer> &util,
    const shared_ptr<typename Poly::Params> ilParams,
    const Matrix<Poly> &pubElemB, uint32_t seed, Matrix<DCRTPoly> *psi);

}  // namespace lbcrypto

#endif
