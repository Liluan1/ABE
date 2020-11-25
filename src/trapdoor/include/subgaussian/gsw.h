// @file gsw.h Provides implementation of the GSW variant described in
// https://eprint.iacr.org/2014/094
// @author TPOC: contact@palisade-crypto.org
//
// @copyright Copyright (c) 2018, New Jersey Institute of Technology (NJIT)
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

#include "lattice/trapdoor.h"
#include "math/discretegaussiangenerator.h"
#include "math/matrix.h"
#include "subgaussian.h"

#ifndef LBCRYPTO_SUBGAUSSIAN_GSW_H
#define LBCRYPTO_SUBGAUSSIAN_GSW_H

namespace lbcrypto {

enum GINVERSEMODE { DETERMINISTIC = 0, RANDOM = 1 };

// dimension Z_q^n-1
template <class Integer>
using GSWSecretKey = Matrix<Integer>;

// dimension Z_q^m \times n
template <class Integer>
using GSWPublicKey = Matrix<Integer>;

// dimension Z_q^n \times n l
template <class Integer>
using GSWCiphertext = Matrix<Integer>;

// dimension Z_q
template <class Integer>
using GSWPlaintext = Integer;

template <class Integer, class Vector>
class GSWCryptoParameters {
 public:
  GSWCryptoParameters()
      : m_n(0),
        m_l(0),
        m_m(0),
        m_base(0),
        m_modulus(Integer(0)),
        m_mode(DETERMINISTIC),
        m_sampler(LatticeSubgaussianUtility<Integer>()) {}

  GSWCryptoParameters(uint32_t n, int64_t base, const Integer &q, double std,
                      GINVERSEMODE mode = DETERMINISTIC)
      : m_n(n), m_base(base), m_modulus(q), m_mode(mode) {
    m_dgg.SetStd(std);
    m_dug.SetModulus(q);
    m_l = std::ceil(q.GetMSB() / log2(base));
    m_m = m_n * m_l;
    m_sampler = LatticeSubgaussianUtility<Integer>(base, q, m_l);
  }

  const DiscreteGaussianGeneratorImpl<Vector> &GetDgg() const { return m_dgg; }

  const DiscreteUniformGeneratorImpl<Vector> &GetDug() const { return m_dug; }

  const LatticeSubgaussianUtility<Integer> &GetSampler() const {
    return m_sampler;
  }

  uint32_t Getn() const { return m_n; }

  uint32_t Getl() const { return m_l; }

  uint32_t Getm() const { return m_m; }

  GINVERSEMODE GetMode() const { return m_mode; }

  int64_t GetBase() const { return m_base; }

  const Integer &GetModulus() const { return m_modulus; }

 private:
  uint32_t m_n;
  uint32_t m_l;
  uint32_t m_m;
  int64_t m_base;
  Integer m_modulus;
  DiscreteGaussianGeneratorImpl<Vector> m_dgg;
  DiscreteUniformGeneratorImpl<Vector> m_dug;

  GINVERSEMODE m_mode;
  LatticeSubgaussianUtility<Integer> m_sampler;
};

template <class Integer, class Vector>
class GSWScheme {
 public:
  void Setup(uint32_t n, uint32_t base, const Integer &q, double std,
             GINVERSEMODE mode = DETERMINISTIC) {
    m_cryptoParams =
        GSWCryptoParameters<Integer, Vector>(n, base, q, std, mode);
  }

  shared_ptr<GSWSecretKey<Integer>> SecretKeyGen() const;

  // shared_ptr<GSWPublicKey<Integer>> PublicKeyGen(const
  // shared_ptr<GSWSecretKey<Integer>>) const;

  shared_ptr<GSWCiphertext<Integer>> Encrypt(
      const GSWPlaintext<Integer> &plaintext,
      const shared_ptr<GSWSecretKey<Integer>>) const;

  GSWPlaintext<Integer> Decrypt(
      const shared_ptr<GSWCiphertext<Integer>> ciphertext,
      const shared_ptr<GSWSecretKey<Integer>>) const;

  shared_ptr<GSWCiphertext<Integer>> EvalAdd(
      const shared_ptr<GSWCiphertext<Integer>>,
      const shared_ptr<GSWCiphertext<Integer>>);

  shared_ptr<GSWCiphertext<Integer>> EvalMult(
      const shared_ptr<GSWCiphertext<Integer>>,
      const shared_ptr<GSWCiphertext<Integer>>);

 private:
  GSWCryptoParameters<Integer, Vector> m_cryptoParams;

  shared_ptr<GSWCiphertext<Integer>> InverseG(
      const shared_ptr<GSWCiphertext<Integer>>);
};

}  // namespace lbcrypto

#endif
