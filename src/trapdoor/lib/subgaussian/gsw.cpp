// @file gsw.cpp Provides implementation of the GSW variant described in
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

#ifndef _LBCRYPTO_SUBGAUSSIAN_GSW_CPP
#define _LBCRYPTO_SUBGAUSSIAN_GSW_CPP

#include "subgaussian/gsw.h"

namespace lbcrypto {

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
shared_ptr<GSWSecretKey<Integer>> GSWScheme<Integer, Vector>::SecretKeyGen()
    const {
  auto dgg_allocator = [&]() {
    return m_cryptoParams.GetDgg().GenerateInteger(m_cryptoParams.GetModulus());
  };
  auto zero_allocator = [&]() { return Integer(0); };

  return std::make_shared<GSWSecretKey<Integer>>(
      zero_allocator, m_cryptoParams.Getn() - 1, 1, dgg_allocator);
}

template <class Integer, class Vector>
shared_ptr<GSWCiphertext<Integer>> GSWScheme<Integer, Vector>::Encrypt(
    const GSWPlaintext<Integer> &plaintext,
    const shared_ptr<GSWSecretKey<Integer>> sk) const {
  auto dgg_allocator = [&]() {
    return m_cryptoParams.GetDgg().GenerateInteger(m_cryptoParams.GetModulus());
  };
  auto uniform_allocator = [&]() {
    return m_cryptoParams.GetDug().GenerateInteger();
  };
  auto zero_allocator = [&]() { return Integer(0); };

  const Integer &modulus = m_cryptoParams.GetModulus();

  Matrix<Integer> cbar(zero_allocator, m_cryptoParams.Getn() - 1,
                       m_cryptoParams.Getm(), uniform_allocator);
  Matrix<Integer> et(zero_allocator, 1, m_cryptoParams.Getm(), dgg_allocator);
  Matrix<Integer> skt = sk->Transpose();
  Matrix<Integer> bt = et.ModSubEq((skt.Mult(cbar)).ModEq(modulus), modulus);
  cbar.VStack(bt);

  Matrix<Integer> g(zero_allocator, m_cryptoParams.Getn(),
                    m_cryptoParams.Getm());
  g = g.GadgetVector(m_cryptoParams.GetBase());

  std::cout << "n = " << m_cryptoParams.Getn() << std::endl;
  std::cout << "m = " << m_cryptoParams.Getm() << std::endl;

  auto c = std::make_shared<Matrix<Integer>>(
      cbar + (g.ScalarMult(plaintext)).ModEq(modulus).ModEq(modulus));

  return c;
}

template <class Integer, class Vector>
GSWPlaintext<Integer> GSWScheme<Integer, Vector>::Decrypt(
    const shared_ptr<GSWCiphertext<Integer>> ciphertext,
    const shared_ptr<GSWSecretKey<Integer>> sk) const {
  const Integer &modulus = m_cryptoParams.GetModulus();
  GSWPlaintext<Integer> mu = Integer(0);
  for (size_t i = 0; i < m_cryptoParams.Getn() - 1; i++) {
    mu += (*sk)(i, 0).ModMul((*ciphertext)(i, m_cryptoParams.Getm() - 2),
                             modulus);
  }
  mu += ((*ciphertext)(m_cryptoParams.Getn() - 1, m_cryptoParams.Getm() - 2));

  mu.ModEq(modulus);

  Integer half = modulus >> 1;
  if (mu > half) mu = modulus - mu;

  return mu.MultiplyAndRound(m_cryptoParams.GetBase(),
                             m_cryptoParams.GetModulus());
}

template <class Integer, class Vector>
shared_ptr<GSWCiphertext<Integer>> GSWScheme<Integer, Vector>::EvalAdd(
    const shared_ptr<GSWCiphertext<Integer>> ct1,
    const shared_ptr<GSWCiphertext<Integer>> ct2) {
  const Integer &modulus = m_cryptoParams.GetModulus();
  auto c = std::make_shared<Matrix<Integer>>((*ct1 + *ct2).ModEq(modulus));

  return c;
}

template <class Integer, class Vector>
shared_ptr<GSWCiphertext<Integer>> GSWScheme<Integer, Vector>::EvalMult(
    const shared_ptr<GSWCiphertext<Integer>> ct1,
    const shared_ptr<GSWCiphertext<Integer>> ct2) {
  const Integer &modulus = m_cryptoParams.GetModulus();
  // std::cout << "before inversion" << std::endl;
  auto ct2Inverse = InverseG(ct2);
  // std::cout << "after inversion" << std::endl;
  auto temp = ct1->Mult(*ct2Inverse);
  // std::cout << "after multiplication" << std::endl;
  auto c = std::make_shared<Matrix<Integer>>(temp.ModEq(modulus));
  // std::cout << "after modulo operation" << std::endl;

  return c;
}

template <class Integer, class Vector>
shared_ptr<GSWCiphertext<Integer>> GSWScheme<Integer, Vector>::InverseG(
    const shared_ptr<GSWCiphertext<Integer>> ct) {
  size_t cols = ct->GetCols();
  size_t rows = ct->GetRows();
  Matrix<Integer> gInverse([&]() { return Integer(0); }, cols, cols);

  std::cout << "cols = " << cols << std::endl;
  std::cout << "rows = " << rows << std::endl;

  if (m_cryptoParams.GetMode() == DETERMINISTIC) {
    for (size_t i = 0; i < cols; i++) {
      for (size_t j = 0; j < rows; j++) {
        auto vector = GetDigits((*ct)(j, i), m_cryptoParams.GetBase(),
                                m_cryptoParams.Getl());
        for (size_t k = 0; k < vector->size(); k++) {
          gInverse(j * m_cryptoParams.Getl() + k, i) = (*vector)[k];
        }
      }
    }
  } else {  // RADOMIZED ROUNDING
    const Integer &q = m_cryptoParams.GetModulus();

    for (size_t i = 0; i < cols; i++) {
      for (size_t j = 0; j < rows; j++) {
        vector<int64_t> vector(m_cryptoParams.Getl());

        m_cryptoParams.GetSampler().InverseG(
            (*ct)(j, i), PseudoRandomNumberGenerator::GetPRNG(), &vector);

        for (size_t k = 0; k < vector.size(); k++) {
          if (vector[k] > 0)
            gInverse(j * m_cryptoParams.Getl() + k, i) = vector[k];
          else
            gInverse(j * m_cryptoParams.Getl() + k, i) =
                q - Integer(-vector[k]);
        }
      }
    }
  }

  return std::make_shared<Matrix<Integer>>(gInverse);
}

}  // namespace lbcrypto

#endif
