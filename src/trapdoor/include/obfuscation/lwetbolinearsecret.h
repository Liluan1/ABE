// @file lwetbolinearsecret.h Implementation of token-based obfuscation of
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

#ifndef LBCRYPTO_OBFUSCATE_LWETBOLINEARSECRET_H
#define LBCRYPTO_OBFUSCATE_LWETBOLINEARSECRET_H

#include <cmath>
#include <memory>
#include <vector>
using std::shared_ptr;
#include <string>
using std::string;

#include "math/backend.h"
#include "math/distrgen.h"
#include "utils/aesutil.h"
#include "utils/inttypes.h"

/**
 * @namespace lbcrypto
 * The namespace of lbcrypto
 */
namespace lbcrypto {

enum TboLinearMode { PRECOMPUTED = 0, AES = 1 };

inline std::ostream &operator<<(std::ostream &s, TboLinearMode mode) {
  switch (mode) {
    case PRECOMPUTED:
      s << "PRECOMPUTED";
      break;
    case AES:
      s << "AES";
      break;
    default:
      s << "UKNOWN";
      break;
  }
  return s;
}

class LWETBOKeys {
 public:
  // constructor for the PRECOMPUTED scenario
  explicit LWETBOKeys(const vector<shared_ptr<NativeVector>> &secretKey,
                      const NativeVector &publicKey,
                      const NativeVector &publicKeyPrecon)
      : m_mode(PRECOMPUTED),
        m_secretKey(secretKey),
        m_publicRandomVector(publicKey),
        m_publicRandomVectorPrecon(publicKeyPrecon),
        m_AESkey(nullptr),
        m_seed(0),
        m_n(0),
        m_modulus(0) {}

  // constructor for the AES PRNG case (secret keys are computed on the fly)
  explicit LWETBOKeys(const NativeVector &publicKey,
                      const NativeVector &publicKeyPrecon,
                      unsigned char *aes_key, uint32_t seed, size_t n,
                      const NativeInteger &modulus)
      : m_mode(AES),
        m_publicRandomVector(publicKey),
        m_publicRandomVectorPrecon(publicKeyPrecon),
        m_AESkey(aes_key),
        m_seed(seed),
        m_n(n),
        m_modulus(modulus) {}

  // either accesses a precomputed secret key or computes it using AES on the
  // fly
  const shared_ptr<NativeVector> GetSecretKey(size_t index) const {
    shared_ptr<NativeVector> v;

    if (m_mode == AES) {
      v = std::make_shared<NativeVector>(NativeVector(m_n, m_modulus));
      unsigned char iv[4] = {0, 1, 2, 3};
      unsigned char counter[16];
      unsigned char result[16];
      int64_t i1, i2;
      AESUtil util(iv, m_AESkey, 32);

      for (size_t i = 0; i < m_n; i += 2) {
        int ctr = m_seed + index * (m_n / 2) + 2 * i;
        util.SplitIntegers(counter, ctr, ctr + 1);
        util.EncryptBlock(counter, result);
        util.CombineBytes(result, i1, i2);
        (*v)[i] = i1 % m_modulus;
        (*v)[i + 1] = i2 % m_modulus;
      }
    } else {
      v = m_secretKey[index];
    }

    return v;
  }

  const NativeVector &GetPublicRandomVector() { return m_publicRandomVector; }

  const NativeVector &GetPublicRandomVectorPrecon() {
    return m_publicRandomVectorPrecon;
  }

  const TboLinearMode GetMode() { return m_mode; }

 private:
  // mode of extracting secret keys: from memory or using AES-CTR on demand
  TboLinearMode m_mode;

  // secret keys for the precomputed scenario
  vector<shared_ptr<NativeVector>> m_secretKey;

  // public key and a precomputed vector for it for the preconditioned NTL
  // modular arithmetic subroutine
  NativeVector m_publicRandomVector;
  NativeVector m_publicRandomVectorPrecon;

  // AES key (nonce)
  unsigned char *m_AESkey;

  // initial value of the counter to be encrypted
  uint32_t m_seed;

  // LWE security parameter (vector size)
  size_t m_n;

  // ciphertext modulus
  NativeInteger m_modulus;
};

/**
 * @brief LWE-based token-based obfuscation of linear functions
 */
class LWETBOLinearSecret {
 public:
  /**
   * Constructor
   *
   * @param N the dimension
   * @param n LWE security parameter
   * @param wmax infinity norm of the weights vector
   * @param pmax infinity norm of input data vector
   * @param numAtt number of attributes
   */
  explicit LWETBOLinearSecret(uint32_t N, uint32_t n, uint32_t wmax,
                              uint32_t pmax, uint32_t numAtt);

  /**
   * Constructor
   *
   * @param N the dimension
   * @param n LWE security parameter
   * @param numAtt in the case of classification
   */
  explicit LWETBOLinearSecret(uint32_t N, uint32_t n, PlaintextModulus p,
                              uint32_t numAtt);

  /**
   * Gets the LWE security parameter
   * @return the LWE security parameter
   */
  uint32_t GetSecurityParameter() const { return m_n; }

  /**
   * Gets the log of the modulus
   * @return the log of the modulus
   */
  uint32_t GetLogModulus() const;

  /**
   * Gets the number of attributes (in classification)
   * @return the number of attributes
   */
  uint32_t GetNumAtt() const { return m_numAtt; }

  /**
   * Gets the weight infinity norm
   * @return the weight norm
   */
  uint32_t GetWeightNorm() const { return m_wmax; }

  /**
   * Gets the "plaintext" modulus p used by the LWE scheme
   * @return p
   */
  PlaintextModulus GetPlaintextModulus() const { return m_p; }

  /**
   * Gets "ciphertext" modulus for the LWE problem
   * @return the ciphertext modulus q
   */
  NativeInteger GetModulus() const { return m_modulus; }

  /**
   * Gets the dimension N
   * @return the dimension N
   */
  uint32_t GetDimension() const { return m_N; }

  /**
   * Generate N random secret vectors Z_q^n and public random vector a - AES
   * mode
   * @return the secret keys and public random vector
   */
  shared_ptr<LWETBOKeys> KeyGen(unsigned char *aes_key, uint32_t seed) const;

  /**
   * Generate N random secret vectors Z_q^n and public random vector a -
   * PRECOMPUTED mode
   * @return the secret keys and public random vector
   */
  shared_ptr<LWETBOKeys> KeyGen() const;

  /**
   * Generate token t = \Sum{w_i s_i} \in Z_q^n for a general case of linear
   * function
   *
   * @param keys secret keys
   * @param input input data vector
   * @return the token
   */
  shared_ptr<NativeVector> TokenGen(
      const vector<shared_ptr<NativeVector>> &keys,
      const vector<NativeInteger> &input) const;

  /**
   * Generate token t = \Sum{x_i s_i} \in Z_q^n for the case when x_i's are 1's
   *
   * @param keys secret keys
   * @param inputIndices indicies where x_i is 1 (0 elsewhere)
   * @return the token
   */
  shared_ptr<NativeVector> TokenGen(shared_ptr<LWETBOKeys> &keys,
                                    const vector<uint32_t> &inputIndices) const;

  /**
   * Generates an encryption of weights (obfuscated program)
   *
   * @param keyPair secret keys + public random vector
   * @param weights the weights vector
   * @return the obfuscated program
   */
  shared_ptr<NativeVector> Obfuscate(
      const shared_ptr<LWETBOKeys> keyPair,
      const vector<NativeInteger> &weights) const;

  /**
   * Evaluates \Sum{w_i x_i} using obfuscated program
   *
   * @param inputIndices indices where input is 1
   * @param ciphertext obfuscated program
   * @param publicRandomVector public random vector
   * @param publicRandomVectorPrecon precomputation for NTL
   * @param token the token for the input data vector
   * @return the result of the summation
   */
  NativeInteger EvaluateClassifier(const vector<uint32_t> &inputIndices,
                                   const shared_ptr<NativeVector> ciphertext,
                                   const NativeVector &publicRandomVector,
                                   const NativeVector &publicRandomVectorPrecon,
                                   const shared_ptr<NativeVector> token) const;

  /**
   * Evaluates \Sum{w_i x_i} using cleartext program
   *
   * @param inputIndices indices where input is 1
   * @param weigts the weights vector
   * @return the result of the summation
   */
  NativeInteger EvaluateClearClassifier(
      const vector<uint32_t> &inputIndices,
      const vector<NativeInteger> weights) const;

 private:
  // Dimension - size of weight/data vectors
  uint32_t m_N;

  // LWE security parameter
  uint32_t m_n;

  // Infinity norm for the weight vector
  uint32_t m_wmax;

  // Plaintext modulus p
  PlaintextModulus m_p;

  // Infinity norm for the input data vector
  uint32_t m_pmax;

  // LWE "ciphertext" modulus q
  NativeInteger m_modulus;

  // number of attributes (when applied to classification)
  uint32_t m_numAtt;

  // Discrete Gaussian distribution for generating the noise in the LWE
  // encryption
  DiscreteGaussianGeneratorImpl<NativeVector> m_dgg;

  /**
   * Method to estimate the modulus
   * Used as a subroutine by constructor LWETBOLinearSecret(uint32_t N, uint32_t
   * n, uint32_t wmax, uint32_t pmax)
   *
   * @return estimated value q of modulus
   */
  double EstimateModulus();

  /**
   * Method to estimate the modulus for the classifier scenario
   * Used as a subroutine by constructor LWETBOLinearSecret(uint32_t N, uint32_t
   * n, uint32_t numAtt);
   *
   * @return estimated value q of modulus
   */
  double EstimateModulusClassifier();
};

}  // namespace lbcrypto

#endif
