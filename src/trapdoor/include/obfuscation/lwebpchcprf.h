// @file lwebpchcprf.h Implementation of constraint-hiding constrained PRFs for
// branching programs as described in https://eprint.iacr.org/2017/143.pdf and
// https://eprint.iacr.org/2018/360.pdf
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

#ifndef LBCRYPTO_OBFUSCATE_LWEBPCHCPRF_H
#define LBCRYPTO_OBFUSCATE_LWEBPCHCPRF_H

#include <cmath>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "lattice/elemparams.h"
#include "lattice/ildcrtparams.h"
#include "lattice/ilelement.h"
#include "lattice/ilparams.h"
#include "lattice/trapdoor.h"
#include "math/backend.h"
#include "math/distrgen.h"
#include "utils/inttypes.h"

/**
 * @namespace lbcrypto
 * The namespace of lbcrypto
 */
namespace lbcrypto {

/**
 * @brief base class of BP CHCPRF scheme
 * @tparam Element a ring element.
 */
template <class Element>
class BPCHCPRF {
 public:
  /**
   * Constructor
   * @param base base for G-sampling
   * @param chunkSize number of bits encoded by one encoding matrix
   * @param length the input length
   * @param n ring dimension
   * @param w GGH15 encoding width
   */
  explicit BPCHCPRF(usint base, usint chunkSize, usint length, usint n,
                    usint w);

  /**
   * Gets the ring dimension
   * @return the ring dimension
   */
  usint GetRingDimension() const;

  /**
   * Gets the log of the modulus
   * @return the log of the modulus
   */
  usint GetLogModulus() const;

  /**
   * Method to generate unconstrained PRF key
   * @return unconstrained PRF key
   */
  const std::pair<vector<vector<Element>>, Matrix<Element>> KeyGen() const;

  /**
   * Method to constrain key by branching program
   * @param key unconstrained PRF key
   * @param M branching program
   * @return constrained key
   */
  const std::pair<Matrix<Element>, vector<vector<Matrix<Element>>>> Constrain(
      const std::pair<vector<vector<Element>>, Matrix<Element>>& key,
      const vector<vector<Matrix<int>>>& M) const;

  /**
   * Method to evaluate PRF using unconstrained PRF key and input
   * @param key unconstrained PRF key
   * @param input PRF input
   * @return PRF output
   */
  shared_ptr<vector<NativePoly>> Evaluate(
      const std::pair<vector<vector<Element>>, Matrix<Element>>& key,
      const string& input) const;

  /**
   * Method to evaluate PRF using constrained PRF key and input
   * @param constrainedKey constrained PRF key
   * @param input PRF input
   * @return PRF output
   */
  shared_ptr<vector<NativePoly>> Evaluate(
      const std::pair<Matrix<Element>, vector<vector<Matrix<Element>>>>&
          constrainedKey,
      const string& input) const;

 protected:
  /**
   * Virtual method for gamma function
   * @param m braching program matrix
   * @param s secret
   * @return matrix to encode
   */
  virtual const Matrix<Element> Gamma(const Matrix<int>& m,
                                      const Element& s) const = 0;

  usint m_base;
  usint m_chunkSize;
  usint m_length;
  usint m_adjustedLength;
  usint m_chunkExponent;
  usint m_w;
  usint m_m;

  shared_ptr<Matrix<Element>> m_J;
  shared_ptr<typename Element::Params> m_elemParams;

  typename Element::DggType m_dgg;
  typename Element::DggType m_dggLargeSigma;

  // used during encryption
  vector<NativeInteger> m_QDivtModq;

  // Stores \frac{t*{Q/q_i}^{-1}/q_i}
  std::vector<double> m_tQHatInvModqDivqFrac;

  // when log2(q_i) >= 45 bits, B = \floor[2^{\ceil{log2(q_i)/2}}
  // Stores \frac{t*{Q/q_i}^{-1}*B/q_i}
  std::vector<double> m_tQHatInvModqBDivqFrac;

  // Stores [\floor{t*{Q/q_i}^{-1}/q_i}]_t
  std::vector<NativeInteger> m_tQHatInvModqDivqModt;
  // Stores NTL precomputations for [\floor{t*{Q/q_i}^{-1}/q_i}]_t
  std::vector<NativeInteger> m_tQHatInvModqDivqModtPrecon;

  // when log2(q_i) >= 45 bits, B = \floor[2^{\ceil{log2(q_i)/2}}
  // Stores [\floor{t*{Q/q_i}^{-1}*B/q_i}]_t
  std::vector<NativeInteger> m_tQHatInvModqBDivqModt;

  // when log2 q_i >= 45 bits, B = \floor[2^{\ceil{log2(q_i)/2}}
  // Stores NTL precomputations for [\floor{t*{Q/q_i}^{-1}*B/q_i}]_t
  std::vector<NativeInteger> m_tQHatInvModqBDivqModtPrecon;

 private:
  /**
   * Method to find modulus estimate given ring dimension n
   * Used as a subroutine by constructor
   * @param n ring dimension
   * @return estimated value q of modulus
   */
  double EstimateRingModulus(usint n) const;

  /**
   * Method to create element parameters for given q and n
   * Used as a subroutine by constructor
   * @param q estimated value of modulus
   * @param n ring dimension
   * @return element parameters
   */
  shared_ptr<typename Element::Params> GenerateElemParams(double q,
                                                          usint n) const;

  /**
   * Method to encode matrix by path A_i -> A_j
   * Used as a subroutine by Constrain
   * @param trapPair trapdoor std::pair for A_i
   * @param A A_j
   * @param matrix matrix to encode
   * @return encoding of matrix
   */
  Matrix<Element> Encode(
      const std::pair<Matrix<Element>, RLWETrapdoorPair<Element>>& trapPair,
      const Matrix<Element>& A, const Matrix<Element>& matrix) const;

  /**
   * Method to transform matrix to PRF output
   * Used as a subroutine by Evaluate
   * @param matrix matrix to transform
   * @return PRF output
   */
  shared_ptr<vector<NativePoly>> TransformMatrixToPRFOutput(
      const Element& input) const;
};

template <>
shared_ptr<typename DCRTPoly::Params> BPCHCPRF<DCRTPoly>::GenerateElemParams(
    double q, usint n) const;

/**
 * @brief BP CHCPRF scheme as described in https://eprint.iacr.org/2017/143.pdf
 * @tparam Element a ring element.
 */
template <class Element>
class CC17Algorithm : public BPCHCPRF<Element> {
 public:
  /**
   * Constructor
   * @param base base for G-sampling
   * @param chunkSize number of bits encoded by one encoding matrix
   * @param length the input length
   * @param n ring dimension
   * @param w branching program width
   */
  explicit CC17Algorithm(usint base, usint chunkSize, usint length, usint n,
                         usint w);

 protected:
  /**
   * Gamma(m, s) = m * s
   * @param m braching program matrix
   * @param s secret
   * @return matrix to encode
   */
  const Matrix<Element> Gamma(const Matrix<int>& m, const Element& s) const;
};

/**
 * @brief BP CHCPRF scheme as described in https://eprint.iacr.org/2018/360.pdf
 * @tparam Element a ring element.
 */
template <class Element>
class CVW18Algorithm : public BPCHCPRF<Element> {
 public:
  /**
   * Constructor
   * @param base base for G-sampling
   * @param chunkSize number of bits encoded by one encoding matrix
   * @param length the input length
   * @param n ring dimension
   * @param v vector in branching program
   */
  explicit CVW18Algorithm(usint base, usint chunkSize, usint length, usint n,
                          const Matrix<int>& v);

 protected:
  /**
   * Gamma(m, s) = diag(s, m * s)
   * @param m braching program matrix
   * @param s secret
   * @return matrix to encode
   */
  const Matrix<Element> Gamma(const Matrix<int>& m, const Element& s) const;
};

template <class Element>
class WitnessEncryption : public BPCHCPRF<Element> {
 public:
  explicit WitnessEncryption(usint base, usint chunkSize, usint n,
                             usint numVariables, usint numClauses);

  const std::pair<Matrix<Element>, vector<vector<Matrix<Element>>>> Encrypt(
      const vector<string>& cnf, usint message) const;

  usint Decrypt(
      const std::pair<Matrix<Element>, vector<vector<Matrix<Element>>>>
          ciphertext,
      const string& x) const;

 protected:
  const Matrix<Element> Gamma(const Matrix<int>& m, const Element& s) const;
};

}  // namespace lbcrypto

#endif
