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
//
// @section DESCRIPTION
//
// This code provides functionality for Key-policy attribute-based encryption
// (KP-ABE). The algorithms and naming conventions can be found from
// this paper: https://eprint.iacr.org/2017/601.pdf

#ifndef TRAPDOOR_LIB_ABE_ABE_RNS_H_
#define TRAPDOOR_LIB_ABE_ABE_RNS_H_

#include <cmath>
#include <memory>
#include <vector>

#include "cryptocontexthelper.h"
#include "lattice/elemparams.h"
#include "lattice/ildcrtparams.h"
#include "lattice/ilelement.h"
#include "lattice/ilparams.h"
#include "lattice/trapdoor.h"
#include "math/backend.h"
#include "math/distrgen.h"
#include "palisade.h"
#include "subgaussian/subgaussian.h"
#include "utils/inttypes.h"
/**
 * @namespace lbcrypto
 * The namespace of lbcrypto
 */
namespace lbcrypto {

class KPABErns {
  using ParmType = typename DCRTPoly::Params;
  using BugType = typename DCRTPoly::BugType;
  using DugType = typename DCRTPoly::DugType;
  using DggType = typename DCRTPoly::DggType;

 public:
  /**
   * Default Constructor
   *
   */
  KPABErns() : m_k(0), m_ell(0), m_N(0), m_q(0), m_m(0), m_base(0) {}

  /**
   * Destructor for releasing dynamic memory
   * used for precomputed psi
   *
   */
  ~KPABErns() {}

  /**
   * Setup function for Private Key Generator (PKG)
   *
   * @param ilParams parameter set
   * @param base is a power of two
   * @param ell total number of attributes
   * @param &dug
   * @param *publicElementB is a matrix where each column corresponds to the
   * public vector of each attribute
   */
  void Setup(const shared_ptr<ParmType> params, int32_t base,
             usint ell,     // number of attributes
             DugType &dug,  // select according to uniform distribution
             Matrix<DCRTPoly> *pubElemB);

  /**
   * Setup function for all parties except the Private Key Generator (PKG)
   *
   * @param ilParams parameter set
   * @param base is a power of two
   * @param ell total number of attributes
   */
  void Setup(const shared_ptr<ParmType> params, int32_t base, const usint ell);

  /**
   * Evaluation function for public vectors publicElementB
   * for the benchmark circuit
   *
   * @param ilParams parameter set
   * @param &publicElementB is a matrix where each column corresponds to the
   * public vector of each attribute
   * @param *evalPubDCRTPoly total number of attributes
   */
  void EvalPK(const shared_ptr<ParmType> params,
              const Matrix<DCRTPoly> &pubElemB,
              Matrix<DCRTPoly> *evalPubElementBf, uint32_t seed = 1);

  /**
   * Evaluation function for public vectors publicElementB
   * for the benchmark circuit
   *
   * @param ilParams parameter set
   * @param &publicElementB is a matrix where each column corresponds to the
   * public vector of each attribute
   * @param x[] array of attributes
   * @param &origCT original ciphertext
   * @param *evalAttribute evaluated value of circuit
   * @param *evalCT evaluated ciphertext value
   */
  void EvalCT(const shared_ptr<ParmType> ilParams,
              const Matrix<DCRTPoly> &pubElemB,
              const usint x[],                 // attributes
              const Matrix<DCRTPoly> &origCT,  // original ciphtertext
              usint *evalAttribute,            // evaluated circuit
              Matrix<DCRTPoly> *evalCT,        // evaluated ciphertext
              uint32_t seed = 1);

  /**
   * Evaluation of a single NAND gate
   * NAND gate is universal,
   * any Boolean function can be constructed from NAND gates
   *
   * @param ilParams parameter set
   * @param &pubElemB0
   * @param &origPubElem original matrix of public vectors for each attribute
   * @param *evalPubElem evaluated value of public element
   */
  /*
   * This is method for evaluating a single NAND gate
   */
  void NANDGateEvalPK(const shared_ptr<ParmType> ilParams,
                      const Matrix<DCRTPoly> &pubElemB,
                      Matrix<DCRTPoly> *evalPubElem, uint32_t seed = 1);

  /**
   * Evaluation of a single NAND gate
   * NAND gate is universal,
   * any Boolean function can be constructed from NAND gates
   *
   * @param ilParams parameter set
   * @param &origPubElem original matrix of public vectors for each attribute
   * @param x[] array of attributes
   * @param &origCT original ciphertext
   * @param *evalAttribute evaluated value of circuit
   * @param *evalCT evaluated ciphertext value
   */
  /*
   * This is method for evaluating a single NAND gate
   */
  void NANDGateEvalCT(const shared_ptr<ParmType> ilParams,
                      const Matrix<DCRTPoly> &origPubElem, const usint x[],
                      const Matrix<DCRTPoly> &origCT, usint *evalAttribute,
                      Matrix<DCRTPoly> *evalCT, uint32_t seed = 1);

  /**
   *Evaluation of simple AND Gate
   *
   * @param ilParams parameter set
   * @param &origPubElementB original matrix of public vectors for each
   *attribute
   * @param *evalPubElementBf evaluated value of public element
   */
  void ANDGateEvalPK(shared_ptr<ParmType> ilParams,
                     const Matrix<DCRTPoly> &origPubElemB,
                     Matrix<DCRTPoly> *evalPubElemBf, uint32_t seed = 1);

  /**
   *Evaluation of simple AND Gate
   *
   * @param ilParams parameter set
   * @param &origPubElemB original matrix of public vectors for each
   *attribute
   * @param x[] array of attributes
   * @param &origCT original ciphertext
   * @param *evalAttribute evaluated value of circuit
   * @param *evalCT evaluated ciphertext value
   */
  void ANDGateEvalCT(const shared_ptr<ParmType> ilParams,
                     const Matrix<DCRTPoly> &origPubElemB,
                     const usint x[2],  // TBA
                     const Matrix<DCRTPoly> &origCT, usint *evalAttribute,
                     Matrix<DCRTPoly> *evalCT, uint32_t seed = 1);

  /**
   * Encrypt Function
   *
   * @param ilParams parameter set
   * @param &pubElementA
   * @param &pubElementB
   * @param &d TODO
   * @param x[] array of attributes
   * @param &pt
   * @param &dgg to generate error terms (Gaussian)
   * @param &dug select according to uniform distribution
   * @param &bug select according to uniform distribution binary
   * @param *ctCin resulting ciphertext Cin as per algorithm
   * @param *ctC1 c1, a separate part of the cipertext as per the algorithm
   */
  void Encrypt(const shared_ptr<ParmType> params,
               const Matrix<DCRTPoly> &pubElemA,
               const Matrix<DCRTPoly> &pubElemB, const DCRTPoly &d,
               const usint x[], const NativePoly &pt, DggType &dgg,
               DugType &dug, BugType &bug, Matrix<DCRTPoly> *ctCin,
               DCRTPoly *ctC1);

  /**
   * KeyGen Function
   *
   * @param params parameter set
   * @param &pubElementA Public parameter $A \in R_q^{1 \times w}$
   * @param &pubElementB Public parameter $B \in R_q^{ell \times k}$
   * @param &beta public key $d \in R_q$  TBA
   * @param &secElemTA Secret parameter $T_H \in R_q^{1 \times k} \times
   * R_q^{1 \times k}$
   * @param &dgg to generate error terms (Gaussian)
   * @param *sk secret key
   */
  void KeyGen(const shared_ptr<ParmType> params,
              const Matrix<DCRTPoly> &pubElemA,
              const Matrix<DCRTPoly> &pubElemB, const DCRTPoly &beta,
              const RLWETrapdoorPair<DCRTPoly> &secElemTA, DggType &dgg,
              Matrix<DCRTPoly> *sk);

  /**
   * Decrypt Function
   *
   * @param params parameter set
   * @param &sk Secret Key
   * @param &ctA ciphertext A as per paper
   * @param &evalCT evaluated ciphertext Cf pertaining to a policy
   * @param &ctC1 ciphertext C1
   * @param *ptext decrypted ciphetext
   */
  void Decrypt(const shared_ptr<ParmType> params, const Matrix<DCRTPoly> &sk,
               const Matrix<DCRTPoly> &ctA, const Matrix<DCRTPoly> &evalCT,
               const DCRTPoly &ctC1, NativePoly *ptext);

 private:
  usint m_k;       // number of bits in the modulus
  usint m_ell;     // number of attributes
  usint m_N;       // ring dimension
  BigInteger m_q;  // modulus
  usint m_m;       // m = k+2
  int32_t m_base;  // base, a power of two

  vector<LatticeSubgaussianUtility<NativeInteger>> m_util;

  // Stores [\floor{Q/t}]_{q_i}
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
};

}  // namespace lbcrypto

#endif /* TRAPDOOR_LIB_ABE_ABE_RNS_H_ */
