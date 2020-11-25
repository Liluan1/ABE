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

#ifndef TRAPDOOR_LIB_ABE_ABE_H_
#define TRAPDOOR_LIB_ABE_ABE_H_

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
/**
 * Setup function for Private Key Generator (PKG)
 * Digit decomposition using higher bases with balanced representation
 * Limits noise growth
 * Temporarily here; but can be made a part of Matrix<Poly> class
 *
 * templated with three Element classes E1, E2, E3)
 * @param ilParams parameter set (of type E1:Params)
 * @param base is a power of two
 * @param k bit size of modulus
 * @param &matrix to be decomposed (of elements E2)
 * @param *psi decomposed matrix (of elements E3)
 **/

/*
 * Input: base
 * Input: vector of (k+2) elements of $R_q$
 * Input: $k = \lceil \log_(base){q} \rceil$; i.e. the digit length of the
 * modulus + 1 (in base) Output: matrix of (k+2)x(k+2) elements of $R_2$ where
 * the coefficients are in balanced representation
 */

template <class Element>
void PolyVec2BalDecom(const shared_ptr<typename Element::Params> ilParams,
                      int32_t base, int k, const Matrix<Element> &pubElemB,
                      Matrix<Element> *psi) {
  usint ringDimesion = ilParams->GetCyclotomicOrder() >> 1;
  usint m = k + 2;
  BigInteger q = ilParams->GetModulus();
  auto big0 = BigInteger(0);
  auto bigBase = BigInteger(base);
  for (usint i = 0; i < m; i++)
    for (usint j = 0; j < m; j++) {
      (*psi)(j, i).SetValuesToZero();
      if ((*psi)(j, i).GetFormat() != Format::COEFFICIENT) {
        (*psi)(j, i).SwitchFormat();
      }
    }
  for (usint ii = 0; ii < m; ii++) {
    int digit_i;
    auto tB = pubElemB(0, ii);
    if (tB.GetFormat() != Format::COEFFICIENT) {
      tB.SwitchFormat();
    }

    for (usint i = 0; i < ringDimesion; i++) {
      auto coeff_i = tB.at(i);
      int j = 0;
      int flip = 0;
      while (coeff_i > big0) {
#ifdef STRANGE_MINGW_COMPILER_ERROR
        // the following line of code, with -O2 or -O3 on, crashes the mingw64
        // compiler replaced with the bracketed code below
        digit_i = coeff_i.GetDigitAtIndexForBase(1, base);
#endif
        {
          digit_i = 0;
          usint newIndex = 1;
          for (int32_t i = 1; i < base; i = i * 2) {
            digit_i += coeff_i.GetBitAtIndex(newIndex) * i;
            newIndex++;
          }
        }

        if (digit_i > (base >> 1)) {
          digit_i = base - digit_i;
          coeff_i = coeff_i + bigBase;  // math backend 2
          (*psi)(j, ii).at(i) = q - BigInteger(digit_i);
        } else if (digit_i == (base >> 1)) {
          if (flip == 0) {
            coeff_i = coeff_i + bigBase;  // math backend 2
            (*psi)(j, ii).at(i) = q - BigInteger(digit_i);
          } else {
            (*psi)(j, ii).at(i) = BigInteger(digit_i);
          }
          flip = flip ^ 1;
        } else {
          (*psi)(j, ii).at(i) = BigInteger(digit_i);
        }

        coeff_i = coeff_i.DividedBy(bigBase);
        j++;
      }
    }
  }
}

enum GaussianMode { SUBGAUSSIAN = 0, NAF = 1 };

/**
 * KPABE class definition template
 * Element is the main ring element used,
 * while Element2 is interpolated ring element.
 * e.g. DCRTPoly is element and Poly is Element2
 */
/*Element is the main ring element used, while Element2 is interpolated ring
 * element. e.g. DCRTPoly is element and Poly is Element2*/
template <class Element, class Element2 = Poly>
class KPABE {
  using ParmType = typename Element::Params;
  using IntType = typename Element::Integer;
  using DugType = typename Element::DugType;
  using DggType = typename Element::DggType;

  using Parm2Type = typename Element2::Params;
  using Int2Type = typename Element2::Integer;

 public:
  /**
   * Default Constructor
   *
   */
  explicit KPABE(GaussianMode mode = NAF);

  /**
   * Destructor for releasing dynamic memory
   * used for precomputed psi
   *
   */
  ~KPABE() {}

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
             Matrix<Element> *pubElemB);

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
   * @param *evalPubElement total number of attributes
   */
  void EvalPK(const shared_ptr<ParmType> params,
              const Matrix<Element> &pubElemB,
              Matrix<Element> *evalPubElementBf, uint32_t seed = 1);

  /**
   * Evaluation function for public vectors publicElementB
   * for the benchmark circuit
   *
   * @param params parameter set
   * @param &publicElementB is a matrix where each column corresponds to the
   * public vector of each attribute
   * @param *evalPubElement total number of attributes
   */
  void EvalPKDCRT(const shared_ptr<ParmType> params,
                  const Matrix<Element> &pubElemB,
                  Matrix<Element> *evalPubElementBf,
                  const shared_ptr<Parm2Type> ilParams, uint32_t seed = 1);

  /**
   * Evaluation function for public vectors publicElementB
   * for the benchmark circuit
   *
   * @param params parameter set
   * @param &publicElementB is a matrix where each column corresponds to the
   * public vector of each attribute
   * @param x[] array of attributes
   * @param &origCT original ciphertext
   * @param *evalAttribute evaluated value of circuit
   * @param *evalCT evaluated ciphertext value
   */
  void EvalCTDCRT(const shared_ptr<ParmType> params,
                  const Matrix<Element> &pubElemB,
                  const usint x[],                // attributes
                  const Matrix<Element> &origCT,  // original ciphtertext
                  usint *evalAttribute,           // evaluated circuit
                  Matrix<Element> *evalCT,        // evaluated ciphertext,
                  const shared_ptr<Parm2Type> ilParams, uint32_t seed = 1);

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
              const Matrix<Element> &pubElemB,
              const usint x[],                // attributes
              const Matrix<Element> &origCT,  // original ciphtertext
              usint *evalAttribute,           // evaluated circuit
              Matrix<Element> *evalCT,        // evaluated ciphertext
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
                      const Matrix<Element> &pubElemB0,
                      const Matrix<Element> &origPubElem,
                      Matrix<Element> *evalPubElem, uint32_t seed = 1);

  /**
   * Evaluation of a single NAND gate
   * NAND gate is universal,
   * any Boolean function can be constructed from NAND gates
   *
   * @param ilParams parameter set
   * @param &ctC0
   * @param x[] array of attributes
   * @param &origPubElem original matrix of public vectors for each attribute
   * @param &origCT original ciphertext
   * @param *evalAttribute evaluated value of circuit
   * @param *evalCT evaluated ciphertext value
   */
  /*
   * This is method for evaluating a single NAND gate
   */
  void NANDGateEvalCT(const shared_ptr<ParmType> ilParams,
                      const Matrix<Element> &ctC0, const usint x[],
                      const Matrix<Element> &origPubElem,
                      const Matrix<Element> &origCT, usint *evalAttribute,
                      Matrix<Element> *evalCT, uint32_t seed = 1);

  /**
   *Evaluation of simple AND Gate
   *
   * @param ilParams parameter set
   * @param &origPubElementB original matrix of public vectors for each
   *attribute
   * @param *evalPubElementBf evaluated value of public element
   */
  void ANDGateEvalPK(shared_ptr<ParmType> ilParams,
                     const Matrix<Element> &origPubElemB,
                     Matrix<Element> *evalPubElemBf, uint32_t seed = 1);
  /**
   *Evaluation of simple AND Gate
   *
   * @param ilParams parameter set
   * @param x[] array of attributes
   * @param &origPubElemB original matrix of public vectors for each attribute
   * @param &origCT original ciphertext
   * @param *evalAttribute evaluated value of circuit
   * @param *evalCT evaluated ciphertext value
   */
  void ANDGateEvalCT(const shared_ptr<ParmType> ilParams,
                     const usint x[2],  // TBA
                     const Matrix<Element> &origPubElemB,
                     const Matrix<Element> &origCT, usint *evalAttribute,
                     Matrix<Element> *evalCT, uint32_t seed = 1);

  /**
   * Encrypt Function
   *
   * @param ilParams parameter set
   * @param &pubElementA
   * @param &pubElementB
   * @param &d
   * @param x[] array of attributes
   * @param &pt
   * @param &dgg to generate error terms (Gaussian)
   * @param &dug select according to uniform distribution
   * @param &bug select according to uniform distribution binary
   * @param *ctCin resulting ciphertext Cin as per algorithm
   * @param *ctC1 c1, a separate part of the cipertext as per the algorithm
   */
  void Encrypt(const shared_ptr<ParmType> params,
               const Matrix<Element> &pubElemA, const Matrix<Element> &pubElemB,
               const Element &d,  // TBA
               const usint x[], const Element &pt, DggType &dgg, DugType &dug,
               BinaryUniformGenerator &bug, Matrix<Element> *ctCin,
               Element *ctC1);

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
              const Matrix<Element> &pubElemA, const Matrix<Element> &pubElemB,
              const Element &beta, const RLWETrapdoorPair<Element> &secElemTA,
              DggType &dgg, Matrix<Element> *sk);

  /**
   * Decrypt Function
   *
   * @param params parameter set
   * @param &sk Secret Key
   * @param &ctA ciphertext A as per paper
   * @param &evalCT evaluated ciphertext Cf pertaining to a policy
   * @param &ctC1 ciphertext C1
   * @param *dtext decrypted ciphetext
   */
  void Decrypt(const shared_ptr<ParmType> params, const Matrix<Element> &sk,
               const Matrix<Element> &ctA, const Matrix<Element> &evalCT,
               const Element &ctC1, Element *dtext);

  /**
   * Decode Function
   *
   * @param *dtext decoded ciphertext
   */
  void Decode(Poly *dtext);

  /**
   * Evaluation of a single NAND gate
   * NAND gate is universal,
   * any Boolean function can be constructed from NAND gates
   *
   * @param params parameter set
   * @param &pubElemB0
   * @param &origPubElem original matrix of public vectors for each attribute
   * @param *evalPubElem evaluated value of public element
   * @param ilParamsConsolidated consolidated params
   */
  /*
   * This is method for evaluating a single NAND gate
   */
  void NANDGateEvalPKDCRT(const shared_ptr<ParmType> params,
                          const Matrix<Element> &pubElemB0,
                          const Matrix<Element> &origPubElem,
                          Matrix<Element> *evalPubElem,
                          const shared_ptr<Parm2Type> ilParamsConsolidated,
                          uint32_t seed = 1);

  /**
   * Evaluation of a single NAND gate
   * NAND gate is universal,
   * any Boolean function can be constructed from NAND gates
   *
   * @param params parameter set
   * @param &ctC0
   * @param x[] array of attributes
   * @param &origPubElem original matrix of public vectors for each attribute
   * @param &origCT original ciphertext
   * @param *evalAttribute evaluated value of circuit
   * @param *evalCT evaluated ciphertext value
   * @param ilParamsConsolidated consolidated params
   */
  /*
   * This is method for evaluating a single NAND gate
   */
  void NANDGateEvalCTDCRT(const shared_ptr<ParmType> params,
                          const Matrix<Element> &ctC0, const usint x[],
                          const Matrix<Element> &origPubElem,
                          const Matrix<Element> &origCT, usint *evalAttribute,
                          Matrix<Element> *evalCT,
                          const shared_ptr<Parm2Type> ilParamsConsolidated,
                          uint32_t seed = 1);

  /**
   *Evaluation of simple Public key AND Gate DCRT
   *
   * @param params parameter set
   * @param &origPubElementB original matrix of public vectors for each
   *attribute
   * @param *evalPubElementBf evaluated value of public element
   * @param ilParamsConsolidated consolidated params
   */
  void ANDGateEvalPKDCRT(const shared_ptr<ParmType> params,
                         const Matrix<Element> &origPubElemB,
                         Matrix<Element> *evalPubElemBf,
                         const shared_ptr<Parm2Type> ilParamsConsolidated,
                         uint32_t seed = 1);
  /**
   *Evaluation of simple Ciphertext AND Gate
   *
   * @param params parameter set
   * @param x[] array of attributes
   * @param &origPubElemB original matrix of public vectors for each attribute
   * @param &origCT original ciphertext
   * @param *evalAttribute evaluated value of circuit
   * @param *evalCT evaluated ciphertext value
   * @param ilParamsConsolidated consolidated params
   */
  void ANDGateEvalCTDCRT(const shared_ptr<ParmType> params,
                         const usint x[2],  // TBA
                         const Matrix<Element> &origPubElemB,
                         const Matrix<Element> &origCT, usint *evalAttribute,
                         Matrix<Element> *evalCT,
                         const shared_ptr<Parm2Type> ilParamsConsolidated,
                         uint32_t seed = 1);

 private:
  usint m_k;       // number of bits of the modulus
  usint m_ell;     // number of attributes
  usint m_N;       // ring dimension
  BigInteger m_q;  // modulus
  usint m_m;       // m = k+2
  int32_t m_base;  // base, a power of two

  GaussianMode m_mode;
  LatticeSubgaussianUtility<Int2Type> m_util;
};

}  // namespace lbcrypto

#endif /* TRAPDOOR_LIB_ABE_ABE_H_ */
