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
// This code provides functionality for KP_ABE. The algorithms and naming
// conventions can be found from this paper:
// https://eprint.iacr.org/2017/601.pdf

#include "abe/kp_abe.h"

namespace lbcrypto {

template <class Element, class Element2>
KPABE<Element, Element2>::KPABE(GaussianMode mode)
    : m_k(0), m_ell(0), m_N(0), m_q(0), m_m(0), m_base(0), m_mode(mode) {}

/*
 * This is a setup function for Private Key Generator (PKG);
 * generates master public key (MPK) and master secret key
 * m_ell is the number of attributes
 */
template <class Element, class Element2>
void KPABE<Element, Element2>::Setup(
    const shared_ptr<ParmType> params, int32_t base,
    usint ell,     // number of attributes
    DugType &dug,  // select according to uniform distribution
    Matrix<Element> *pubElemB) {
  m_N = params->GetCyclotomicOrder() >> 1;
  BigInteger q(params->GetModulus());
  m_q = q;
  m_base = base;

  double val = q.ConvertToDouble();
  double logTwo = log(val - 1.0) / log(base) + 1.0;

  if (m_mode == NAF)
    m_k = (usint)floor(logTwo) + 1;
  else
    m_k = (long)ceil(log2(val) / log2(base));

  m_m = m_k + 2;

  m_ell = ell;

  m_util = LatticeSubgaussianUtility<Int2Type>(m_base, q, m_k);

  for (usint i = 0; i < (*pubElemB).GetRows(); i++)
    for (usint j = 0; j < (*pubElemB).GetCols(); j++) {
      if ((*pubElemB)(i, j).GetFormat() != Format::COEFFICIENT)
        (*pubElemB)(i, j).SwitchFormat();
      (*pubElemB)(i, j) = Element(dug, params, Format::COEFFICIENT);
      (*pubElemB)(i, j)
          .SwitchFormat();  // always kept in Format::EVALUATION format
    }
}

/**
 * This setup function is used by users; namely senders and receivers
 * Initialize private members of the object such as modulus, cyclotomic order,
 * etc. m_ell is the number of attributes
 */
template <class Element, class Element2>
void KPABE<Element, Element2>::Setup(const shared_ptr<ParmType> params,
                                     int32_t base, const usint ell) {
  m_N = params->GetCyclotomicOrder() >> 1;
  BigInteger q(params->GetModulus());
  m_q = q;
  m_base = base;

  double val = q.ConvertToDouble();
  double logTwo = log(val - 1.0) / log(base) + 1.0;

  if (m_mode == NAF)
    m_k = (usint)floor(logTwo) + 1;
  else
    m_k = (long)ceil(log2(val) / log2(base));

  m_m = m_k + 2;

  m_ell = ell;

  m_util = LatticeSubgaussianUtility<Int2Type>(m_base, q, m_k);
}

/*
 * Given public parameters, attribute values and ciphertexts corresponding to
 * attributes, computes the ciphertext and the public key evalPubElement for the
 * circuit of attributes m_ell is the number of attributes and the circuit is
 * assumed to be a binary tree of NAND gates Thus, m_ell must be a power of two
 */
template <class Element, class Element2>
void KPABE<Element, Element2>::EvalPK(const shared_ptr<ParmType> params,
                                      const Matrix<Element> &pubElemB,
                                      Matrix<Element> *evalPubElemBf,
                                      uint32_t seed) {
  auto zero_alloc = Element::Allocator(params, Format::EVALUATION);

  usint gateCnt = m_ell - 1;

  Matrix<Element> psi(zero_alloc, m_m,
                      m_m);  // Needed for bit decomposition matrices
  // w stands for wire
  Matrix<Element> wpublicElementB(
      zero_alloc, gateCnt,
      m_m);  // Bis associated with internal wires of the circuit
  // Temporary variables for bit decomposition operation
  Matrix<Element> negpublicElementB(zero_alloc, 1,
                                    m_m);  // Format::EVALUATION (NTT domain)
  std::vector<Element> digitsC1(m_m);

  // Input level of the circuit
  usint t = m_ell >> 1;  // the number of the gates in the first level (the
                         // number of input gates)
  // looping to evaluate and calculate w, wB, wC and R for all first level input
  // gates
  for (usint i = 0; i < t; i++) {
    for (usint j = 0; j < m_m; j++)  // Negating Bis for bit decomposition
      negpublicElementB(0, j) = pubElemB(2 * i + 1, j).Negate();

    if (m_mode == NAF)
      PolyVec2BalDecom<Element>(params, m_base, m_k, negpublicElementB, &psi);
    else
      InverseRingVector<Element>(m_util, params, negpublicElementB, seed, &psi);

    psi.SwitchFormat();

    /* Psi^T*C2 and B2*Psi */
    for (usint j = 0; j < m_m;
         j++) {  // the following two for loops are for vector matrix
                 // multiplication (a.k.a B(i+1) * BitDecompose(-Bi) and  gamma
                 // (0, 2) (for the second attribute of the circuit) *
                 // bitDecompose(-B))
      wpublicElementB(i, j) =
          pubElemB(2 * i + 2, 0) * psi(0, j);  // B2 * BD(-Bi)
      for (usint k = 1; k < m_m; k++) {
        wpublicElementB(i, j) += pubElemB(2 * i + 2, k) * psi(k, j);
      }
    }

    for (usint j = 0; j < m_m; j++) {
      wpublicElementB(i, j) = pubElemB(0, j) - wpublicElementB(i, j);
    }
  }

  /* For internal wires of the circuit.
   * Depth 0 refers to the circuit level where the input gates are located.
   * Thus, we start with depth 1
   */
  usint depth = log2(m_ell);
  for (usint d = 1; d < depth; d++) {
    usint inStart =
        m_ell -
        (m_ell >> (d - 1));  // Starting index for the input wires in level d
    usint outStart =
        m_ell - (m_ell >> d);  // Starting index for the output wires in level d
    usint gCntinLeveld = m_ell >> (d + 1);  // number of gates in level d

    for (usint i = 0; i < gCntinLeveld; i++) {
      for (usint j = 0; j < m_m; j++)
        negpublicElementB(0, j) = wpublicElementB(inStart + 2 * i, j).Negate();

      if (m_mode == NAF)
        PolyVec2BalDecom<Element>(params, m_base, m_k, negpublicElementB, &psi);
      else
        InverseRingVector<Element>(m_util, params, negpublicElementB, seed,
                                   &psi);

      psi.SwitchFormat();

      for (usint j = 0; j < m_m; j++) {
        wpublicElementB(outStart + i, j) =
            wpublicElementB(inStart + 2 * i + 1, 0) * psi(0, j);  // B2 * Psi
        for (usint k = 1; k < m_m; k++) {
          wpublicElementB(outStart + i, j) +=
              wpublicElementB(inStart + 2 * i + 1, k) * psi(k, j);  // B2 * Psi
        }
      }

      for (usint j = 0; j < m_m; j++) {
        wpublicElementB(outStart + i, j) =
            pubElemB(0, j) - wpublicElementB(outStart + i, j);
      }
    }
  }

  for (usint j = 0; j < m_m; j++) {
    (*evalPubElemBf)(0, j) = wpublicElementB(gateCnt - 1, j);
  }
}

/*
 * Given public parameters, attribute values and ciphertexts corresponding to
 * attributes, computes the ciphertext and the public key evalPubElement for the
 * circuit of attributes m_ell is the number of attributes and the circuit is
 * assumed to be a binary tree of NAND gates Thus, m_ell must be a power of two
 */
template <class Element, class Element2>
void KPABE<Element, Element2>::EvalPKDCRT(
    const shared_ptr<ParmType> params, const Matrix<Element> &pubElemB,
    Matrix<Element> *evalPubElemBf,
    const shared_ptr<Parm2Type> ilParamsConsolidated, uint32_t seed) {
  auto zero_alloc = Element::Allocator(params, Format::EVALUATION);

  usint gateCnt = m_ell - 1;

  // Needed for bit decomposition matrices
  Matrix<Element> psi(zero_alloc, m_m, m_m);
  // w stands for wire
  Matrix<Element> wpublicElementB(
      zero_alloc, gateCnt,
      m_m);  // Bis associated with internal wires of the circuit
  // Temporary variables for bit decomposition operation
  Matrix<Element> negPubElemB(zero_alloc, 1,
                              m_m);  // Format::EVALUATION (NTT domain)
  std::vector<Element> digitsC1(m_m);

  // Added for dcrt
  auto zero_alloc_el2 =
      Element2::Allocator(ilParamsConsolidated, Format::COEFFICIENT);

  // Needed for Bit Decomposition (BD) matrices
  Matrix<Element2> psiPoly(zero_alloc_el2, m_m, m_m);

  Matrix<Element2> negBPolyMatrix(zero_alloc_el2, 1, m_m);

  // Input level of the circuit
  usint t = m_ell >> 1;  // the number of the gates in the first level (the
                         // number of input gates)

  // looping to evaluate and calculate w, wB, wC
  // and R for all first level input gates
  for (usint i = 0; i < t; i++) {
    for (usint j = 0; j < m_m; j++)  // Negating Bis for bit decomposition
      negPubElemB(0, j) = pubElemB(2 * i + 1, j).Negate();

    // DCRT INTERPOLATE

    for (usint i = 0; i < negPubElemB.GetRows(); i++) {
      for (usint j = 0; j < negPubElemB.GetCols(); j++) {
        negBPolyMatrix(i, j) = negPubElemB(i, j).CRTInterpolate();
      }
    }

    if (m_mode == NAF)
      PolyVec2BalDecom<Element2>(ilParamsConsolidated, m_base, m_k,
                                 negBPolyMatrix, &psiPoly);
    else
      InverseRingVector<Element2>(m_util, ilParamsConsolidated, negBPolyMatrix,
                                  seed, &psiPoly);

    // DCRT CREATE

    for (usint i = 0; i < psiPoly.GetRows(); i++) {
      for (usint j = 0; j < psiPoly.GetCols(); j++) {
        Element temp(psiPoly(i, j), params);
        psi(i, j) = temp;
      }
    }

    psi.SwitchFormat();

    /* Psi^T*C2 and B2*Psi */
    for (usint j = 0; j < m_m;
         j++) {  // the following two for loops are for vector matrix
                 // multiplication (a.k.a B(i+1) * BitDecompose(-Bi) and  gamma
                 // (0, 2) (for the second attribute of the circuit) *
                 // bitDecompose(-B))
      wpublicElementB(i, j) =
          pubElemB(2 * i + 2, 0) * psi(0, j);  // B2 * BD(-Bi)
      for (usint k = 1; k < m_m; k++) {
        wpublicElementB(i, j) += pubElemB(2 * i + 2, k) * psi(k, j);
      }
    }

    for (usint j = 0; j < m_m; j++) {
      wpublicElementB(i, j) = pubElemB(0, j) - wpublicElementB(i, j);
    }
  }

  /* For internal wires of the circuit.
   * Depth 0 refers to the circuit level where the input gates are located.
   * Thus, we start with depth 1
   */

  Matrix<Element2> negBPolyMatrix_Two(zero_alloc_el2, 1, m_m);
  // Needed for Bit Decomposition (BD) matrices
  Matrix<Element2> psiPoly_Two(zero_alloc_el2, m_m, m_m);

  usint depth = log2(m_ell);
  for (usint d = 1; d < depth; d++) {
    // Starting index for the input wires in level d
    usint inStart = m_ell - (m_ell >> (d - 1));
    // Starting index for the output wires in level d
    usint outStart = m_ell - (m_ell >> d);
    usint gCntinLeveld = m_ell >> (d + 1);  // number of gates in level d

    for (usint i = 0; i < gCntinLeveld; i++) {
      for (usint j = 0; j < m_m; j++)
        negPubElemB(0, j) = wpublicElementB(inStart + 2 * i, j).Negate();

      // DCRT INTERPOLATE

      for (usint i = 0; i < negBPolyMatrix_Two.GetRows(); i++) {
        for (usint j = 0; j < negBPolyMatrix_Two.GetCols(); j++) {
          //    std::cout << "i:" << i << "j:" << j << std::endl;
          negBPolyMatrix_Two(i, j) = negPubElemB(i, j).CRTInterpolate();
        }
      }

      if (m_mode == NAF)
        PolyVec2BalDecom<Element2>(ilParamsConsolidated, m_base, m_k,
                                   negBPolyMatrix_Two, &psiPoly_Two);
      else
        InverseRingVector<Element2>(m_util, ilParamsConsolidated,
                                    negBPolyMatrix_Two, seed, &psiPoly_Two);

      // DCRT CREATE

      for (usint i = 0; i < psiPoly_Two.GetRows(); i++) {
        for (usint j = 0; j < psiPoly_Two.GetCols(); j++) {
          Element temp(psiPoly_Two(i, j), params);
          psi(i, j) = temp;
        }
      }

      psi.SwitchFormat();

      for (usint j = 0; j < m_m; j++) {
        wpublicElementB(outStart + i, j) =
            wpublicElementB(inStart + 2 * i + 1, 0) * psi(0, j);  // B2 * Psi
        for (usint k = 1; k < m_m; k++) {
          wpublicElementB(outStart + i, j) +=
              wpublicElementB(inStart + 2 * i + 1, k) * psi(k, j);  // B2 * Psi
        }
      }

      for (usint j = 0; j < m_m; j++) {
        wpublicElementB(outStart + i, j) =
            pubElemB(0, j) - wpublicElementB(outStart + i, j);
      }
    }
  }

  for (usint j = 0; j < m_m; j++) {
    (*evalPubElemBf)(0, j) = wpublicElementB(gateCnt - 1, j);
  }
}

/*
 * Given public parameters, attribute values and ciphertexts corresponding to
 * attributes, computes the ciphertext and the public key Bf for the circuit of
 * attributes m_ell is the number of attributes and the circuit is assumed to be
 * a binary tree of NAND gates Thus, m_ell must be a power of two
 */
template <class Element, class Element2>
void KPABE<Element, Element2>::EvalCT(const shared_ptr<ParmType> ilParams,
                                      const Matrix<Element> &pubElemB,
                                      const usint x[],  // Attributes
                                      const Matrix<Element> &origCT,
                                      usint *evalAttributes,
                                      Matrix<Element> *evalCT, uint32_t seed) {
  // Part pertaining to A (does not change)
  for (usint i = 0; i < m_m; i++) (*evalCT)(0, i) = origCT(0, i);

  auto zero_alloc = Element::Allocator(ilParams, Format::EVALUATION);

  usint gateCnt = m_ell - 1;
  Matrix<Element> psi(zero_alloc, m_m,
                      m_m);  // Needed for Bit Decomposition (BD) matrices

  // w stands for Wire
  Matrix<Element> wPublicElementB(
      zero_alloc, gateCnt,
      m_m);  // Bis associated with internal wires of the circuit
  Matrix<Element> wCT(
      zero_alloc, gateCnt,
      m_m);  // Ciphertexts associated with internal wires of the circuit

  // Attribute values associated with internal wires of the circuit
  std::vector<usint> wX(gateCnt);

  // Temporary variables for bit decomposition operation
  Matrix<Element> negB(zero_alloc, 1, m_m);  // Format::EVALUATION (NTT domain)
  std::vector<Poly> digitsC1(m_m);

  // Input level of the circuit
  usint t = m_ell >> 1;  // the number of the gates in the first level (the
                         // number of input gates)

  // looping to evaluate and calculate w, wB, wC
  // and R for all first level input gates
  for (usint i = 0; i < t; i++) {
    wX[i] =
        x[0] - x[2 * i + 1] * x[2 * i + 2];  // calculating binary wire value

    for (usint j = 0; j < m_m; j++)  // Negating Bis for bit decomposition
      negB(0, j) = pubElemB(2 * i + 1, j).Negate();

    if (m_mode == NAF)
      PolyVec2BalDecom<Element>(ilParams, m_base, m_k, negB, &psi);
    else
      InverseRingVector<Element>(m_util, ilParams, negB, seed, &psi);

    psi.SwitchFormat();

    /*Starting computation for a NAND circuit*/
    /* x2 * C1 */
    for (usint j = 0; j < m_m; j++) {
      if (x[2 * i + 2] != 0)
        wCT(i, j) = origCT(2 * i + 1, j);
      else
        wCT(i, j).SetValuesToZero();
    }

    /* Psi^T*C2 and B2*Psi */
    for (usint j = 0; j < m_m;
         j++) {  // the following two for loops are for vector matrix
                 // multiplication (a.k.a B(i+1) * BitDecompose(-Bi) and  gamma
                 // (0, 2) (for the second attribute of the circuit) *
                 // bitDecompose(-B))
      wPublicElementB(i, j) =
          pubElemB(2 * i + 2, 0) * psi(0, j);         // B2 * BD(-Bi)
      wCT(i, j) += psi(0, j) * origCT(2 * i + 2, 0);  // BD(-Bi)*C2
      for (usint k = 1; k < m_m; k++) {
        wPublicElementB(i, j) += pubElemB(2 * i + 2, k) * psi(k, j);
        wCT(i, j) += psi(k, j) * origCT(2 * i + 2, k);
      }
    }

    /* B0 - B2*R and C0 - x2*C1 - C2*R */
    for (usint j = 0; j < m_m; j++) {
      wPublicElementB(i, j) = pubElemB(0, j) - wPublicElementB(i, j);
      wCT(i, j) = origCT(0, j) - wCT(i, j);  // C0 - x2*C1 - R*C2
    }
  }

  /* For internal wires of the circuit.
   * Depth 0 refers to the circuit level where the input gates are located.
   * Thus, we start with depth 1
   */
  usint depth = log2(m_ell);
  for (usint d = 1; d < depth; d++) {
    // Starting index for the input wires in level d
    usint InStart = m_ell - (m_ell >> (d - 1));
    // Starting index for the output wires in level d
    usint OutStart = m_ell - (m_ell >> d);
    usint gCntinLeveld = m_ell >> (d + 1);  // number of gates in level d

    for (usint i = 0; i < gCntinLeveld; i++) {
      wX[OutStart + i] = x[0] - wX[InStart + 2 * i] * wX[InStart + 2 * i + 1];

      for (usint j = 0; j < m_m; j++)
        negB(0, j) = wPublicElementB(InStart + 2 * i, j).Negate();

      if (m_mode == NAF)
        PolyVec2BalDecom<Element>(ilParams, m_base, m_k, negB, &psi);
      else
        InverseRingVector<Element>(m_util, ilParams, negB, seed, &psi);

      psi.SwitchFormat();

      // x2*C1
      for (usint j = 0; j < m_m; j++) {
        if (wX[InStart + 2 * i + 1] != 0)
          wCT(OutStart + i, j) = wCT(InStart + 2 * i, j);
        else
          wCT(OutStart + i, j).SetValuesToZero();
      }

      for (usint j = 0; j < m_m; j++) {
        wPublicElementB(OutStart + i, j) =
            wPublicElementB(InStart + 2 * i + 1, 0) * psi(0, j);  // B2 * psi
        wCT(OutStart + i, j) +=
            psi(0, j) * wCT(InStart + 2 * i + 1, 0);  // psi * C2
        for (usint k = 1; k < m_m; k++) {
          wPublicElementB(OutStart + i, j) +=
              wPublicElementB(InStart + 2 * i + 1, k) * psi(k, j);  // B2 * psi
          wCT(OutStart + i, j) +=
              psi(k, j) * wCT(InStart + 2 * i + 1, k);  // psi * C2
        }
      }

      for (usint j = 0; j < m_m; j++) {
        wPublicElementB(OutStart + i, j) =
            pubElemB(0, j) - wPublicElementB(OutStart + i, j);
        wCT(OutStart + i, j) = origCT(0, j) - wCT(OutStart + i, j);
      }
    }
  }

  for (usint j = 0; j < m_m; j++) {
    (*evalCT)(0, j) = wCT(gateCnt - 1, j);
  }

  (*evalAttributes) = wX[gateCnt - 1];
}

/*
 * Given public parameters, attribute values and ciphertexts corresponding to
 * attributes, computes the ciphertext and the public key Bf for the circuit of
 * attributes m_ell is the number of attributes and the circuit is assumed to be
 * a binary tree of NAND gates Thus, m_ell must be a power of two
 */
template <class Element, class Element2>
void KPABE<Element, Element2>::EvalCTDCRT(
    const shared_ptr<ParmType> params, const Matrix<Element> &pubElemB,
    const usint x[],  // Attributes
    const Matrix<Element> &origCT, usint *evalAttributes,
    Matrix<Element> *evalCT, const shared_ptr<Parm2Type> ilParamsConsolidated,
    uint32_t seed) {
  // Part pertaining to A (does not change)
  for (usint i = 0; i < m_m; i++) (*evalCT)(0, i) = origCT(0, i);

  auto zero_alloc = Element::Allocator(params, Format::EVALUATION);

  usint gateCnt = m_ell - 1;
  Matrix<Element> psi(zero_alloc, m_m, m_m);
  // w stands for Wire
  Matrix<Element> wPublicElementB(
      zero_alloc, gateCnt,
      m_m);  // Bis associated with internal wires of the circuit
  Matrix<Element> wCT(
      zero_alloc, gateCnt,
      m_m);  // Ciphertexts associated with internal wires of the circuit

  // Attribute values associated with internal wires of the circuit
  std::vector<usint> wX(gateCnt);

  // Temporary variables for bit decomposition operation
  Matrix<Element> negB(zero_alloc, 1, m_m);  // Format::EVALUATION (NTT domain)
  std::vector<Poly> digitsC1(m_m);

  // Input level of the circuit
  usint t = m_ell >> 1;  // the number of the gates in the first level (the
                         // number of input gates)

  auto zero_alloc_el2 =
      Element2::Allocator(ilParamsConsolidated, Format::COEFFICIENT);

  Matrix<Element2> psiPoly(zero_alloc_el2, m_m,
                           m_m);  // Needed for Bit Decomposition (BD) matrices

  Matrix<Element2> negBPolyMatrix(zero_alloc_el2, 1, m_m);

  // looping to evaluate and calculate w, wB, wC
  // and R for all first level input gates
  for (usint i = 0; i < t; i++) {
    // calculating binary wire value
    wX[i] = x[0] - x[2 * i + 1] * x[2 * i + 2];

    for (usint j = 0; j < m_m; j++)  // Negating Bis for bit decomposition
      negB(0, j) = pubElemB(2 * i + 1, j).Negate();

    for (usint i = 0; i < negB.GetRows(); i++) {
      for (usint j = 0; j < negB.GetCols(); j++) {
        negBPolyMatrix(i, j) = negB(i, j).CRTInterpolate();
      }
    }

    if (m_mode == NAF)
      PolyVec2BalDecom<Element2>(ilParamsConsolidated, m_base, m_k,
                                 negBPolyMatrix, &psiPoly);
    else
      InverseRingVector<Element2>(m_util, ilParamsConsolidated, negBPolyMatrix,
                                  seed, &psiPoly);

    // DCRT CREATE

    for (usint i = 0; i < psiPoly.GetRows(); i++) {
      for (usint j = 0; j < psiPoly.GetCols(); j++) {
        Element temp(psiPoly(i, j), params);
        psi(i, j) = temp;
      }
    }

    psi.SwitchFormat();
    /*Starting computation for a NAND circuit*/
    /* x2 * C1 */
    for (usint j = 0; j < m_m; j++) {
      if (x[2 * i + 2] != 0)
        wCT(i, j) = origCT(2 * i + 1, j);
      else
        wCT(i, j).SetValuesToZero();
    }

    /* Psi^T*C2 and B2*Psi */
    for (usint j = 0; j < m_m;
         j++) {  // the following two for loops are for vector matrix
                 // multiplication (a.k.a B(i+1) * BitDecompose(-Bi) and  gamma
                 // (0, 2) (for the second attribute of the circuit) *
                 // bitDecompose(-B))
      wPublicElementB(i, j) =
          pubElemB(2 * i + 2, 0) * psi(0, j);         // B2 * BD(-Bi)
      wCT(i, j) += psi(0, j) * origCT(2 * i + 2, 0);  // BD(-Bi)*C2
      for (usint k = 1; k < m_m; k++) {
        wPublicElementB(i, j) += pubElemB(2 * i + 2, k) * psi(k, j);
        wCT(i, j) += psi(k, j) * origCT(2 * i + 2, k);
      }
    }

    /* B0 - B2*R and C0 - x2*C1 - C2*R */
    for (usint j = 0; j < m_m; j++) {
      wPublicElementB(i, j) = pubElemB(0, j) - wPublicElementB(i, j);
      wCT(i, j) = origCT(0, j) - wCT(i, j);  // C0 - x2*C1 - R*C2
    }
  }

  /* For internal wires of the circuit.
   * Depth 0 refers to the circuit level where the input gates are located.
   * Thus, we start with depth 1
   */

  Matrix<Element2> psiPoly_two(
      zero_alloc_el2, m_m, m_m);  // Needed for Bit Decomposition (BD) matrices

  Matrix<Element2> negBPolyMatrix_two(zero_alloc_el2, 1, m_m);

  usint depth = log2(m_ell);
  for (usint d = 1; d < depth; d++) {
    usint InStart =
        m_ell -
        (m_ell >> (d - 1));  // Starting index for the input wires in level d
    usint OutStart =
        m_ell - (m_ell >> d);  // Starting index for the output wires in level d
    usint gCntinLeveld = m_ell >> (d + 1);  // number of gates in level d

    for (usint i = 0; i < gCntinLeveld; i++) {
      wX[OutStart + i] = x[0] - wX[InStart + 2 * i] * wX[InStart + 2 * i + 1];

      for (usint j = 0; j < m_m; j++)
        negB(0, j) = wPublicElementB(InStart + 2 * i, j).Negate();

      for (usint i = 0; i < negB.GetRows(); i++) {
        for (usint j = 0; j < negB.GetCols(); j++) {
          negBPolyMatrix_two(i, j) = negB(i, j).CRTInterpolate();
        }
      }

      if (m_mode == NAF)
        PolyVec2BalDecom<Element2>(ilParamsConsolidated, m_base, m_k,
                                   negBPolyMatrix_two, &psiPoly_two);
      else
        InverseRingVector<Element2>(m_util, ilParamsConsolidated,
                                    negBPolyMatrix_two, seed, &psiPoly_two);

      for (usint i = 0; i < psiPoly_two.GetRows(); i++) {
        for (usint j = 0; j < psiPoly_two.GetCols(); j++) {
          Element temp(psiPoly_two(i, j), params);
          psi(i, j) = temp;
        }
      }

      psi.SwitchFormat();

      // x2*C1
      for (usint j = 0; j < m_m; j++) {
        if (wX[InStart + 2 * i + 1] != 0)
          wCT(OutStart + i, j) = wCT(InStart + 2 * i, j);
        else
          wCT(OutStart + i, j).SetValuesToZero();
      }

      for (usint j = 0; j < m_m; j++) {
        wPublicElementB(OutStart + i, j) =
            wPublicElementB(InStart + 2 * i + 1, 0) * psi(0, j);  // B2 * psi
        wCT(OutStart + i, j) +=
            psi(0, j) * wCT(InStart + 2 * i + 1, 0);  // psi * C2
        for (usint k = 1; k < m_m; k++) {
          wPublicElementB(OutStart + i, j) +=
              wPublicElementB(InStart + 2 * i + 1, k) * psi(k, j);  // B2 * psi
          wCT(OutStart + i, j) +=
              psi(k, j) * wCT(InStart + 2 * i + 1, k);  // psi * C2
        }
      }

      for (usint j = 0; j < m_m; j++) {
        wPublicElementB(OutStart + i, j) =
            pubElemB(0, j) - wPublicElementB(OutStart + i, j);
        wCT(OutStart + i, j) = origCT(0, j) - wCT(OutStart + i, j);
      }
    }
  }

  for (usint j = 0; j < m_m; j++) {
    (*evalCT)(0, j) = wCT(gateCnt - 1, j);
  }

  (*evalAttributes) = wX[gateCnt - 1];
}

/* The encryption function takes public parameters A, B, and d, attribute values
 * x and the plaintext pt and generates the ciphertext pair c0 and c1 Note that
 * B is two dimensional array of ring elements (matrix); Each row corresponds
 * B_i for i = 0, 1, ... ell, where ell is the number of attributes
 */
template <class Element, class Element2>
void KPABE<Element, Element2>::Encrypt(
    const shared_ptr<ParmType> params, const Matrix<Element> &pubElemA,
    const Matrix<Element> &pubElemB,
    const Element &d,  // TBA
    const usint x[], const Element &ptext,
    DggType &dgg,  // to generate error terms (Gaussian)
    DugType &dug,  // select according to uniform distribution
    BinaryUniformGenerator
        &bug,                // select according to uniform distribution binary
    Matrix<Element> *ctCin,  // value set in this function
    Element *ctC1            // value set in this function
) {
  // compute c1 first
  Element s(dug, params, Format::COEFFICIENT);
  s.SwitchFormat();

  Element qHalf(params, Format::COEFFICIENT, true);
  qHalf += (m_q >> 1);
  qHalf.SwitchFormat();
  qHalf.AddILElementOne();

  Element err1 = Element(dgg, params, Format::COEFFICIENT);
  err1.SwitchFormat();

  *ctC1 = s * d + err1 + ptext * qHalf;

  // ***
  // Compute Cin
  auto zero_alloc = Element::Allocator(params, Format::EVALUATION);
  Matrix<Element> g = Matrix<Element>(zero_alloc, 1, m_k).GadgetVector(m_base);

  Matrix<Element> errA(Element::MakeDiscreteGaussianCoefficientAllocator(
                           params, Format::EVALUATION, SIGMA),
                       1, m_m);
  Matrix<Element> errCin(zero_alloc, 1, m_m);

  for (usint j = 0; j < m_m; j++) {
    (*ctCin)(0, j) = pubElemA(0, j) * s + errA(0, j);
  }
  for (usint i = 1; i < m_ell + 2; i++) {
    // Si values
    for (usint si = 0; si < m_m; si++) {
      errCin(0, si).SetValuesToZero();
      for (usint sj = 0; sj < m_m; sj++) {
        if (bug.GenerateInteger() == BigInteger(1))
          errCin(0, si) += errA(0, sj);
        else
          errCin(0, si) -= errA(0, sj);
      }
    }

    for (usint j = 0; j < m_k; j++) {
      if (x[i - 1] != 0)
        (*ctCin)(i, j) = (g(0, j) + pubElemB(i - 1, j)) * s + errCin(0, j);
      else
        (*ctCin)(i, j) = pubElemB(i - 1, j) * s + errCin(0, j);
    }
    (*ctCin)(i, m_m - 2) = pubElemB(i - 1, m_m - 2) * s + errCin(0, m_m - 2);
    (*ctCin)(i, m_m - 1) = pubElemB(i - 1, m_m - 1) * s + errCin(0, m_m - 1);
  }
}

/* Given public parameter d and a public key B,
it generates the corresponding secret key: skA for A and skB for B */
/* Note that only PKG can call this fcuntion as it needs the trapdoor T_A */
template <class Element, class Element2>
void KPABE<Element, Element2>::KeyGen(
    const shared_ptr<ParmType> params,
    const Matrix<Element>
        &pubElemA,  // Public parameter $A \in R_q^{1 \times w}$
    const Matrix<Element>
        &evalPubElemBf,  // Public parameter $B \in R_q^{ell \times k}$
    const Element &publicElemBeta,  // public key $d \in R_q$
    const RLWETrapdoorPair<Element>
        &secElemTA,  // Secret parameter $T_H \in R_q^{1 \times k} \times R_q^{1
                     // \times k}$
    DggType &dgg,    // to generate error terms (Gaussian)
    Matrix<Element> *sk  // Secret key
) {
  double s = SPECTRAL_BOUND(m_N, m_m - 2, m_base);
  //    Matrix<Element>
  // skB(Element::MakeDiscreteGaussianCoefficientAllocator(params,
  // Format::EVALUATION, SIGMA), m_m, 1);
  Matrix<Element> skB(Element::MakeDiscreteGaussianCoefficientAllocator(
                          params, Format::EVALUATION, s),
                      m_m, 1);

  Element newChallenge(params, Format::EVALUATION, true);
  for (usint j = 0; j < m_m; j++)
    newChallenge += (evalPubElemBf(0, j) * skB(j, 0));

  newChallenge = publicElemBeta - newChallenge;

  double c = (m_base + 1) * SIGMA;

  DggType dggLargeSigma = DggType(sqrt(s * s - c * c));

  Matrix<Element> skA(Element::Allocator(params, Format::EVALUATION), m_m, 1);
  skA = RLWETrapdoorUtility<Element>::GaussSamp(
      m_N, m_k, pubElemA, secElemTA, newChallenge, dgg, dggLargeSigma, m_base);

  for (usint i = 0; i < m_m; i++) (*sk)(0, i) = skA(i, 0);
  for (usint i = 0; i < m_m; i++) (*sk)(1, i) = skB(i, 0);
}

/*
 * Decryption function takes the ciphertext pair and the secret keys
 * and yields the decrypted plaintext in Format::COEFFICIENT form
 */
template <class Element, class Element2>
void KPABE<Element, Element2>::Decrypt(const shared_ptr<ParmType> params,
                                       const Matrix<Element> &sk,
                                       const Matrix<Element> &ctA,
                                       const Matrix<Element> &evalCT,
                                       const Element &ctC1, Element *dtext) {
  *dtext = ctA(0, 0) * sk(0, 0);
  for (usint i = 1; i < m_m; i++) *dtext += ctA(0, i) * sk(0, i);

  for (usint i = 0; i < m_m; i++) *dtext += evalCT(0, i) * sk(1, i);

  *dtext = ctC1 - *dtext;
  dtext->SwitchFormat();
}

/*
 * Decryption function takes the ciphertext pair and the secret keys
 * and yields the decrypted plaintext in Format::COEFFICIENT form
 */
template <class Element, class Element2>
void KPABE<Element, Element2>::Decode(Poly *dtext) {
  BigInteger dec, threshold = m_q >> 2, qHalf = m_q >> 1;

  for (usint i = 0; i < m_N; i++) {
    dec = dtext->at(i);

    if (dec > qHalf) dec = m_q - dec;

    if (dec > threshold)
      dtext->at(i) = BigInteger(1);
    else
      dtext->at(i) = BigInteger(0);
  }
}

/*
 * This is method for evaluating a single NAND gate
 */
template <class Element, class Element2>
void KPABE<Element, Element2>::NANDGateEvalCT(
    const shared_ptr<ParmType> ilParams, const Matrix<Element> &ctC0,
    const usint x[], const Matrix<Element> &origPubElem,
    const Matrix<Element> &origCT, usint *evalAttribute,
    Matrix<Element> *evalCT, uint32_t seed) {
  auto zero_alloc = Element::Allocator(ilParams, Format::EVALUATION);

  Matrix<Element> psi(zero_alloc, m_m, m_m);

  Matrix<Element> negB(zero_alloc, 1, m_m);  // EVALUATE (NTT domain)

  (*evalAttribute) = 1 - x[0] * x[1];  // Boolean output

  /* -B1 */
  for (usint j = 0; j < m_m; j++)  // Negating B1 for bit decomposition
    negB(0, j) = origPubElem(0, j).Negate();

  if (m_mode == NAF)
    PolyVec2BalDecom<Element>(ilParams, m_base, m_k, negB, &psi);
  else
    InverseRingVector<Element>(m_util, ilParams, negB, seed, &psi);

  psi.SwitchFormat();

  /* x2*C1 */
  for (usint i = 0; i < m_m; i++) {
    if (x[1] != 0)
      (*evalCT)(0, i) = origCT(0, i);
    else
      (*evalCT)(0, i).SetValuesToZero();
  }

  /* B2*Psi; Psi*C2 */
  for (usint i = 0; i < m_m; i++) {
    (*evalCT)(0, i) += psi(0, i) * origCT(1, 0);
    for (usint j = 1; j < m_m; j++) {
      (*evalCT)(0, i) += psi(j, i) * origCT(1, j);
    }
  }

  for (usint i = 0; i < m_m; i++) {
    (*evalCT)(0, i) = ctC0(0, i) - (*evalCT)(0, i);
  }
}

/*
 * This is method for evaluating a single NAND gate
 */
template <class Element, class Element2>
void KPABE<Element, Element2>::NANDGateEvalPK(
    const shared_ptr<ParmType> ilParams, const Matrix<Element> &pubElemB0,
    const Matrix<Element> &origPubElem, Matrix<Element> *evalPubElem,
    uint32_t seed) {
  auto zero_alloc = Element::Allocator(ilParams, Format::EVALUATION);

  Matrix<Element> psi(zero_alloc, m_m, m_m);

  Matrix<Element> negB(zero_alloc, 1, m_m);  // EVALUATE (NTT domain)

  /* -B1 */
  for (usint j = 0; j < m_m; j++)  // Negating B1 for bit decomposition
    negB(0, j) = origPubElem(0, j).Negate();

  if (m_mode == NAF)
    PolyVec2BalDecom<Element>(ilParams, m_base, m_k, negB, &psi);
  else
    InverseRingVector<Element>(m_util, ilParams, negB, seed, &psi);

  psi.SwitchFormat();

  /* B2*Psi; Psi*C2 */
  for (usint i = 0; i < m_m; i++) {
    (*evalPubElem)(0, i) = origPubElem(1, 0) * psi(0, i);
    for (usint j = 1; j < m_m; j++) {
      (*evalPubElem)(0, i) += origPubElem(1, j) * psi(j, i);
    }
  }

  for (usint i = 0; i < m_m; i++) {
    (*evalPubElem)(0, i) = pubElemB0(0, i) - (*evalPubElem)(0, i);
  }
}

template <class Element, class Element2>
void KPABE<Element, Element2>::NANDGateEvalPKDCRT(
    const shared_ptr<ParmType> params, const Matrix<Element> &pubElemB0,
    const Matrix<Element> &origPubElem, Matrix<Element> *evalPubElem,
    const shared_ptr<Parm2Type> ilParamsConsolidated, uint32_t seed) {
  auto zero_alloc = Element::Allocator(params, Format::EVALUATION);

  Matrix<Element> psi(zero_alloc, m_m, m_m);

  Matrix<Element> negB(zero_alloc, 1, m_m);  // EVALUATE (NTT domain)

  auto zero_alloc_poly =
      Element2::Allocator(ilParamsConsolidated, Format::COEFFICIENT);

  Matrix<Element2> psiPoly(zero_alloc_poly, m_m,
                           m_m);  // Needed for Bit Decomposition (BD) matrices

  Matrix<Element2> negBPolyMatrix(zero_alloc_poly, 1, m_m);

  /* -B1 */
  for (usint j = 0; j < m_m; j++)  // Negating B1 for bit decomposition
    negB(0, j) = origPubElem(0, j).Negate();

  for (usint i = 0; i < negB.GetRows(); i++) {
    for (usint j = 0; j < negB.GetCols(); j++) {
      negBPolyMatrix(i, j) = negB(i, j).CRTInterpolate();
    }
  }

  if (m_mode == NAF)
    PolyVec2BalDecom<Element2>(ilParamsConsolidated, m_base, m_k,
                               negBPolyMatrix, &psiPoly);
  else
    InverseRingVector<Element2>(m_util, ilParamsConsolidated, negBPolyMatrix,
                                seed, &psiPoly);

  // DCRT CREATE

  for (usint i = 0; i < psiPoly.GetRows(); i++) {
    for (usint j = 0; j < psiPoly.GetCols(); j++) {
      Element temp(psiPoly(i, j), params);
      psi(i, j) = temp;
    }
  }

  psi.SwitchFormat();

  /* B2*Psi; Psi*C2 */
  for (usint i = 0; i < m_m; i++) {
    (*evalPubElem)(0, i) = origPubElem(1, 0) * psi(0, i);
    for (usint j = 1; j < m_m; j++) {
      (*evalPubElem)(0, i) += origPubElem(1, j) * psi(j, i);
    }
  }

  for (usint i = 0; i < m_m; i++) {
    (*evalPubElem)(0, i) = pubElemB0(0, i) - (*evalPubElem)(0, i);
  }
}

template <class Element, class Element2>
void KPABE<Element, Element2>::NANDGateEvalCTDCRT(
    const shared_ptr<ParmType> params, const Matrix<Element> &ctC0,
    const usint x[], const Matrix<Element> &origPubElem,
    const Matrix<Element> &origCT, usint *evalAttribute,
    Matrix<Element> *evalCT, const shared_ptr<Parm2Type> ilParamsConsolidated,
    uint32_t seed) {
  auto zero_alloc = Element::Allocator(params, Format::EVALUATION);

  Matrix<Element> psi(zero_alloc, m_m, m_m);

  Matrix<Element> negB(zero_alloc, 1, m_m);  // EVALUATE (NTT domain)

  auto zero_alloc_poly =
      Element2::Allocator(ilParamsConsolidated, Format::COEFFICIENT);

  // Needed for Bit Decomposition (BD) matrices
  Matrix<Element2> psiPoly(zero_alloc_poly, m_m, m_m);

  Matrix<Element2> negBPolyMatrix(zero_alloc_poly, 1, m_m);

  (*evalAttribute) = 1 - x[0] * x[1];  // Boolean output

  /* -B1 */
  for (usint j = 0; j < m_m; j++)  // Negating B1 for bit decomposition
    negB(0, j) = origPubElem(0, j).Negate();

  for (usint i = 0; i < negB.GetRows(); i++) {
    for (usint j = 0; j < negB.GetCols(); j++) {
      negBPolyMatrix(i, j) = negB(i, j).CRTInterpolate();
    }
  }

  if (m_mode == NAF)
    PolyVec2BalDecom<Element2>(ilParamsConsolidated, m_base, m_k,
                               negBPolyMatrix, &psiPoly);
  else
    InverseRingVector<Element2>(m_util, ilParamsConsolidated, negBPolyMatrix,
                                seed, &psiPoly);

  // DCRT CREATE

  for (usint i = 0; i < psiPoly.GetRows(); i++) {
    for (usint j = 0; j < psiPoly.GetCols(); j++) {
      Element temp(psiPoly(i, j), params);
      psi(i, j) = temp;
    }
  }

  psi.SwitchFormat();

  /* x2*C1 */
  for (usint i = 0; i < m_m; i++) {
    if (x[1] != 0)
      (*evalCT)(0, i) = origCT(0, i);
    else
      (*evalCT)(0, i).SetValuesToZero();
  }

  /* B2*Psi; Psi*C2 */
  for (usint i = 0; i < m_m; i++) {
    (*evalCT)(0, i) += psi(0, i) * origCT(1, 0);
    for (usint j = 1; j < m_m; j++) {
      (*evalCT)(0, i) += psi(j, i) * origCT(1, j);
    }
  }

  for (usint i = 0; i < m_m; i++) {
    (*evalCT)(0, i) = ctC0(0, i) - (*evalCT)(0, i);
  }
}

template <class Element, class Element2>
void KPABE<Element, Element2>::ANDGateEvalPK(
    shared_ptr<ParmType> ilParams, const Matrix<Element> &origPubElemB,
    Matrix<Element> *evalPubElemBf, uint32_t seed) {
  auto zero_alloc = Element::Allocator(ilParams, Format::EVALUATION);
  Matrix<Element> psi(zero_alloc, m_m, m_m);
  Matrix<Element> negB(zero_alloc, 1, m_m);  // EVALUATE (NTT domain)

  /* -B1 */
  for (usint j = 0; j < m_m; j++) {  // Negating B1 for bit decomposition
    negB(0, j) = origPubElemB(0, j).Negate();
  }

  if (m_mode == NAF)
    PolyVec2BalDecom<Element>(ilParams, m_base, m_k, negB, &psi);
  else
    InverseRingVector<Element>(m_util, ilParams, negB, seed, &psi);

  psi.SwitchFormat();

  /* B2*Psi; Psi*C2 */
  for (usint i = 0; i < m_m; i++) {
    (*evalPubElemBf)(0, i) = origPubElemB(1, 0) * psi(0, i);
    for (usint j = 1; j < m_m; j++) {
      (*evalPubElemBf)(0, i) += origPubElemB(1, j) * psi(j, i);
    }
  }
}

template <class Element, class Element2>
void KPABE<Element, Element2>::ANDGateEvalCT(
    shared_ptr<ParmType> ilParams, const usint x[],
    const Matrix<Element> &origPubElemB, const Matrix<Element> &origCT,
    usint *evalAttribute, Matrix<Element> *evalCT, uint32_t seed) {
  auto zero_alloc = Element::Allocator(ilParams, Format::EVALUATION);
  Matrix<Element> psi(zero_alloc, m_m, m_m);
  Matrix<Element> negB(zero_alloc, 1, m_m);  // EVALUATE (NTT domain)
  (*evalAttribute) = x[0] * x[1];            // Boolean output
  /* -B1 */
  for (usint j = 0; j < m_m; j++) {  // Negating B1 for bit decomposition
    negB(0, j) = origPubElemB(0, j).Negate();
  }

  if (m_mode == NAF)
    PolyVec2BalDecom<Element>(ilParams, m_base, m_k, negB, &psi);
  else
    InverseRingVector<Element>(m_util, ilParams, negB, seed, &psi);

  psi.SwitchFormat();
  /* x2*C1 */
  for (usint i = 0; i < m_m; i++) {
    if (x[1] != 0)
      (*evalCT)(0, i) = origCT(0, i);
    else
      (*evalCT)(0, i).SetValuesToZero();
  }
  /* B2*Psi; Psi*C2 */
  for (usint i = 0; i < m_m; i++) {
    (*evalCT)(0, i) += psi(0, i) * origCT(1, 0);
    for (usint j = 1; j < m_m; j++) {
      (*evalCT)(0, i) += psi(j, i) * origCT(1, j);
    }
  }
}

template <class Element, class Element2>
void KPABE<Element, Element2>::ANDGateEvalPKDCRT(
    const shared_ptr<ParmType> params, const Matrix<Element> &origPubElemB,
    Matrix<Element> *evalPubElemBf,
    const shared_ptr<Parm2Type> ilParamsConsolidated, uint32_t seed) {
  auto zero_alloc = Element::Allocator(params, Format::EVALUATION);
  Matrix<Element> psi(zero_alloc, m_m, m_m);
  Matrix<Element> negB(zero_alloc, 1, m_m);  // EVALUATE (NTT domain)

  auto zero_alloc_poly =
      Element2::Allocator(ilParamsConsolidated, Format::COEFFICIENT);

  // Needed for Bit Decomposition (BD) matrices
  Matrix<Element2> psiPoly(zero_alloc_poly, m_m, m_m);

  Matrix<Element2> negBPolyMatrix(zero_alloc_poly, 1, m_m);

  /* -B1 */
  for (usint j = 0; j < m_m; j++) {  // Negating B1 for bit decomposition
    negB(0, j) = origPubElemB(0, j).Negate();
  }

  for (usint i = 0; i < negB.GetRows(); i++) {
    for (usint j = 0; j < negB.GetCols(); j++) {
      negBPolyMatrix(i, j) = negB(i, j).CRTInterpolate();
    }
  }

  if (m_mode == NAF)
    PolyVec2BalDecom<Element2>(ilParamsConsolidated, m_base, m_k,
                               negBPolyMatrix, &psiPoly);
  else
    InverseRingVector<Element2>(m_util, ilParamsConsolidated, negBPolyMatrix,
                                seed, &psiPoly);

  // DCRT CREATE

  for (usint i = 0; i < psiPoly.GetRows(); i++) {
    for (usint j = 0; j < psiPoly.GetCols(); j++) {
      Element temp(psiPoly(i, j), params);
      psi(i, j) = temp;
    }
  }

  psi.SwitchFormat();

  /* B2*Psi; Psi*C2 */
  for (usint i = 0; i < m_m; i++) {
    (*evalPubElemBf)(0, i) = origPubElemB(1, 0) * psi(0, i);
    for (usint j = 1; j < m_m; j++) {
      (*evalPubElemBf)(0, i) += origPubElemB(1, j) * psi(j, i);
    }
  }
}

template <class Element, class Element2>
void KPABE<Element, Element2>::ANDGateEvalCTDCRT(
    const shared_ptr<ParmType> params,
    const usint x[2],  // TBA
    const Matrix<Element> &origPubElemB, const Matrix<Element> &origCT,
    usint *evalAttribute, Matrix<Element> *evalCT,
    const shared_ptr<Parm2Type> ilParamsConsolidated, uint32_t seed) {
  auto zero_alloc = Element::Allocator(params, Format::EVALUATION);
  Matrix<Element> psi(zero_alloc, m_m, m_m);
  Matrix<Element> negB(zero_alloc, 1, m_m);  // EVALUATE (NTT domain)
  (*evalAttribute) = x[0] * x[1];            // Boolean output

  auto zero_alloc_poly =
      Element2::Allocator(ilParamsConsolidated, Format::COEFFICIENT);

  // Needed for Bit Decomposition (BD) matrices
  Matrix<Element2> psiPoly(zero_alloc_poly, m_m, m_m);

  Matrix<Element2> negBPolyMatrix(zero_alloc_poly, 1, m_m);

  /* -B1 */
  for (usint j = 0; j < m_m; j++) {  // Negating B1 for bit decomposition
    negB(0, j) = origPubElemB(0, j).Negate();
  }

  for (usint i = 0; i < negB.GetRows(); i++) {
    for (usint j = 0; j < negB.GetCols(); j++) {
      negBPolyMatrix(i, j) = negB(i, j).CRTInterpolate();
    }
  }

  if (m_mode == NAF)
    PolyVec2BalDecom<Element2>(ilParamsConsolidated, m_base, m_k,
                               negBPolyMatrix, &psiPoly);
  else
    InverseRingVector<Element2>(m_util, ilParamsConsolidated, negBPolyMatrix,
                                seed, &psiPoly);

  // DCRT CREATE

  for (usint i = 0; i < psiPoly.GetRows(); i++) {
    for (usint j = 0; j < psiPoly.GetCols(); j++) {
      Element temp(psiPoly(i, j), params);
      psi(i, j) = temp;
    }
  }

  psi.SwitchFormat();

  /* x2*C1 */
  for (usint i = 0; i < m_m; i++) {
    if (x[1] != 0)
      (*evalCT)(0, i) = origCT(0, i);
    else
      (*evalCT)(0, i).SetValuesToZero();
  }
  /* B2*Psi; Psi*C2 */
  for (usint i = 0; i < m_m; i++) {
    (*evalCT)(0, i) += psi(0, i) * origCT(1, 0);
    for (usint j = 1; j < m_m; j++) {
      (*evalCT)(0, i) += psi(j, i) * origCT(1, j);
    }
  }
}

}  // namespace lbcrypto
