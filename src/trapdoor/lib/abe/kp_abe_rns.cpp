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
// This code provides functionality for KP_ABE. The algorithms and naming
// conventions can be found from this paper:
// https://eprint.iacr.org/2017/601.pdf

#define PROFILE

#include "abe/kp_abe_rns.h"
#include "utils/debug.h"

namespace lbcrypto {

/*
 * This is a setup function for Private Key Generator (PKG);
 * generates master public key (MPK) and master secret key
 * m_ell is the number of attributes
 */

void KPABErns::Setup(const shared_ptr<ParmType> params, int32_t base,
                     usint ell,     // number of attributes
                     DugType &dug,  // select according to uniform distribution
                     Matrix<DCRTPoly> *pubElemB) {
  Setup(params, base, ell);

  for (usint i = 0; i < (*pubElemB).GetRows(); i++)
    for (usint j = 0; j < (*pubElemB).GetCols(); j++) {
      (*pubElemB)(i, j) = DCRTPoly(dug, params, Format::COEFFICIENT);
      (*pubElemB)(i, j)
          .SwitchFormat();  // always kept in Format::EVALUATION format
    }
}

/**
 * This setup function is used by users; namely senders and receivers
 * Initialize private members of the object such as modulus, cyclotomic order,
 * etc. m_ell is the number of attributes
 */

void KPABErns::Setup(const shared_ptr<ParmType> params, int32_t base,
                     const usint ell) {
  m_N = params->GetCyclotomicOrder() >> 1;
  BigInteger Q(params->GetModulus());
  m_q = Q;
  m_base = base;

  size_t size = params->GetParams().size();

  size_t digitCount =
      (long)ceil(log2(params->GetParams()[0]->GetModulus().ConvertToDouble()) /
                 log2(base));
  m_k = digitCount * size;

  m_m = m_k + 2;

  m_ell = ell;

  vector<NativeInteger> moduli(size);
  vector<NativeInteger> roots(size);
  for (size_t i = 0; i < size; i++) {
    moduli[i] = params->GetParams()[i]->GetModulus();
    roots[i] = params->GetParams()[i]->GetRootOfUnity();
  }

  for (size_t i = 0; i < size; i++)
    m_util.push_back(LatticeSubgaussianUtility<NativeInteger>(m_base, moduli[i],
                                                              digitCount));

  NativeInteger t(2);

  const BigInteger QDivt = Q.DividedBy(t);

  m_QDivtModq.resize(size);
  for (size_t i = 0; i < size; i++) {
    BigInteger qi(moduli[i].ConvertToInt());
    BigInteger QDivtModqi = QDivt.Mod(qi);
    m_QDivtModq[i] = NativeInteger(QDivtModqi.ConvertToInt());
  }

  usint qMSB = moduli[0].GetMSB();
  usint sizeQMSB = GetMSB64(size);
  m_tQHatInvModqDivqModt.resize(size);
  m_tQHatInvModqDivqModtPrecon.resize(size);
  m_tQHatInvModqDivqFrac.resize(size);
  if (qMSB + sizeQMSB < 52) {
    for (size_t i = 0; i < size; i++) {
      BigInteger qi(moduli[i].ConvertToInt());
      BigInteger tQHatInvModqi =
          ((Q.DividedBy(qi)).ModInverse(qi) * BigInteger(t));
      BigInteger tQHatInvModqDivqi = tQHatInvModqi.DividedBy(qi);
      m_tQHatInvModqDivqModt[i] = tQHatInvModqDivqi.Mod(t).ConvertToInt();
      m_tQHatInvModqDivqModtPrecon[i] =
          m_tQHatInvModqDivqModt[i].PrepModMulConst(t);

      int64_t numerator = tQHatInvModqi.Mod(qi).ConvertToInt();
      int64_t denominator = moduli[i].ConvertToInt();
      m_tQHatInvModqDivqFrac[i] =
          static_cast<double>(numerator) / static_cast<double>(denominator);
    }
  } else {
    m_tQHatInvModqBDivqModt.resize(size);
    m_tQHatInvModqBDivqModtPrecon.resize(size);
    m_tQHatInvModqBDivqFrac.resize(size);
    usint qMSBHf = qMSB >> 1;
    for (size_t i = 0; i < size; i++) {
      BigInteger qi(moduli[i].ConvertToInt());
      BigInteger tQHatInvModqi =
          ((Q.DividedBy(qi)).ModInverse(qi) * BigInteger(t));
      BigInteger tQHatInvModqDivqi = tQHatInvModqi.DividedBy(qi);
      m_tQHatInvModqDivqModt[i] = tQHatInvModqDivqi.Mod(t).ConvertToInt();
      m_tQHatInvModqDivqModtPrecon[i] =
          m_tQHatInvModqDivqModt[i].PrepModMulConst(t);

      int64_t numerator = tQHatInvModqi.Mod(qi).ConvertToInt();
      int64_t denominator = moduli[i].ConvertToInt();
      m_tQHatInvModqDivqFrac[i] =
          static_cast<double>(numerator) / static_cast<double>(denominator);

      tQHatInvModqi.LShiftEq(qMSBHf);
      tQHatInvModqDivqi = tQHatInvModqi.DividedBy(qi);
      m_tQHatInvModqBDivqModt[i] = tQHatInvModqDivqi.Mod(t).ConvertToInt();
      m_tQHatInvModqBDivqModtPrecon[i] =
          m_tQHatInvModqBDivqModt[i].PrepModMulConst(t);

      numerator = tQHatInvModqi.Mod(qi).ConvertToInt();
      m_tQHatInvModqBDivqFrac[i] =
          static_cast<double>(numerator) / static_cast<double>(denominator);
    }
  }
}

/*
 * Given public parameters, attribute values and ciphertexts corresponding to
 * attributes, computes the ciphertext and the public key evalPubDCRTPoly for
 * the circuit of attributes m_ell is the number of attributes and the circuit
 * is assumed to be a binary tree of NAND gates Thus, m_ell must be a power of
 * two
 */

void KPABErns::EvalPK(const shared_ptr<ParmType> params,
                      const Matrix<DCRTPoly> &pubElemB,
                      Matrix<DCRTPoly> *evalPubElemBf, uint32_t seed) {
  TimeVar t1, t2;
  double offTime(0.0);
  double subgaussianTime(0.0);

  auto zero_alloc = DCRTPoly::Allocator(params, Format::EVALUATION);

  usint gateCnt = m_ell - 1;

  // Matrix<DCRTPoly> psi(zero_alloc, m_m, m_m); // Needed for bit decomposition
  // matrices
  // w stands for wire
  Matrix<DCRTPoly> wpublicElementB(
      zero_alloc, gateCnt,
      m_m);  // Bis associated with internal wires of the circuit
  // Temporary variables for bit decomposition operation
  Matrix<DCRTPoly> negPubElemB(zero_alloc, 1,
                               m_m);  // Format::EVALUATION (NTT domain)
  // std::vector<DCRTPoly> digitsC1(m_m);

  // Input level of the circuit
  usint t = m_ell >> 1;  // the number of the gates in the first level (the
                         // number of input gates)

  // looping to evaluate and calculate w, wB, wC
  // and R for all first level input gates
  for (usint i = 0; i < t; i++) {
#pragma omp parallel for schedule(dynamic)
    for (usint j = 0; j < m_m; j++) {  // Negating Bis for bit decomposition
      negPubElemB(0, j) = pubElemB(2 * i + 1, j).Negate();
      negPubElemB(0, j).SwitchFormat();
    }

    TIC(t1);
    TIC(t2);
    auto psi = InverseRingVectorDCRT(m_util, negPubElemB, 1);
    subgaussianTime += TOC_US(t1);

    psi->SwitchFormat();
    offTime += TOC_US(t2);

    /* Psi^T*C2 and B2*Psi */
#pragma omp parallel for schedule(dynamic)
    for (usint j = 0; j < m_m;
         j++) {  // the following two for loops are for vector matrix
                 // multiplication (a.k.a B(i+1) * BitDecompose(-Bi) and  gamma
                 // (0, 2) (for the second attribute of the circuit) *
                 // bitDecompose(-B))
      wpublicElementB(i, j) =
          pubElemB(2 * i + 2, 0) * (*psi)(0, j);  // B2 * BD(-Bi)
      for (usint k = 1; k < m_m; k++) {
        wpublicElementB(i, j) += pubElemB(2 * i + 2, k) * (*psi)(k, j);
      }
    }

#pragma omp parallel for schedule(dynamic)
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
    // Starting index for the input wires in level d
    usint inStart = m_ell - (m_ell >> (d - 1));
    // Starting index for the output wires in level d
    usint outStart = m_ell - (m_ell >> d);
    usint gCntinLeveld = m_ell >> (d + 1);  // number of gates in level d

    for (usint i = 0; i < gCntinLeveld; i++) {
#pragma omp parallel for schedule(dynamic)
      for (usint j = 0; j < m_m; j++) {
        negPubElemB(0, j) = wpublicElementB(inStart + 2 * i, j).Negate();
        negPubElemB(0, j).SwitchFormat();
      }

      TIC(t1);
      TIC(t2);
      auto psi = InverseRingVectorDCRT(m_util, negPubElemB, 1);
      subgaussianTime += TOC_US(t1);

      psi->SwitchFormat();
      offTime += TOC_US(t2);

#pragma omp parallel for schedule(dynamic)
      for (usint j = 0; j < m_m; j++) {
        wpublicElementB(outStart + i, j) =
            wpublicElementB(inStart + 2 * i + 1, 0) * (*psi)(0, j);  // B2 * Psi
        for (usint k = 1; k < m_m; k++) {
          wpublicElementB(outStart + i, j) +=
              wpublicElementB(inStart + 2 * i + 1, k) *
              (*psi)(k, j);  // B2 * Psi
        }
      }

#pragma omp parallel for schedule(dynamic)
      for (usint j = 0; j < m_m; j++) {
        wpublicElementB(outStart + i, j) =
            pubElemB(0, j) - wpublicElementB(outStart + i, j);
      }
    }
  }

#pragma omp parallel for schedule(dynamic)
  for (usint j = 0; j < m_m; j++) {
    (*evalPubElemBf)(0, j) = wpublicElementB(gateCnt - 1, j);
  }

  // std::cerr << "Computation of G^(-1):\t" << subgaussianTime/1000 <<
  // std::endl; std::cerr << "Computation of G^(-1) with NTT:\t" << offTime/1000
  // << std::endl;
}

/*
 * Given public parameters, attribute values and ciphertexts corresponding to
 * attributes, computes the ciphertext and the public key Bf for the circuit of
 * attributes m_ell is the number of attributes and the circuit is assumed to be
 * a binary tree of NAND gates Thus, m_ell must be a power of two
 */

void KPABErns::EvalCT(const shared_ptr<ParmType> params,
                      const Matrix<DCRTPoly> &pubElemB,
                      const usint x[],  // Attributes
                      const Matrix<DCRTPoly> &origCT, usint *evalAttributes,
                      Matrix<DCRTPoly> *evalCT, uint32_t seed) {
  // Part pertaining to A (does not change)
  for (usint i = 0; i < m_m; i++) (*evalCT)(0, i) = origCT(0, i);

  auto zero_alloc = DCRTPoly::Allocator(params, Format::EVALUATION);

  usint gateCnt = m_ell - 1;
  // Matrix<DCRTPoly> psi(zero_alloc, m_m, m_m);
  // w stands for Wire
  Matrix<DCRTPoly> wPublicElementB(
      zero_alloc, gateCnt,
      m_m);  // Bis associated with internal wires of the circuit
  Matrix<DCRTPoly> wCT(
      zero_alloc, gateCnt,
      m_m);  // Ciphertexts associated with internal wires of the circuit

  // Attribute values associated with internal wires of the circuit
  std::vector<usint> wX(gateCnt);

  // Temporary variables for bit decomposition operation
  Matrix<DCRTPoly> negB(zero_alloc, 1, m_m);  // Format::EVALUATION (NTT domain)

  // Input level of the circuit
  usint t = m_ell >> 1;  // the number of the gates in the first level (the
                         // number of input gates)

  // looping to evaluate and calculate w, wB, wC
  // and R for all first level input gates
  for (usint i = 0; i < t; i++) {
    wX[i] =
        x[0] - x[2 * i + 1] * x[2 * i + 2];  // calculating binary wire value

#pragma omp parallel for schedule(dynamic)
    for (usint j = 0; j < m_m; j++) {  // Negating Bis for bit decomposition
      negB(0, j) = pubElemB(2 * i + 1, j).Negate();
      negB(0, j).SwitchFormat();
    }

    auto psi = InverseRingVectorDCRT(m_util, negB, 1);

    psi->SwitchFormat();
    /*Starting computation for a NAND circuit*/
    /* x2 * C1 */
#pragma omp parallel for schedule(dynamic)
    for (usint j = 0; j < m_m; j++) {
      if (x[2 * i + 2] != 0)
        wCT(i, j) = origCT(2 * i + 1, j);
      else
        wCT(i, j).SetValuesToZero();
    }

    /* Psi^T*C2 and B2*Psi */
#pragma omp parallel for schedule(dynamic)
    for (usint j = 0; j < m_m;
         j++) {  // the following two for loops are for vector matrix
                 // multiplication (a.k.a B(i+1) * BitDecompose(-Bi) and  gamma
                 // (0, 2) (for the second attribute of the circuit) *
                 // bitDecompose(-B))
      wPublicElementB(i, j) =
          pubElemB(2 * i + 2, 0) * (*psi)(0, j);         // B2 * BD(-Bi)
      wCT(i, j) += (*psi)(0, j) * origCT(2 * i + 2, 0);  // BD(-Bi)*C2
      for (usint k = 1; k < m_m; k++) {
        wPublicElementB(i, j) += pubElemB(2 * i + 2, k) * (*psi)(k, j);
        wCT(i, j) += (*psi)(k, j) * origCT(2 * i + 2, k);
      }
    }

    /* B0 - B2*R and C0 - x2*C1 - C2*R */
#pragma omp parallel for schedule(dynamic)
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
    usint InStart =
        m_ell -
        (m_ell >> (d - 1));  // Starting index for the input wires in level d
    usint OutStart =
        m_ell - (m_ell >> d);  // Starting index for the output wires in level d
    usint gCntinLeveld = m_ell >> (d + 1);  // number of gates in level d

    for (usint i = 0; i < gCntinLeveld; i++) {
      wX[OutStart + i] = x[0] - wX[InStart + 2 * i] * wX[InStart + 2 * i + 1];

#pragma omp parallel for schedule(dynamic)
      for (usint j = 0; j < m_m; j++) {
        negB(0, j) = wPublicElementB(InStart + 2 * i, j).Negate();
        negB(0, j).SwitchFormat();
      }

      auto psi = InverseRingVectorDCRT(m_util, negB, 1);

      psi->SwitchFormat();

      // x2*C1
#pragma omp parallel for schedule(dynamic)
      for (usint j = 0; j < m_m; j++) {
        if (wX[InStart + 2 * i + 1] != 0)
          wCT(OutStart + i, j) = wCT(InStart + 2 * i, j);
        else
          wCT(OutStart + i, j).SetValuesToZero();
      }

#pragma omp parallel for schedule(dynamic)
      for (usint j = 0; j < m_m; j++) {
        wPublicElementB(OutStart + i, j) =
            wPublicElementB(InStart + 2 * i + 1, 0) * (*psi)(0, j);  // B2 * psi
        wCT(OutStart + i, j) +=
            (*psi)(0, j) * wCT(InStart + 2 * i + 1, 0);  // psi * C2
        for (usint k = 1; k < m_m; k++) {
          wPublicElementB(OutStart + i, j) +=
              wPublicElementB(InStart + 2 * i + 1, k) *
              (*psi)(k, j);  // B2 * psi
          wCT(OutStart + i, j) +=
              (*psi)(k, j) * wCT(InStart + 2 * i + 1, k);  // psi * C2
        }
      }

#pragma omp parallel for schedule(dynamic)
      for (usint j = 0; j < m_m; j++) {
        wPublicElementB(OutStart + i, j) =
            pubElemB(0, j) - wPublicElementB(OutStart + i, j);
        wCT(OutStart + i, j) = origCT(0, j) - wCT(OutStart + i, j);
      }
    }
  }

#pragma omp parallel for schedule(dynamic)
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

void KPABErns::Encrypt(
    const shared_ptr<ParmType> params, const Matrix<DCRTPoly> &pubElemA,
    const Matrix<DCRTPoly> &pubElemB,
    const DCRTPoly &d,  // TBA
    const usint x[], const NativePoly &ptext,
    DggType &dgg,             // to generate error terms (Gaussian)
    DugType &dug,             // select according to uniform distribution
    BugType &bug,             // select according to uniform distribution binary
    Matrix<DCRTPoly> *ctCin,  // value set in this function
    DCRTPoly *ctC1            // value set in this function
) {
  // compute c1 first
  DCRTPoly s(dug, params, Format::COEFFICIENT);
  s.SwitchFormat();

  DCRTPoly dtext(ptext, params);
  dtext.SwitchFormat();

  const std::vector<NativeInteger> &deltaTable = m_QDivtModq;

  DCRTPoly err1 = DCRTPoly(dgg, params, Format::COEFFICIENT);
  err1.SwitchFormat();

  *ctC1 = s * d + err1 + dtext.Times(deltaTable);

  // ***
  // Compute Cin
  auto zero_alloc = DCRTPoly::Allocator(params, Format::EVALUATION);
  Matrix<DCRTPoly> g =
      Matrix<DCRTPoly>(zero_alloc, 1, m_k).GadgetVector(m_base);

  Matrix<DCRTPoly> errA(DCRTPoly::MakeDiscreteGaussianCoefficientAllocator(
                            params, Format::EVALUATION, SIGMA),
                        1, m_m);
  Matrix<DCRTPoly> errCin(zero_alloc, 1, m_m);

#pragma omp parallel for schedule(dynamic)
  for (usint j = 0; j < m_m; j++) {
    (*ctCin)(0, j) = pubElemA(0, j) * s + errA(0, j);
  }

  for (usint i = 1; i < m_ell + 2; i++) {
    // Si values
#pragma omp parallel for schedule(dynamic)
    for (usint si = 0; si < m_m; si++) {
      errCin(0, si).SetValuesToZero();
      for (usint sj = 0; sj < m_m; sj++) {
        if (bug.GenerateInteger() == NativeInteger(1))
          errCin(0, si) += errA(0, sj);
        else
          errCin(0, si) -= errA(0, sj);
      }
    }

#pragma omp parallel for schedule(dynamic)
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

void KPABErns::KeyGen(
    const shared_ptr<ParmType> params,
    const Matrix<DCRTPoly>
        &pubElemA,  // Public parameter $A \in R_q^{1 \times w}$
    const Matrix<DCRTPoly>
        &evalPubElemBf,  // Public parameter $B \in R_q^{ell \times k}$
    const DCRTPoly &publicElemBeta,  // public key $d \in R_q$
    const RLWETrapdoorPair<DCRTPoly>
        &secElemTA,  // Secret parameter $T_H \in R_q^{1 \times k} \times R_q^{1
                     // \times k}$
    DggType &dgg,    // to generate error terms (Gaussian)
    Matrix<DCRTPoly> *sk  // Secret key
) {
  double s = SPECTRAL_BOUND(m_N, m_m - 2, m_base);
  //    Matrix<DCRTPoly>
  // skB(DCRTPoly::MakeDiscreteGaussianCoefficientAllocator(params,
  // Format::EVALUATION, SIGMA), m_m, 1);
  Matrix<DCRTPoly> skB(DCRTPoly::MakeDiscreteGaussianCoefficientAllocator(
                           params, Format::EVALUATION, s),
                       m_m, 1);

  DCRTPoly newChallenge(params, Format::EVALUATION, true);
  for (usint j = 0; j < m_m; j++)
    newChallenge += (evalPubElemBf(0, j) * skB(j, 0));

  newChallenge = publicElemBeta - newChallenge;

  double c = (m_base + 1) * SIGMA;

  DggType dggLargeSigma = DggType(sqrt(s * s - c * c));

  Matrix<DCRTPoly> skA(DCRTPoly::Allocator(params, Format::EVALUATION), m_m, 1);
  skA = RLWETrapdoorUtility<DCRTPoly>::GaussSamp(
      m_N, m_k, pubElemA, secElemTA, newChallenge, dgg, dggLargeSigma, m_base);

#pragma omp parallel for schedule(dynamic)
  for (usint i = 0; i < m_m; i++) (*sk)(0, i) = skA(i, 0);
#pragma omp parallel for schedule(dynamic)
  for (usint i = 0; i < m_m; i++) (*sk)(1, i) = skB(i, 0);
}

/*
 * Decryption function takes the ciphertext pair and the secret keys
 * and yields the decrypted plaintext in Format::COEFFICIENT form
 */

void KPABErns::Decrypt(const shared_ptr<ParmType> params,
                       const Matrix<DCRTPoly> &sk, const Matrix<DCRTPoly> &ctA,
                       const Matrix<DCRTPoly> &evalCT, const DCRTPoly &ctC1,
                       NativePoly *ptext) {
  DCRTPoly dtext = ctA(0, 0) * sk(0, 0);
  for (usint i = 1; i < m_m; i++) dtext += ctA(0, i) * sk(0, i);

  for (usint i = 0; i < m_m; i++) dtext += evalCT(0, i) * sk(1, i);

  dtext = ctC1 - dtext;
  dtext.SwitchFormat();

  // this is the resulting vector of coefficients;
  *ptext = dtext.ScaleAndRound(
      NativeInteger(2), m_tQHatInvModqDivqModt, m_tQHatInvModqDivqModtPrecon,
      m_tQHatInvModqBDivqModt, m_tQHatInvModqBDivqModtPrecon,
      m_tQHatInvModqDivqFrac, m_tQHatInvModqBDivqFrac);
}

void KPABErns::NANDGateEvalPK(const shared_ptr<ParmType> params,
                              const Matrix<DCRTPoly> &pubElemB,
                              Matrix<DCRTPoly> *evalPubElem, uint32_t seed) {
  auto zero_alloc = DCRTPoly::Allocator(params, Format::EVALUATION);

  Matrix<DCRTPoly> negB(zero_alloc, 1, m_m);  // EVALUATE (NTT domain)

  /* -B1 */
  for (usint j = 0; j < m_m; j++)  // Negating B1 for bit decomposition
    negB(0, j) = pubElemB(1, j).Negate();

  auto psi = InverseRingVectorDCRT(m_util, negB, 1);

  psi->SwitchFormat();

  /* B2*Psi; Psi*C2 */
  for (usint i = 0; i < m_m; i++) {
    (*evalPubElem)(0, i) = pubElemB(2, 0) * (*psi)(0, i);
    for (usint j = 1; j < m_m; j++) {
      (*evalPubElem)(0, i) += pubElemB(2, j) * (*psi)(j, i);
    }
  }

  for (usint i = 0; i < m_m; i++) {
    (*evalPubElem)(0, i) = pubElemB(0, i) - (*evalPubElem)(0, i);
  }
}

void KPABErns::NANDGateEvalCT(const shared_ptr<ParmType> params,
                              const Matrix<DCRTPoly> &origPubElem,
                              const usint x[], const Matrix<DCRTPoly> &origCT,
                              usint *evalAttribute, Matrix<DCRTPoly> *evalCT,
                              uint32_t seed) {
  auto zero_alloc = DCRTPoly::Allocator(params, Format::EVALUATION);

  Matrix<DCRTPoly> negB(zero_alloc, 1, m_m);  // EVALUATE (NTT domain)

  (*evalAttribute) = 1 - x[0] * x[1];  // Boolean output

  /* -B1 */
  for (usint j = 0; j < m_m; j++)  // Negating B1 for bit decomposition
    negB(0, j) = origPubElem(1, j).Negate();

  auto psi = InverseRingVectorDCRT(m_util, negB, 1);

  psi->SwitchFormat();

  /* x2*C1 */
  for (usint i = 0; i < m_m; i++) {
    if (x[1] != 0)
      (*evalCT)(0, i) = origCT(1, i);
    else
      (*evalCT)(0, i).SetValuesToZero();
  }

  /* B2*Psi; Psi*C2 */
  for (usint i = 0; i < m_m; i++) {
    (*evalCT)(0, i) += (*psi)(0, i) * origCT(2, 0);
    for (usint j = 1; j < m_m; j++) {
      (*evalCT)(0, i) += (*psi)(j, i) * origCT(2, j);
    }
  }

  for (usint i = 0; i < m_m; i++) {
    (*evalCT)(0, i) = origCT(0, i) - (*evalCT)(0, i);
  }
}

void KPABErns::ANDGateEvalPK(const shared_ptr<ParmType> params,
                             const Matrix<DCRTPoly> &origPubElemB,
                             Matrix<DCRTPoly> *evalPubElemBf, uint32_t seed) {
  auto zero_alloc = DCRTPoly::Allocator(params, Format::EVALUATION);
  Matrix<DCRTPoly> negB(zero_alloc, 1, m_m);  // EVALUATE (NTT domain)

  /* -B1 */
  for (usint j = 0; j < m_m; j++) {  // Negating B1 for bit decomposition
    negB(0, j) = origPubElemB(1, j).Negate();
  }

  auto psi = InverseRingVectorDCRT(m_util, negB, 1);

  psi->SwitchFormat();

  /* B2*Psi; Psi*C2 */
  for (usint i = 0; i < m_m; i++) {
    (*evalPubElemBf)(0, i) = origPubElemB(2, 0) * (*psi)(0, i);
    for (usint j = 1; j < m_m; j++) {
      (*evalPubElemBf)(0, i) += origPubElemB(2, j) * (*psi)(j, i);
    }
  }
}

void KPABErns::ANDGateEvalCT(const shared_ptr<ParmType> params,
                             const Matrix<DCRTPoly> &origPubElemB,
                             const usint x[2],  // TBA
                             const Matrix<DCRTPoly> &origCT,
                             usint *evalAttribute, Matrix<DCRTPoly> *evalCT,
                             uint32_t seed) {
  auto zero_alloc = DCRTPoly::Allocator(params, Format::EVALUATION);
  Matrix<DCRTPoly> negB(zero_alloc, 1, m_m);  // EVALUATE (NTT domain)
  (*evalAttribute) = x[0] * x[1];             // Boolean output

  /* -B1 */
  for (usint j = 0; j < m_m; j++) {  // Negating B1 for bit decomposition
    negB(0, j) = origPubElemB(1, j).Negate();
  }

  auto psi = InverseRingVectorDCRT(m_util, negB, 1);

  psi->SwitchFormat();

  /* x2*C1 */
  for (usint i = 0; i < m_m; i++) {
    if (x[1] != 0)
      (*evalCT)(0, i) = origCT(0, i);
    else
      (*evalCT)(0, i).SetValuesToZero();
  }
  /* B2*Psi; Psi*C2 */
  for (usint i = 0; i < m_m; i++) {
    (*evalCT)(0, i) += (*psi)(0, i) * origCT(1, 0);
    for (usint j = 1; j < m_m; j++) {
      (*evalCT)(0, i) += (*psi)(j, i) * origCT(1, j);
    }
  }
}

}  // namespace lbcrypto
