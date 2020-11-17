// @file lweconjunctionchcprf.cpp Implementation of conjunction
// constraint-hiding constrained PRFs as described in
// https://eprint.iacr.org/2017/143.pdf
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

#ifndef LBCRYPTO_OBFUSCATE_LWECONJUNCTIONCHCPRF_CPP
#define LBCRYPTO_OBFUSCATE_LWECONJUNCTIONCHCPRF_CPP

#include "obfuscation/lweconjunctionchcprf.h"

namespace lbcrypto {

template <class Element>
LWEConjunctionCHCPRFAlgorithm<Element>::LWEConjunctionCHCPRFAlgorithm(
    usint base, usint chunkSize, usint length, usint n)
    : m_base(base),
      m_chunkSize(chunkSize),
      m_length(length),
      m_adjustedLength(length / chunkSize),
      m_chunkExponent(1 << m_chunkSize),
      m_dgg(SIGMA),
      m_A(std::make_shared<std::vector<Matrix<Element>>>()),
      m_T(std::make_shared<std::vector<RLWETrapdoorPair<Element>>>()) {
  // Generate ring parameters
  double qEst = EstimateRingModulus(n);
  m_elemParams = GenerateElemParams(qEst, n);

  double modulus = m_elemParams->GetModulus().ConvertToDouble();

  // Initialize m_dggLargeSigma
  usint k = floor(log2(modulus - 1.0) + 1.0);
  usint m = ceil(k / log2(base)) + 2;

  double c = (base + 1) * SIGMA;
  double s = SPECTRAL_BOUND(n, m - 2, base);

  if (sqrt(s * s - c * c) <= KARNEY_THRESHOLD)
    m_dggLargeSigma = typename Element::DggType(sqrt(s * s - c * c));
  else
    m_dggLargeSigma = m_dgg;

  size_t size = m_elemParams->GetParams().size();

  BigInteger Q(m_elemParams->GetModulus());

  vector<NativeInteger> moduli(size);
  vector<NativeInteger> roots(size);
  for (size_t i = 0; i < size; i++) {
    moduli[i] = m_elemParams->GetParams()[i]->GetModulus();
    roots[i] = m_elemParams->GetParams()[i]->GetRootOfUnity();
  }

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

  // Generate encoding keys
  EncodingParamsGen();
}

template <class Element>
usint LWEConjunctionCHCPRFAlgorithm<Element>::GetRingDimension() const {
  return m_elemParams->GetRingDimension();
}

template <class Element>
usint LWEConjunctionCHCPRFAlgorithm<Element>::GetLogModulus() const {
  double q = m_elemParams->GetModulus().ConvertToDouble();
  usint logModulus = floor(log2(q - 1.0) + 1.0);
  return logModulus;
}

template <class Element>
shared_ptr<vector<vector<Element>>>
LWEConjunctionCHCPRFAlgorithm<Element>::KeyGen() {
  auto s = std::make_shared<vector<vector<Element>>>();

  for (usint i = 0; i < m_adjustedLength; i++) {
    vector<Element> s_i;

    for (usint k = 0; k < m_chunkExponent; k++) {
      Element s_ik = Element(m_dgg, m_elemParams, Format::COEFFICIENT);
      s_ik.SwitchFormat();

      s_i.push_back(s_ik);
    }

    s->push_back(s_i);
  }

  return s;
}

template <class Element>
shared_ptr<vector<vector<shared_ptr<Matrix<Element>>>>>
LWEConjunctionCHCPRFAlgorithm<Element>::Constrain(
    const shared_ptr<vector<vector<Element>>> s, const std::string &pattern) {
  auto D = std::make_shared<vector<vector<shared_ptr<Matrix<Element>>>>>();

  for (usint i = 0; i < m_adjustedLength; i++) {
    // current chunk of cleartext pattern
    std::string chunk = pattern.substr(i * m_chunkSize, m_chunkSize);

    // build a chunk mask that maps "10??" to "1100" - ones correspond to
    // non-wildcard character
    std::string chunkTemp = replaceChar(chunk, '0', '1');
    chunkTemp = replaceChar(chunkTemp, '?', '0');
    // store the mask as integer for bitwise operations
    usint chunkMask = std::stoi(chunkTemp, nullptr, 2);

    // build a chunk target that maps "10??" to "1000" - replacing wildcard
    // character by 0
    chunkTemp = replaceChar(chunk, '?', '0');
    usint chunkTarget = std::stoi(chunkTemp, nullptr, 2);

    vector<shared_ptr<Matrix<Element>>> D_i;

    for (usint k = 0; k < m_chunkExponent; k++) {
      Element s_ik = (*s)[i][k];

      if ((k & chunkMask) != chunkTarget) {
        s_ik = Element(m_dgg, m_elemParams, Format::COEFFICIENT);
        s_ik.SwitchFormat();
      }

      shared_ptr<Matrix<Element>> D_ik = Encode(i, i + 1, s_ik);
      D_i.push_back(D_ik);
    }

    D->push_back(D_i);
  }

  return D;
}

template <class Element>
shared_ptr<vector<NativePoly>> LWEConjunctionCHCPRFAlgorithm<Element>::Evaluate(
    const shared_ptr<vector<vector<Element>>> s,
    const std::string &input) const {
  Element yCurrent;

  for (usint i = 0; i < m_adjustedLength; i++) {
    std::string chunk = input.substr(i * m_chunkSize, m_chunkSize);
    int k = std::stoi(chunk, nullptr, 2);

    if (i == 0)
      yCurrent = (*s)[i][k];
    else
      yCurrent *= (*s)[i][k];
  }

  Element y = (*m_A)[m_adjustedLength](0, 1) * yCurrent;

  return TransformMatrixToPRFOutput(y);
}

template <class Element>
shared_ptr<vector<NativePoly>> LWEConjunctionCHCPRFAlgorithm<Element>::Evaluate(
    const shared_ptr<vector<vector<shared_ptr<Matrix<Element>>>>> D,
    const std::string &input) const {
  Matrix<Element> y = (*m_A)[0];

  for (usint i = 0; i < m_adjustedLength; i++) {
    std::string chunk = input.substr(i * m_chunkSize, m_chunkSize);
    int k = std::stoi(chunk, nullptr, 2);

    y = y * *(*D)[i][k];
  }

  return TransformMatrixToPRFOutput(y(0, 1));
}

template <class Element>
double LWEConjunctionCHCPRFAlgorithm<Element>::EstimateRingModulus(usint n) {
  // smoothing parameter - also standard deviation for noise Elementnomials
  double sigma = SIGMA;

  // assurance measure
  double alpha = 36;

  // empirical parameter
  double beta = 6;

  // Bound of the Gaussian error Elementnomial
  double Berr = sigma * sqrt(alpha);

  // probability of hitting the "danger" zone that affects the rounding result
  double Pe = 1 << 20;

  uint32_t length = m_adjustedLength;
  uint32_t base = m_base;

  // Correctness constraint
  auto qCorrectness = [&](uint32_t n, uint32_t m, uint32_t k) -> double {
    return 1024 * Pe * Berr *
           pow(sqrt(m * n) * beta * SPECTRAL_BOUND(n, m - 2, base), length - 1);
  };

  double qPrev = 1e6;
  double q = 0;
  usint k = 0;
  usint m = 0;

  // initial value
  k = floor(log2(qPrev - 1.0) + 1.0);
  m = ceil(k / log2(base)) + 2;
  q = qCorrectness(n, m, k);

  // get a more accurate value of q
  while (std::abs(q - qPrev) > 0.001 * q) {
    qPrev = q;
    k = floor(log2(qPrev - 1.0) + 1.0);
    m = ceil(k / log2(base)) + 2;
    q = qCorrectness(n, m, k);
  }

  return q;
}

template <class Element>
void LWEConjunctionCHCPRFAlgorithm<Element>::EncodingParamsGen() {
  const shared_ptr<typename Element::Params> params = m_elemParams;
  usint base = m_base;
  usint stddev = m_dgg.GetStd();

  for (size_t i = 0; i <= m_adjustedLength; i++) {
    std::pair<Matrix<Element>, RLWETrapdoorPair<Element>> trapPair =
        RLWETrapdoorUtility<Element>::TrapdoorGen(params, stddev,
                                                  base);  // TODO remove stddev

    m_A->push_back(trapPair.first);
    m_T->push_back(trapPair.second);
  }
}

template <class Element>
shared_ptr<Matrix<Element>> LWEConjunctionCHCPRFAlgorithm<Element>::Encode(
    usint i, usint j, const Element &elem) {
  const Matrix<Element> Ai = (*m_A)[i];
  const Matrix<Element> Aj = (*m_A)[j];
  const RLWETrapdoorPair<Element> Ti = (*m_T)[i];

  size_t m = Ai.GetCols();
  size_t k = m - 2;
  size_t n = GetRingDimension();
  auto zero_alloc = Element::Allocator(elem.GetParams(), Format::EVALUATION);

  // generate a row vector of discrete Gaussian ring elements
  // YSP this can be done using discrete Gaussian allocator later - after the
  // dgg allocator is updated to use the same dgg instance DBC all the following
  // have insignificant timing
  Matrix<Element> ej(zero_alloc, 1, m);
#pragma omp parallel for
  for (size_t i = 0; i < m; i++) {
    ej(0, i) = Element(m_dgg, elem.GetParams(), Format::COEFFICIENT);
    ej(0, i).SwitchFormat();
  }

  const Matrix<Element> &bj = Aj.ScalarMult(elem) + ej;

  auto result = std::make_shared<Matrix<Element>>(zero_alloc, m, m);

// DBC: this loop takes all the time in encode
// TODO (dcousins): move gaussj generation out of the loop to enable
// parallelisation
#pragma omp parallel for schedule(dynamic)
  for (size_t i = 0; i < m; i++) {
    // the following takes approx 250 msec
    const Matrix<Element> &gaussj = RLWETrapdoorUtility<Element>::GaussSamp(
        n, k, Ai, Ti, bj(0, i), m_dgg, m_dggLargeSigma, m_base);
    // the following takes no time
    for (size_t j = 0; j < m; j++) {
      (*result)(j, i) = gaussj(j, 0);
    }
  }

  return result;
}

template <class Element>
shared_ptr<vector<NativePoly>>
LWEConjunctionCHCPRFAlgorithm<Element>::TransformMatrixToPRFOutput(
    const Element &input) const {
  auto result = std::make_shared<vector<NativePoly>>(1);

  auto element = input;
  element.SwitchFormat();

  (*result)[0] = element.ScaleAndRound(
      NativeInteger(2), m_tQHatInvModqDivqModt, m_tQHatInvModqDivqModtPrecon,
      m_tQHatInvModqBDivqModt, m_tQHatInvModqBDivqModtPrecon,
      m_tQHatInvModqDivqFrac, m_tQHatInvModqBDivqFrac);

  return result;
}

}  // namespace lbcrypto

#endif
