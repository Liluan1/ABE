// @file TODO
// @author TPOC: contact@palisade-crypto.org
//
// @copyright Copyright (c) TODO

#define PROFILE  // define this to enable PROFILELOG and TIC/TOC

#include <fstream>
#include <iostream>
#include "gsw.h"
#include "subgaussian.h"

#include "utils/debug.h"

using namespace lbcrypto;

shared_ptr<Matrix<DCRTPoly>> InverseG(const Matrix<DCRTPoly> &A, usint base);

int main() {
  usint n = 1024;  // cyclotomic order
  size_t kRes = 60;

  size_t depth = 2;

  usint base = 2;

  size_t size = 2;

  std::cout << "n: " << n << std::endl;

  // double sigma = SIGMA;

  std::vector<NativeInteger> moduli;
  std::vector<NativeInteger> roots_Of_Unity;

  // makes sure the first integer is less than 2^60-1 to take advangate of NTL
  // optimizations
  NativeInteger firstInteger = FirstPrime<NativeInteger>(kRes, 2 * n);

  NativeInteger q = PreviousPrime<NativeInteger>(firstInteger, 2 * n);
  moduli.push_back(q);
  roots_Of_Unity.push_back(RootOfUnity<NativeInteger>(2 * n, moduli[0]));

  NativeInteger prevQ = q;
  for (size_t i = 1; i < size; i++) {
    prevQ = lbcrypto::PreviousPrime<NativeInteger>(prevQ, 2 * n);
    NativeInteger nextRootOfUnity(RootOfUnity<NativeInteger>(2 * n, prevQ));
    moduli.push_back(prevQ);
    roots_Of_Unity.push_back(nextRootOfUnity);
  }

  auto ilDCRTParams =
      std::make_shared<ILDCRTParams<BigInteger>>(2 * n, moduli, roots_Of_Unity);

  ChineseRemainderTransformFTT<NativeVector>::PreCompute(roots_Of_Unity, 2 * n,
                                                         moduli);

  std::cout << "k: " << ilDCRTParams->GetModulus().GetMSB() << std::endl;

  size_t digitCount = (long)ceil(
      log2(ilDCRTParams->GetParams()[0]->GetModulus().ConvertToDouble()) /
      log2(base));
  size_t k = digitCount * ilDCRTParams->GetParams().size();

  std::cout << "digit count = " << digitCount << std::endl;
  std::cout << "k = " << k << std::endl;

  size_t m = k + 2;

  vector<LatticeSubgaussianUtility<NativeInteger>> util;

  for (size_t i = 0; i < moduli.size(); i++)
    util.push_back(
        LatticeSubgaussianUtility<NativeInteger>(base, moduli[i], digitCount));

  DCRTPoly::DggType dgg = DCRTPoly::DggType(SIGMA);
  DCRTPoly::DugType dug = DCRTPoly::DugType();
  DCRTPoly::BugType bug = DCRTPoly::BugType();

  auto zero_alloc = DCRTPoly::Allocator(ilDCRTParams, Format::EVALUATION);
  auto uniform_alloc =
      DCRTPoly::MakeDiscreteUniformAllocator(ilDCRTParams, Format::COEFFICIENT);
  auto gaussian_alloc = DCRTPoly::MakeDiscreteGaussianCoefficientAllocator(
      ilDCRTParams, Format::EVALUATION, SIGMA);

  Matrix<DCRTPoly> E1(zero_alloc, 1, m, gaussian_alloc);
  Matrix<DCRTPoly> E2(zero_alloc, 1, m, gaussian_alloc);
  Matrix<DCRTPoly> A(zero_alloc, 1, m, uniform_alloc);

  for (size_t i = 0; i < depth; i++) {
    auto gInverse = InverseRingVectorDCRT(util, A, 1);
    gInverse->SwitchFormat();
    E1 = E1 * (*gInverse);

    auto temp = E1;
    temp.SwitchFormat();
    std::cout << "level: " << i + 1 << "; norm: " << temp.Norm() << std::endl;
  }

  for (size_t i = 0; i < depth; i++) {
    auto gInverse = InverseG(A, base);
    // std::cerr << (*gInverse)(0,0).GetFormat() << std::endl;
    gInverse->SwitchFormat();
    E2 = E2 * (*gInverse);

    auto temp = E2;
    temp.SwitchFormat();
    std::cout << "level: " << i + 1 << "; norm: " << temp.Norm() << std::endl;
  }

  return 0;
}

shared_ptr<Matrix<DCRTPoly>> InverseG(const Matrix<DCRTPoly> &A, usint base) {
  usint n = A(0, 0).GetRingDimension();

  size_t k =
      (long)ceil(log2(A(0, 0).GetModulus().ConvertToDouble()) / log2(base));

  usint m = A.GetCols();

  auto zero_alloc_poly =
      DCRTPoly::Allocator(A(0, 0).GetParams(), Format::COEFFICIENT);
  auto psi = std::make_shared<Matrix<DCRTPoly>>(zero_alloc_poly, m, m);

  for (size_t i = 0; i < A.GetCols(); i++) {
    Poly temp = A(0, i).CRTInterpolate();
    for (size_t j = 0; j < n; j++) {
      std::vector<int64_t> digits = *(GetDigits(temp[j], base, k));
      for (size_t v = 0; v < A(0, 0).GetNumOfElements(); v++) {
        for (size_t p = 0; p < k; p++)
          (*psi)(p, i).ElementAtIndex(v)[j] = digits[p];
      }
    }
  }

  return psi;
}
