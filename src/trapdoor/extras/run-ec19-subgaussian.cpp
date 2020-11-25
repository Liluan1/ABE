// @file TODO
// @author TPOC: contact@palisade-crypto.org
//
// @copyright Copyright (c) TODO

#define PROFILE  // define this to enable PROFILELOG and TIC/TOC

// native libs
// #include <sys/resource.h>
#include <math.h>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>

#include "utils/debug.h"

#include "abe/kp_abe.h"
#include "subgaussian/subgaussian.h"

using namespace std;
using namespace lbcrypto;

void RunFigure1();
bool RunFigure2A();
bool RunFigure2B();

int main() {
  // RunFigure1();
  RunFigure2A();
  RunFigure2B();

  return 0;
}

void RunFigure1() {
  TimeVar t1;  // for TIC TOC

  double timeEval;

  usint n = 1024;
  size_t kRes = 60;

  NativeInteger firstInteger = FirstPrime<NativeInteger>(kRes, 2 * n);

  NativeInteger q = PreviousPrime<NativeInteger>(firstInteger, 2 * n);

  const size_t count = 100000;

  DiscreteUniformGeneratorImpl<NativeVector> dug;
  dug.SetModulus(q);

  NativeVector randomVector = dug.GenerateVector(count);

  for (size_t b = 2; b < 1073741825; b = b * 2) {
    size_t k = (long)ceil(log2(q.ConvertToDouble()) / log2(b));

    LatticeSubgaussianUtility<NativeInteger> sampler(b, q, k);

    vector<int64_t> nativeOutput(k);

    TIC(t1);  // start timer for total time
    for (size_t i = 0; i < count; i++)
      sampler.InverseG(randomVector[i], PseudoRandomNumberGenerator::GetPRNG(),
                       &nativeOutput);
    timeEval = TOC_US(t1);

    std::cout << "log2(base): " << log2(b)
              << "; Number of samples: " << 1e6 / timeEval * count << std::endl;
  }

  return;
}

bool RunFigure2A() {
  usint ringDimension = 4096;

  usint n = ringDimension * 2;  // cyclotomic order

  usint kRes = 60;

  usint base = 1 << 20;

  for (size_t j = 1; j < 11; j++) {
    size_t size = j;

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
      prevQ = PreviousPrime<NativeInteger>(prevQ, 2 * n);
      NativeInteger nextRootOfUnity(RootOfUnity<NativeInteger>(2 * n, prevQ));
      moduli.push_back(prevQ);
      roots_Of_Unity.push_back(nextRootOfUnity);
    }

    auto params =
        std::make_shared<ILDCRTParams<BigInteger>>(n, moduli, roots_Of_Unity);

    ChineseRemainderTransformFTT<NativeVector>::PreCompute(roots_Of_Unity,
                                                           2 * n, moduli);

    uint64_t digitCount = (long)ceil(log2(q.ConvertToDouble()) / log2(base));

    // std::cout << "digit count = " << size*digitCount << std::endl;

    usint m = moduli.size() * digitCount + 2;

    auto zero_alloc = DCRTPoly::Allocator(params, Format::COEFFICIENT);
    auto zero_alloc_eval = DCRTPoly::Allocator(params, Format::EVALUATION);

    Matrix<DCRTPoly> matrixTobeDecomposed(zero_alloc, 1, m);

    DiscreteGaussianGenerator dgg = DiscreteGaussianGenerator(SIGMA);
    DCRTPoly::DugType dug = DCRTPoly::DugType();

    for (usint i = 0; i < matrixTobeDecomposed.GetRows(); i++) {
      for (usint j = 0; j < matrixTobeDecomposed.GetCols(); j++) {
        matrixTobeDecomposed(i, j) = DCRTPoly(dug, params, Format::COEFFICIENT);
        // matrixTobeDecomposed(i, j).SwitchFormat(); // always kept in
        // Format::EVALUATION format
      }
    }

    Matrix<DCRTPoly> results(zero_alloc_eval, 1, m);
    Matrix<DCRTPoly> g(zero_alloc_eval, 1, m);

    size_t bk = 1;

    for (size_t k = 0; k < digitCount; k++) {
      for (size_t i = 0; i < moduli.size(); i++) {
        NativePoly temp(params->GetParams()[i]);
        temp = bk;
        g(0, k + i * digitCount).SetElementAtIndex(i, temp);
      }
      bk *= base;
    }

    // std::cout << g << std::endl;

    // std::cout << "g size = " << g.GetCols() << std::endl;

    TimeVar t1;  // for TIC TOC
    double timeEval(0.0);

    std::vector<LatticeSubgaussianUtility<NativeInteger>> sampler;

    for (size_t i = 0; i < size; i++)
      sampler.push_back(LatticeSubgaussianUtility<NativeInteger>(
          base, moduli[i], digitCount));

    TIC(t1);
    for (size_t k = 0; k < 10; k++) {
      auto psi = InverseRingVectorDCRT(sampler, matrixTobeDecomposed, 1);
    }
    timeEval += TOC_US(t1);

    auto psi = InverseRingVectorDCRT(sampler, matrixTobeDecomposed, 1);

    psi->SwitchFormat();
    results = g * (*psi);
    matrixTobeDecomposed.SwitchFormat();

    for (usint i = 0; i < results.GetRows(); i++) {
      for (usint j = 0; j < results.GetCols(); j++) {
        if (results(i, j) != matrixTobeDecomposed(i, j)) {
          std::cout << "index i = " << i << "; index j = " << j << std::endl;
          // std::cout<< results(i,j) <<std::endl;
          // std::cout<< matrixTobeDecomposed(i,j) <<std::endl;
          return false;
        }
      }
    }

    std::cout << "size =\t" << size
              << "\t; PALISADE polynomial sampling rate:\t"
              << 1e6 / timeEval * 10 * m << std::endl;
  }

  return true;
}

bool RunFigure2B() {
  usint ringDimension = 4096;

  usint n = ringDimension * 2;  // cyclotomic order

  usint kRes = 60;

  usint base = 1 << 20;

  for (size_t j = 1; j < 11; j++) {
    size_t size = j;

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
      prevQ = PreviousPrime<NativeInteger>(prevQ, 2 * n);
      NativeInteger nextRootOfUnity(RootOfUnity<NativeInteger>(2 * n, prevQ));
      moduli.push_back(prevQ);
      roots_Of_Unity.push_back(nextRootOfUnity);
    }

    auto params =
        std::make_shared<ILDCRTParams<BigInteger>>(n, moduli, roots_Of_Unity);
    auto ilParams =
        std::make_shared<ILParams>(n, params->GetModulus(), BigInteger(1));

    BigInteger bigModulus = params->GetModulus();

    int64_t digitCount =
        (long)ceil(log2(bigModulus.ConvertToDouble()) / log2(base));

    // std::cout << "digit count = " << digitCount << std::endl;

    usint m = digitCount + 2;

    auto zero_alloc = DCRTPoly::Allocator(params, Format::COEFFICIENT);
    auto zero_alloc_eval = DCRTPoly::Allocator(params, Format::EVALUATION);
    auto zero_alloc_poly = Poly::Allocator(ilParams, Format::COEFFICIENT);

    Matrix<DCRTPoly> matrixTobeDecomposed(zero_alloc, 1, m);

    DiscreteGaussianGenerator dgg = DiscreteGaussianGenerator(SIGMA);
    DCRTPoly::DugType dug = DCRTPoly::DugType();

    for (usint i = 0; i < matrixTobeDecomposed.GetRows(); i++) {
      for (usint j = 0; j < matrixTobeDecomposed.GetCols(); j++) {
        matrixTobeDecomposed(i, j) = DCRTPoly(dug, params, Format::COEFFICIENT);
        // matrixTobeDecomposed(i, j).SwitchFormat(); // always kept in
        // Format::EVALUATION format
      }
    }

    TimeVar t1;  // for TIC TOC
    double timeEval(0.0);

    Matrix<DCRTPoly> psiDCRT(zero_alloc, m, m);
    Matrix<DCRTPoly> psi(zero_alloc, m, m);

    LatticeSubgaussianUtility<BigInteger> sampler(base, bigModulus, digitCount);

    Matrix<Poly> matrixDecomposePoly(zero_alloc_poly, 1, m);

    TIC(t1);
    for (size_t k = 0; k < 10; k++) {
      for (usint i = 0; i < m; i++) {
        matrixDecomposePoly(0, i) = matrixTobeDecomposed(0, i).CRTInterpolate();
      }

      InverseRingVectorSpecial<Poly>(sampler, ilParams, matrixDecomposePoly, 1,
                                     &psi);
    }
    timeEval += TOC_US(t1);

    std::cout << "size =\t" << size
              << "\t; PALISADE polynomial sampling rate:\t"
              << 1e6 / timeEval * 10 * m << std::endl;
  }

  return true;
}
