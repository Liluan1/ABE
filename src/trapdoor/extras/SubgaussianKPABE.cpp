// @file TODO
// @author TPOC: contact@palisade-crypto.org
//
// @copyright Copyright (c) TODO

#define PROFILE  // define this to enable PROFILELOG and TIC/TOC

#include <fstream>
#include <iostream>
#include "abe/kp_abe.h"
#include "subgaussian/subgaussian.h"

#include "utils/debug.h"

using namespace lbcrypto;

void KPABEBenchMarkCircuit(int32_t base, usint k, usint ringDimension,
                           usint iter);
template <class Element, class Element2>
bool TestDCRTVecDecompose(int32_t base, usint k, usint ringDimension,
                          size_t type);
int KPABE_BenchmarkCircuitTestDCRT(usint iter, int32_t base);
usint EvalNANDTree(usint *x, usint ell);
void KPABE_NANDGATE(int32_t base, usint k, usint ringDimension);
void KPABE_NANDGATEDCRT(int32_t base, usint k, usint ringDimension);
void KPABEANDGate(int32_t base, usint k, usint ringDimension);
void KPABEANDGateDCRT(int32_t base, usint k, usint ringDimension);

void KPABE_NANDGATE_RANDOM(int32_t base, usint k, usint ringDimension);
void KPABE_NANDGATEDCRT_RANDOM(int32_t base, usint k, usint ringDimension);

bool TestDCRTVecDecomposeOptimized(int32_t base, usint k, usint ringDimension);

int main() {
  if (TestDCRTVecDecomposeOptimized(16, 51, 1024))
    std::cout << "CRT decomposition test: SUCCESS" << std::endl;
  else
    std::cout << "CRT digit decomposition test: FAILURE" << std::endl;
  /*
          //KPABE_BenchmarkCircuitTestDCRT(4, 32);
          if (TestDCRTVecDecompose<DCRTPoly, Poly>(256,51,2048,1))
                  std::cout << "NAF digit decomposition test: SUCCESS" <<
     std::endl; else std::cout << "NAF digit decomposition test: FAILURE" <<
     std::endl;

          if (TestDCRTVecDecompose<DCRTPoly, Poly>(257,51,2048,2))
                  std::cout << "Randomized digit decomposition test: SUCCESS" <<
     std::endl; else std::cout << "Randomized digit decomposition test: FAILURE"
     << std::endl;

          //KPABEBenchMarkCircuit(2, 51, 2048, 100);

          std::cout << "\nNAF Poly test" << std::endl;

          KPABE_NANDGATE(32,51,2048);

          std::cout << "\nSUBGAUSSIAN Poly test" << std::endl;

          KPABE_NANDGATE_RANDOM(33,51,2048);

          std::cout << "\nNAF DCRTPoly test" << std::endl;

          KPABE_NANDGATEDCRT(16, 8, 2048);

          std::cout << "\nSUBGAUSSIAN DCRTPoly test" << std::endl;

          KPABE_NANDGATEDCRT_RANDOM(17, 8, 2048);
  */
  // KPABEANDGate(32,51,2048);
  // KPABEANDGateDCRT(16, 8, 2048);
  return 0;
}

void KPABEBenchMarkCircuit(int32_t base, usint k, usint ringDimension,
                           usint iter) {
  usint n = ringDimension * 2;  // cyclotomic order
  usint ell = 4;                // No of attributes

  BigInteger q = BigInteger(1) << (k - 1);
  q = lbcrypto::FirstPrime<BigInteger>(k, n);
  BigInteger rootOfUnity(RootOfUnity(n, q));

  double val = q.ConvertToDouble();
  double logTwo = log(val - 1.0) / log(base) + 1.0;
  size_t k_ = (usint)floor(logTwo) + 1;

  usint m = k_ + 2;

  auto ilParams = std::make_shared<ILParams>(n, q, rootOfUnity);

  auto zero_alloc = Poly::Allocator(ilParams, Format::COEFFICIENT);

  DiscreteGaussianGenerator dgg = DiscreteGaussianGenerator(SIGMA);
  Poly::DugType dug = Poly::DugType();
  dug.SetModulus(q);
  BinaryUniformGenerator bug = BinaryUniformGenerator();

  // Precompuations for FTT
  ChineseRemainderTransformFTT<BigVector>::PreCompute(rootOfUnity, n, q);

  // Trapdoor Generation
  std::pair<Matrix<Poly>, RLWETrapdoorPair<Poly>> trapdoorA =
      RLWETrapdoorUtility<Poly>::TrapdoorGen(
          ilParams, SIGMA, base, true);  // A.first is the public element

  Poly pubElemBeta(dug, ilParams, Format::EVALUATION);

  Matrix<Poly> publicElementB(zero_alloc, ell + 1, m);
  Matrix<Poly> ctCin(zero_alloc, ell + 2, m);
  Poly c1(dug, ilParams, Format::EVALUATION);

  KPABE<Poly, Poly> pkg, sender, receiver;

  pkg.Setup(ilParams, base, ell, dug, &publicElementB);
  sender.Setup(ilParams, base, ell);
  receiver.Setup(ilParams, base, ell);

  usint x[] = {1, 1,
               1};  // array of attributes, everything is set to 1 for NAND gate
                    // evaluation, values set based on experimental results

  usint y;

  // plaintext
  Poly ptext(ilParams, Format::COEFFICIENT, true);

  // circuit outputs
  Matrix<Poly> evalBf(Poly::Allocator(ilParams, Format::EVALUATION), 1,
                      m);  // evaluated Bs
  Matrix<Poly> evalCf(Poly::Allocator(ilParams, Format::EVALUATION), 1,
                      m);  // evaluated Cs
  Matrix<Poly> ctCA(Poly::Allocator(ilParams, Format::EVALUATION), 1, m);  // CA

  for (usint i = 0; i < iter; i++) {
    // secret key corresponding to the circuit output
    Matrix<Poly> sk(zero_alloc, 2, m);

    // decrypted text
    Poly dtext(ilParams, Format::EVALUATION, true);
    // Encrypt a uniformly randomly selected message ptext (in ptext in $R_2$)
    ptext.SetValues(bug.GenerateVector(ringDimension, q), Format::COEFFICIENT);
    ptext.SwitchFormat();
    sender.Encrypt(ilParams, trapdoorA.first, publicElementB, pubElemBeta,
                   &x[0], ptext, dgg, dug, bug, &ctCin,
                   &c1);  // Cin and c1 are the ciphertext

    ctCA = ctCin.ExtractRow(0);  // CA is A^T * s + e 0,A

    receiver.EvalCT(ilParams, publicElementB, &x[0],
                    ctCin.ExtractRows(1, ell + 1), &y, &evalCf);

    pkg.EvalPK(ilParams, publicElementB, &evalBf);
    pkg.KeyGen(ilParams, trapdoorA.first, evalBf, pubElemBeta, trapdoorA.second,
               dgg, &sk);

    Poly t(pubElemBeta);
    t.SetValuesToZero();

    for (usint i = 0; i < m; i++) {
      t += (trapdoorA.first(0, i) * sk(0, i));
      t += (evalBf(0, i) * sk(1, i));
    }

    receiver.Decrypt(ilParams, sk, ctCA, evalCf, c1, &dtext);
    receiver.Decode(&dtext);

    ptext.SwitchFormat();

    if (ptext.GetValues() == dtext.GetValues()) {
      std::cout << "Decrypted Properly" << std::endl;
    }
  }
}

int KPABE_BenchmarkCircuitTestDCRT(usint iter, int32_t base) {
  usint ringDimension = 2048;   // ring dimension
  usint n = ringDimension * 2;  // cyclotomic order
                                //  usint k = 21;
  usint ell = 4;                // No of attributes

  //  NativeInteger q = NativeInteger::ONE << (k - 1);
  //  q = lbcrypto::FirstPrime<NativeInteger>(k, n);

  NativeInteger q("2101249");

  NativeInteger rootOfUnity(RootOfUnity(n, q));

  double val = q.ConvertToDouble();
  double logTwo = log(val - 1.0) / log(base) + 1.0;
  size_t k_ = (usint)floor(logTwo) + 1;  //  (+1) is For NAF
  std::cout << "q: " << q << std::endl;
  std::cout << "modulus length: " << k_ << std::endl;
  std::cout << "root of unity: " << rootOfUnity << std::endl;
  std::cout << "Standard deviation: " << SIGMA << std::endl;

  //  NativeInteger nextQ = NativeInteger::ONE << (k-1);
  //  nextQ = lbcrypto::NextPrime<NativeInteger>(q, n);
  //  std::cout << "nextQ: " << nextQ << std::endl;

  NativeInteger nextQ("2236417");

  NativeInteger nextRootOfUnity(RootOfUnity<NativeInteger>(n, nextQ));

  //  NativeInteger nextQ2 = NativeInteger::ONE << (k-1);
  //  nextQ2 = lbcrypto::NextPrime<NativeInteger>(nextQ, n);

  NativeInteger nextQ2("2277377");
  NativeInteger nextRootOfUnity2(RootOfUnity<NativeInteger>(n, nextQ2));

  usint m = 3 * k_ + 2;

  std::vector<NativeInteger> moduli;
  std::vector<NativeInteger> roots_Of_Unity;
  moduli.reserve(3);
  roots_Of_Unity.reserve(3);

  moduli.push_back(q);
  moduli.push_back(nextQ);
  moduli.push_back(nextQ2);

  roots_Of_Unity.push_back(rootOfUnity);
  roots_Of_Unity.push_back(nextRootOfUnity);
  roots_Of_Unity.push_back(nextRootOfUnity2);

  BigInteger bigModulus =
      BigInteger("2101249") * BigInteger("2236417") * BigInteger("2277377");

  BigInteger bigRootOfUnity(RootOfUnity(n, bigModulus));

  BinaryUniformGenerator bug = BinaryUniformGenerator();

  auto ilDCRTParams =
      std::make_shared<ILDCRTParams<BigInteger>>(n, moduli, roots_Of_Unity);
  auto ilParamsConsolidated =
      std::make_shared<ILParams>(n, bigModulus, bigRootOfUnity);

  auto zero_alloc = DCRTPoly::Allocator(ilDCRTParams, Format::COEFFICIENT);

  DCRTPoly::DggType dgg = DCRTPoly::DggType(SIGMA);
  DCRTPoly::DugType dug = DCRTPoly::DugType();

  // Precompuations for FTT
  //  ChineseRemainderTransformFTT<NativeInteger,
  // BigVector>::PreCompute(rootOfUnity, n, q);

  // Trapdoor Generation
  std::pair<Matrix<DCRTPoly>, RLWETrapdoorPair<DCRTPoly>> trapdoorA =
      RLWETrapdoorUtility<DCRTPoly>::TrapdoorGen(
          ilDCRTParams, SIGMA, base, true);  // A.first is the public element

  DCRTPoly pubElemBeta(dug, ilDCRTParams, Format::EVALUATION);

  Matrix<DCRTPoly> publicElementB(zero_alloc, ell + 1, m);
  Matrix<DCRTPoly> ctCin(zero_alloc, ell + 2, m);
  DCRTPoly c1(dug, ilDCRTParams, Format::EVALUATION);

  KPABE<DCRTPoly, Poly> pkg, sender, receiver;

  pkg.Setup(ilDCRTParams, base, ell, dug, &publicElementB);
  sender.Setup(ilDCRTParams, base, ell);
  receiver.Setup(ilDCRTParams, base, ell);

  // Attribute values all are set to 1 for NAND gate evaluation
  std::vector<usint> x(ell + 1);

  usint found = 0;
  while (found == 0) {
    for (usint i = 1; i < ell + 1; i++) x[i] = rand() & 0x1;
    if (EvalNANDTree(&x[1], ell) == 0) found = 1;
  }

  usint y;

  double avg_keygen, avg_eval, avg_enc, avg_dec;
  avg_keygen = avg_eval = avg_enc = avg_dec = 0.0;

  // plaintext
  for (usint i = 0; i < iter; i++) {
    Poly ptext1(ilParamsConsolidated, Format::COEFFICIENT, true);
    ptext1.SetValues(bug.GenerateVector(ringDimension, bigModulus),
                     Format::COEFFICIENT);

    DCRTPoly ptext(ptext1, ilDCRTParams);

    // circuit outputs
    Matrix<DCRTPoly> evalBf(
        DCRTPoly::Allocator(ilDCRTParams, Format::EVALUATION), 1,
        m);  // evaluated Bs
    Matrix<DCRTPoly> evalCf(
        DCRTPoly::Allocator(ilDCRTParams, Format::EVALUATION), 1,
        m);  // evaluated Cs
    Matrix<DCRTPoly> ctCA(DCRTPoly::Allocator(ilDCRTParams, Format::EVALUATION),
                          1, m);  // CA

    // secret key corresponding to the circuit output
    Matrix<DCRTPoly> sk(zero_alloc, 2, m);

    // decrypted text
    DCRTPoly dtext(ilDCRTParams, Format::EVALUATION, true);

    double start, finish;

    ptext.SwitchFormat();
    start = currentDateTime();
    sender.Encrypt(ilDCRTParams, trapdoorA.first, publicElementB, pubElemBeta,
                   &x[0], ptext, dgg, dug, bug, &ctCin,
                   &c1);  // Cin and c1 are the ciphertext
    finish = currentDateTime();
    avg_enc += (finish - start);

    ctCA = ctCin.ExtractRow(0);  // CA is A^T * s + e 0,A

    start = currentDateTime();
    receiver.EvalCTDCRT(ilDCRTParams, publicElementB, &x[0],
                        ctCin.ExtractRows(1, ell + 1), &y, &evalCf,
                        ilParamsConsolidated);

    finish = currentDateTime();
    avg_eval += (finish - start);

    start = currentDateTime();
    pkg.EvalPKDCRT(ilDCRTParams, publicElementB, &evalBf, ilParamsConsolidated);
    pkg.KeyGen(ilDCRTParams, trapdoorA.first, evalBf, pubElemBeta,
               trapdoorA.second, dgg, &sk);
    finish = currentDateTime();
    avg_keygen += (finish - start);
    //  CheckSecretKeyKPDCRT(m, trapdoorA.first, evalBf, sk, pubElemBeta);

    start = currentDateTime();
    receiver.Decrypt(ilDCRTParams, sk, ctCA, evalCf, c1, &dtext);

    Poly dtextPoly(dtext.CRTInterpolate());

    receiver.Decode(&dtextPoly);

    finish = currentDateTime();
    avg_dec += (finish - start);

    if (ptext1.GetValues() != dtextPoly.GetValues()) {
      std::cout << "Decryption fails at iteration: " << i << std::endl;
      return 0;
    }
  }

  std::cout << "Encryption is successful after " << iter << " iterations!\n";
  std::cout << "Average key generation time : "
            << "\t" << (avg_keygen) / iter << " ms" << std::endl;
  std::cout << "Average evaluation time : "
            << "\t" << (avg_eval) / iter << " ms" << std::endl;
  std::cout << "Average encryption time : "
            << "\t" << (avg_enc) / iter << " ms" << std::endl;
  std::cout << "Average decryption time : "
            << "\t" << (avg_dec) / iter << " ms" << std::endl;

  return 0;
}

template <class Element, class Element2>
bool TestDCRTVecDecompose(int32_t base, usint k, usint ringDimension,
                          size_t type) {
  usint n = ringDimension * 2;  // cyclotomic order

  NativeInteger q = NativeInteger(1) << (k - 1);
  q = lbcrypto::FirstPrime<NativeInteger>(k, n);
  NativeInteger rootOfUnity(RootOfUnity<NativeInteger>(n, q));

  NativeInteger nextQ = NativeInteger(1) << (k - 1);
  nextQ = lbcrypto::NextPrime<NativeInteger>(q, n);
  NativeInteger nextRootOfUnity(RootOfUnity<NativeInteger>(n, nextQ));

  std::vector<NativeInteger> moduli;
  std::vector<NativeInteger> roots_Of_Unity;
  moduli.reserve(2);
  roots_Of_Unity.reserve(2);

  moduli.push_back(q);
  moduli.push_back(nextQ);

  roots_Of_Unity.push_back(rootOfUnity);
  roots_Of_Unity.push_back(nextRootOfUnity);

  BigInteger bigModulus = BigInteger(q) * BigInteger(nextQ);

  BigInteger bigRootOfUnity(RootOfUnity(n, bigModulus));

  auto params =
      std::make_shared<ILDCRTParams<BigInteger>>(n, moduli, roots_Of_Unity);
  auto ilParams = std::make_shared<ILParams>(n, bigModulus, rootOfUnity);

  int64_t digitCount =
      (long)ceil(log2(bigModulus.ConvertToDouble()) / log2(base));

  std::cout << "digit count = " << digitCount << std::endl;

  usint m = digitCount + 2;

  auto zero_alloc_poly = Element2::Allocator(ilParams, Format::COEFFICIENT);
  auto zero_alloc = Element::Allocator(params, Format::COEFFICIENT);
  auto zero_alloc_eval = DCRTPoly::Allocator(params, Format::EVALUATION);

  Matrix<DCRTPoly> matrixTobeDecomposed(zero_alloc, 1, m);

  DiscreteGaussianGenerator dgg = DiscreteGaussianGenerator(SIGMA);
  DCRTPoly::DugType dug = DCRTPoly::DugType();

  for (usint i = 0; i < matrixTobeDecomposed.GetRows(); i++) {
    for (usint j = 0; j < matrixTobeDecomposed.GetCols(); j++) {
      matrixTobeDecomposed(i, j) = Element(dug, params, Format::COEFFICIENT);
      matrixTobeDecomposed(i, j)
          .SwitchFormat();  // always kept in Format::EVALUATION format
    }
  }

  Matrix<DCRTPoly> results(zero_alloc_eval, 1, m);
  Matrix<DCRTPoly> g =
      Matrix<DCRTPoly>(zero_alloc_eval, 1, m).GadgetVector(base);

  std::cout << "g size = " << g.GetCols() << std::endl;

  Matrix<DCRTPoly> psiDCRT(zero_alloc, m, m);
  Matrix<Poly> psi(zero_alloc_poly, m, m);

  Matrix<Poly> matrixDecomposePoly(zero_alloc_poly, 1, m);

  for (usint i = 0; i < m; i++) {
    matrixDecomposePoly(0, i) = matrixTobeDecomposed(0, i).CRTInterpolate();
  }

  TimeVar t1;  // for TIC TOC
  double timeEval;

  if (type == 1) {
    TIC(t1);
    lbcrypto::PolyVec2BalDecom<Poly>(ilParams, base, digitCount,
                                     matrixDecomposePoly, &psi);
    timeEval = TOC(t1);
  } else {
    int64_t digitCount =
        (long)ceil(log2(bigModulus.ConvertToDouble()) / log2(base));
    LatticeSubgaussianUtility<BigInteger> sampler(base, bigModulus, digitCount);
    TIC(t1);
    InverseRingVector<Poly>(sampler, ilParams, matrixDecomposePoly, 1, &psi);
    timeEval = TOC(t1);
  }

  for (usint i = 0; i < psi.GetRows(); i++) {
    for (usint j = 0; j < psi.GetCols(); j++) {
      Element temp(psi(i, j), params);
      psiDCRT(i, j) = temp;
    }
  }
  psiDCRT.SwitchFormat();
  results = g * psiDCRT;

  for (usint i = 0; i < results.GetRows(); i++) {
    for (usint j = 0; j < results.GetCols(); j++) {
      if (results(i, j) != matrixTobeDecomposed(i, j)) {
        std::cout << "index i = " << i << "; index j = " << j << std::endl;
        std::cout << results(i, j) << std::endl;
        std::cout << matrixTobeDecomposed(i, j) << std::endl;
        return false;
      }
    }
  }

  std::cout << "PALISADE digit decomposition time: " << timeEval
            << " microseconds" << std::endl;

  return true;
}

bool TestDCRTVecDecomposeOptimized(int32_t base, usint k, usint ringDimension) {
  usint n = ringDimension * 2;  // cyclotomic order

  size_t size = 4;

  std::vector<NativeInteger> moduli;
  std::vector<NativeInteger> roots_Of_Unity;

  NativeInteger q = NativeInteger(1) << (k - 1);
  q = lbcrypto::FirstPrime<NativeInteger>(k, n);
  NativeInteger rootOfUnity(RootOfUnity<NativeInteger>(n, q));
  moduli.push_back(q);
  roots_Of_Unity.push_back(rootOfUnity);

  NativeInteger nextQ = q;
  for (size_t i = 1; i < size; i++) {
    nextQ = lbcrypto::NextPrime<NativeInteger>(nextQ, n);
    NativeInteger nextRootOfUnity(RootOfUnity<NativeInteger>(n, nextQ));
    moduli.push_back(nextQ);
    roots_Of_Unity.push_back(nextRootOfUnity);
  }

  auto params =
      std::make_shared<ILDCRTParams<BigInteger>>(n, moduli, roots_Of_Unity);

  uint64_t digitCount = (long)ceil(log2(q.ConvertToDouble()) / log2(base));

  std::cout << "digit count = " << size * digitCount << std::endl;

  usint m = moduli.size() * digitCount + 2;

  auto zero_alloc = DCRTPoly::Allocator(params, Format::COEFFICIENT);
  auto zero_alloc_eval = DCRTPoly::Allocator(params, Format::EVALUATION);

  Matrix<DCRTPoly> matrixTobeDecomposed(zero_alloc, 1, m);

  DiscreteGaussianGenerator dgg = DiscreteGaussianGenerator(SIGMA);
  DCRTPoly::DugType dug = DCRTPoly::DugType();

  for (usint i = 0; i < matrixTobeDecomposed.GetRows(); i++) {
    for (usint j = 0; j < matrixTobeDecomposed.GetCols(); j++) {
      matrixTobeDecomposed(i, j) = DCRTPoly(dug, params, Format::COEFFICIENT);
      // always kept in Format::EVALUATION format
      matrixTobeDecomposed(i, j).SwitchFormat();
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

  std::cout << "g size = " << g.GetCols() << std::endl;

  TimeVar t1;  // for TIC TOC
  double timeEval;

  std::vector<LatticeSubgaussianUtility<NativeInteger>> sampler;

  for (size_t i = 0; i < size; i++)
    sampler.push_back(
        LatticeSubgaussianUtility<NativeInteger>(base, moduli[i], digitCount));

  TIC(t1);
  auto psi = InverseRingVectorDCRT(sampler, matrixTobeDecomposed, 1);
  timeEval = TOC(t1);

  psi->SwitchFormat();
  results = g * (*psi);

  for (usint i = 0; i < results.GetRows(); i++) {
    for (usint j = 0; j < results.GetCols(); j++) {
      if (results(i, j) != matrixTobeDecomposed(i, j)) {
        std::cout << "index i = " << i << "; index j = " << j << std::endl;
        std::cout << results(i, j) << std::endl;
        std::cout << matrixTobeDecomposed(i, j) << std::endl;
        return false;
      }
    }
  }

  std::cout << "PALISADE digit decomposition time: " << timeEval
            << " microseconds" << std::endl;

  return true;
}

usint EvalNANDTree(usint *x, usint ell) {
  usint y;

  if (ell == 2) {
    y = 1 - x[0] * x[1];
    return y;
  } else {
    ell >>= 1;
    y = 1 - (EvalNANDTree(&x[0], ell) * EvalNANDTree(&x[ell], ell));
  }
  return y;
}

void KPABE_NANDGATE(int32_t base, usint k, usint ringDimension) {
  usint n = ringDimension * 2;
  usint ell = 2;  // No of attributes for NAND gate

  BigInteger q = BigInteger(1) << (k - 1);
  q = lbcrypto::FirstPrime<BigInteger>(k, n);
  BigInteger rootOfUnity(RootOfUnity(n, q));

  double val = q.ConvertToDouble();
  double logTwo = log(val - 1.0) / log(base) + 1.0;
  size_t k_ = (usint)floor(logTwo) + 1;

  usint m = k_ + 2;

  auto ilParams = std::make_shared<ILParams>(n, q, rootOfUnity);

  auto zero_alloc = Poly::Allocator(ilParams, Format::COEFFICIENT);

  DiscreteGaussianGenerator dgg = DiscreteGaussianGenerator(SIGMA);
  Poly::DugType dug = Poly::DugType();
  dug.SetModulus(q);
  BinaryUniformGenerator bug = BinaryUniformGenerator();

  // Precompuations for FTT
  ChineseRemainderTransformFTT<BigVector>::PreCompute(rootOfUnity, n, q);

  // Trapdoor Generation
  std::pair<Matrix<Poly>, RLWETrapdoorPair<Poly>> A =
      RLWETrapdoorUtility<Poly>::TrapdoorGen(ilParams, SIGMA, base, true);

  Poly pubElemBeta(dug, ilParams, Format::EVALUATION);

  Matrix<Poly> publicElementB(zero_alloc, ell + 1, m);
  Matrix<Poly> ctCin(zero_alloc, ell + 2, m);
  Poly c1(dug, ilParams, Format::EVALUATION);

  KPABE<Poly, Poly> pkg, sender, receiver;

  pkg.Setup(ilParams, base, ell, dug, &publicElementB);
  sender.Setup(ilParams, base, ell);
  receiver.Setup(ilParams, base, ell);

  // Attribute values all are set to 1 for NAND gate evaluation
  std::vector<usint> x(ell + 1);
  x[0] = x[1] = x[2] = 1;
  usint y;

  // plain text in $R_2$
  Poly ptext(ilParams, Format::COEFFICIENT, true);

  // circuit outputs
  Matrix<Poly> pubElemBf(Poly::Allocator(ilParams, Format::EVALUATION), 1, m);
  Matrix<Poly> ctCf(Poly::Allocator(ilParams, Format::EVALUATION), 1, m);
  Matrix<Poly> ctCA(Poly::Allocator(ilParams, Format::EVALUATION), 1, m);

  // Secret key for the output of the circuit
  Matrix<Poly> sk(zero_alloc, 2, m);

  // text after the decryption
  Poly dtext(ilParams, Format::EVALUATION, true);

  // Encrypt a uniformly randomly selected message ptext (in ptext in $R_2$)
  ptext.SetValues(bug.GenerateVector(ringDimension, q), Format::COEFFICIENT);
  ptext.SwitchFormat();

  sender.Encrypt(ilParams, A.first, publicElementB, pubElemBeta, &x[0], ptext,
                 dgg, dug, bug, &ctCin, &c1);

  ctCA = ctCin.ExtractRow(0);

  receiver.KPABE::NANDGateEvalPK(ilParams, publicElementB.ExtractRow(0),
                                 publicElementB.ExtractRows(1, 2), &pubElemBf);

  receiver.KPABE::NANDGateEvalCT(ilParams, ctCin.ExtractRow(1), &x[1],
                                 publicElementB.ExtractRows(1, 2),
                                 ctCin.ExtractRows(2, 3), &y, &ctCf);

  pkg.KeyGen(ilParams, A.first, pubElemBf, pubElemBeta, A.second, dgg, &sk);

  receiver.Decrypt(ilParams, sk, ctCA, ctCf, c1, &dtext);
  receiver.Decode(&dtext);

  ptext.SwitchFormat();
  if (ptext.GetValues() == dtext.GetValues()) {
    std::cout << "Success" << std::endl;
  }
}

void KPABE_NANDGATE_RANDOM(int32_t base, usint k, usint ringDimension) {
  usint n = ringDimension * 2;
  usint ell = 2;  // No of attributes for NAND gate

  BigInteger q = BigInteger(1) << (k - 1);
  q = lbcrypto::FirstPrime<BigInteger>(k, n);
  BigInteger rootOfUnity(RootOfUnity(n, q));

  // double val = q.ConvertToDouble();
  // double logTwo = log(val-1.0)/log(base)+1.0;

  size_t k_ = (long)ceil(log2(q.ConvertToDouble()) / log2(base));

  usint m = k_ + 2;

  auto ilParams = std::make_shared<ILParams>(n, q, rootOfUnity);

  auto zero_alloc = Poly::Allocator(ilParams, Format::COEFFICIENT);

  DiscreteGaussianGenerator dgg = DiscreteGaussianGenerator(SIGMA);
  Poly::DugType dug = Poly::DugType();
  dug.SetModulus(q);
  BinaryUniformGenerator bug = BinaryUniformGenerator();

  // Precompuations for FTT
  ChineseRemainderTransformFTT<BigVector>::PreCompute(rootOfUnity, n, q);

  // Trapdoor Generation
  std::pair<Matrix<Poly>, RLWETrapdoorPair<Poly>> A =
      RLWETrapdoorUtility<Poly>::TrapdoorGen(ilParams, SIGMA, base, false);

  Poly pubElemBeta(dug, ilParams, Format::EVALUATION);

  Matrix<Poly> publicElementB(zero_alloc, ell + 1, m);
  Matrix<Poly> ctCin(zero_alloc, ell + 2, m);
  Poly c1(dug, ilParams, Format::EVALUATION);

  KPABE<Poly, Poly> pkg(SUBGAUSSIAN), sender(SUBGAUSSIAN),
      receiver(SUBGAUSSIAN);

  pkg.Setup(ilParams, base, ell, dug, &publicElementB);
  sender.Setup(ilParams, base, ell);
  receiver.Setup(ilParams, base, ell);

  // Attribute values all are set to 1 for NAND gate evaluation
  std::vector<usint> x(ell + 1);
  x[0] = x[1] = x[2] = 1;
  usint y;

  // plain text in $R_2$
  Poly ptext(ilParams, Format::COEFFICIENT, true);

  // circuit outputs
  Matrix<Poly> pubElemBf(Poly::Allocator(ilParams, Format::EVALUATION), 1, m);
  Matrix<Poly> ctCf(Poly::Allocator(ilParams, Format::EVALUATION), 1, m);
  Matrix<Poly> ctCA(Poly::Allocator(ilParams, Format::EVALUATION), 1, m);

  // Secret key for the output of the circuit
  Matrix<Poly> sk(zero_alloc, 2, m);

  // text after the decryption
  Poly dtext(ilParams, Format::EVALUATION, true);

  // Encrypt a uniformly randomly selected message ptext (in ptext in $R_2$)
  ptext.SetValues(bug.GenerateVector(ringDimension, q), Format::COEFFICIENT);
  ptext.SwitchFormat();

  sender.Encrypt(ilParams, A.first, publicElementB, pubElemBeta, &x[0], ptext,
                 dgg, dug, bug, &ctCin, &c1);

  ctCA = ctCin.ExtractRow(0);

  receiver.KPABE::NANDGateEvalPK(ilParams, publicElementB.ExtractRow(0),
                                 publicElementB.ExtractRows(1, 2), &pubElemBf,
                                 1);

  receiver.KPABE::NANDGateEvalCT(ilParams, ctCin.ExtractRow(1), &x[1],
                                 publicElementB.ExtractRows(1, 2),
                                 ctCin.ExtractRows(2, 3), &y, &ctCf, 1);

  pkg.KeyGen(ilParams, A.first, pubElemBf, pubElemBeta, A.second, dgg, &sk);

  receiver.Decrypt(ilParams, sk, ctCA, ctCf, c1, &dtext);

  receiver.Decode(&dtext);

  ptext.SwitchFormat();
  if (ptext.GetValues() == dtext.GetValues()) {
    std::cout << "Success" << std::endl;
  } else {
    std::cout << "Failure" << std::endl;
  }
}

void KPABE_NANDGATEDCRT(int32_t base, usint k, usint ringDimension) {
  usint n = ringDimension * 2;  // cyclotomic order

  //  usint k = 21;
  usint ell = 4;  // No of attributes

  //  NativeInteger q = NativeInteger::ONE << (k - 1);
  //  q = lbcrypto::FirstPrime<NativeInteger>(k, n);

  NativeInteger q("2101249");

  NativeInteger rootOfUnity(RootOfUnity(n, q));

  //  NativeInteger rootOfUnity("794438271477401");

  double val = q.ConvertToDouble();
  double logTwo = log(val - 1.0) / log(base) + 1.0;
  size_t k_ = (usint)floor(logTwo) + 1;  //  (+1) is For NAF
  std::cout << "q: " << q << std::endl;
  std::cout << "modulus length: " << k_ << std::endl;
  std::cout << "root of unity: " << rootOfUnity << std::endl;
  std::cout << "Standard deviation: " << SIGMA << std::endl;

  //  NativeInteger nextQ = NativeInteger::ONE << (k-1);
  //  nextQ = lbcrypto::NextPrime<NativeInteger>(q, n);
  //  std::cout << "nextQ: " << nextQ << std::endl;

  NativeInteger nextQ("2236417");

  NativeInteger nextRootOfUnity(RootOfUnity<NativeInteger>(n, nextQ));

  //  NativeInteger nextQ2 = NativeInteger::ONE << (k-1);
  //  nextQ2 = lbcrypto::NextPrime<NativeInteger>(nextQ, n);

  NativeInteger nextQ2("2277377");
  NativeInteger nextRootOfUnity2(RootOfUnity<NativeInteger>(n, nextQ2));

  usint m = 3 * k_ + 2;

  std::vector<NativeInteger> moduli;
  std::vector<NativeInteger> roots_Of_Unity;
  moduli.reserve(3);
  roots_Of_Unity.reserve(3);

  moduli.push_back(q);
  moduli.push_back(nextQ);
  moduli.push_back(nextQ2);

  roots_Of_Unity.push_back(rootOfUnity);
  roots_Of_Unity.push_back(nextRootOfUnity);
  roots_Of_Unity.push_back(nextRootOfUnity2);

  BigInteger bigModulus =
      BigInteger("2101249") * BigInteger("2236417") * BigInteger("2277377");

  BigInteger bigRootOfUnity(RootOfUnity(n, bigModulus));

  BinaryUniformGenerator bug = BinaryUniformGenerator();

  auto ilDCRTParams =
      std::make_shared<ILDCRTParams<BigInteger>>(n, moduli, roots_Of_Unity);
  auto ilParamsConsolidated =
      std::make_shared<ILParams>(n, bigModulus, bigRootOfUnity);

  auto zero_alloc = DCRTPoly::Allocator(ilDCRTParams, Format::COEFFICIENT);

  DCRTPoly::DggType dgg = DCRTPoly::DggType(SIGMA);
  DCRTPoly::DugType dug = DCRTPoly::DugType();

  // Trapdoor Generation
  std::pair<Matrix<DCRTPoly>, RLWETrapdoorPair<DCRTPoly>> trapdoorA =
      RLWETrapdoorUtility<DCRTPoly>::TrapdoorGen(
          ilDCRTParams, SIGMA, base, true);  // A.first is the public element

  DCRTPoly pubElemBeta(dug, ilDCRTParams, Format::EVALUATION);

  Matrix<DCRTPoly> publicElementB(zero_alloc, ell + 1, m);
  Matrix<DCRTPoly> ctCin(zero_alloc, ell + 2, m);
  DCRTPoly c1(dug, ilDCRTParams, Format::EVALUATION);

  KPABE<DCRTPoly, Poly> pkg, sender, receiver;

  pkg.Setup(ilDCRTParams, base, ell, dug, &publicElementB);
  sender.Setup(ilDCRTParams, base, ell);
  receiver.Setup(ilDCRTParams, base, ell);

  // Attribute values all are set to 1 for NAND gate evaluation
  std::vector<usint> x(ell + 1);
  x[0] = x[1] = x[2] = 1;
  usint y;

  Poly ptext1(ilParamsConsolidated, Format::COEFFICIENT, true);
  ptext1.SetValues(bug.GenerateVector(ringDimension, bigModulus),
                   Format::COEFFICIENT);

  DCRTPoly ptext(ptext1, ilDCRTParams);

  // circuit outputs
  // evaluated Bs
  Matrix<DCRTPoly> pubElemBf(
      DCRTPoly::Allocator(ilDCRTParams, Format::EVALUATION), 1, m);
  // evaluated Cs
  Matrix<DCRTPoly> ctCf(DCRTPoly::Allocator(ilDCRTParams, Format::EVALUATION),
                        1, m);
  // CA
  Matrix<DCRTPoly> ctCA(DCRTPoly::Allocator(ilDCRTParams, Format::EVALUATION),
                        1, m);

  // secret key corresponding to the circuit output
  Matrix<DCRTPoly> sk(zero_alloc, 2, m);

  // decrypted text
  DCRTPoly dtext(ilDCRTParams, Format::EVALUATION, true);

  ptext.SwitchFormat();

  sender.Encrypt(ilDCRTParams, trapdoorA.first, publicElementB, pubElemBeta,
                 &x[0], ptext, dgg, dug, bug, &ctCin, &c1);

  ctCA = ctCin.ExtractRow(0);

  receiver.NANDGateEvalPKDCRT(ilDCRTParams, publicElementB.ExtractRow(0),
                              publicElementB.ExtractRows(1, 2), &pubElemBf,
                              ilParamsConsolidated);

  receiver.NANDGateEvalCTDCRT(ilDCRTParams, ctCin.ExtractRow(1), &x[1],
                              publicElementB.ExtractRows(1, 2),
                              ctCin.ExtractRows(2, 3), &y, &ctCf,
                              ilParamsConsolidated);

  pkg.KeyGen(ilDCRTParams, trapdoorA.first, pubElemBf, pubElemBeta,
             trapdoorA.second, dgg, &sk);

  receiver.Decrypt(ilDCRTParams, sk, ctCA, ctCf, c1, &dtext);

  Poly dtextPoly(dtext.CRTInterpolate());

  receiver.Decode(&dtextPoly);

  ptext.SwitchFormat();
  if (ptext1.GetValues() == dtextPoly.GetValues()) {
    std::cout << "Success" << std::endl;
  }
}

void KPABE_NANDGATEDCRT_RANDOM(int32_t base, usint k, usint ringDimension) {
  usint n = ringDimension * 2;  // cyclotomic order

  //  usint k = 21;
  usint ell = 4;  // No of attributes

  //  NativeInteger q = NativeInteger::ONE << (k - 1);
  //  q = lbcrypto::FirstPrime<NativeInteger>(k, n);

  NativeInteger q("2101249");

  NativeInteger rootOfUnity(RootOfUnity(n, q));

  //  NativeInteger rootOfUnity("794438271477401");

  double val = q.ConvertToDouble();
  double logTwo = log(val - 1.0) / log(base) + 1.0;
  size_t k_ = (usint)floor(logTwo) + 1;  //  (+1) is For NAF
  std::cout << "q: " << q << std::endl;
  std::cout << "modulus length: " << k_ << std::endl;
  std::cout << "root of unity: " << rootOfUnity << std::endl;
  std::cout << "Standard deviation: " << SIGMA << std::endl;

  //  NativeInteger nextQ = NativeInteger::ONE << (k-1);
  //  nextQ = lbcrypto::NextPrime<NativeInteger>(q, n);
  //  std::cout << "nextQ: " << nextQ << std::endl;

  NativeInteger nextQ("2236417");

  NativeInteger nextRootOfUnity(RootOfUnity<NativeInteger>(n, nextQ));

  //  NativeInteger nextQ2 = NativeInteger::ONE << (k-1);
  //  nextQ2 = lbcrypto::NextPrime<NativeInteger>(nextQ, n);

  NativeInteger nextQ2("2277377");
  NativeInteger nextRootOfUnity2(RootOfUnity<NativeInteger>(n, nextQ2));

  usint m = 3 * k_ + 2;

  std::vector<NativeInteger> moduli;
  std::vector<NativeInteger> roots_Of_Unity;
  moduli.reserve(3);
  roots_Of_Unity.reserve(3);

  moduli.push_back(q);
  moduli.push_back(nextQ);
  moduli.push_back(nextQ2);

  roots_Of_Unity.push_back(rootOfUnity);
  roots_Of_Unity.push_back(nextRootOfUnity);
  roots_Of_Unity.push_back(nextRootOfUnity2);

  BigInteger bigModulus =
      BigInteger("2101249") * BigInteger("2236417") * BigInteger("2277377");

  BigInteger bigRootOfUnity(RootOfUnity(n, bigModulus));

  BinaryUniformGenerator bug = BinaryUniformGenerator();

  auto ilDCRTParams =
      std::make_shared<ILDCRTParams<BigInteger>>(n, moduli, roots_Of_Unity);
  auto ilParamsConsolidated =
      std::make_shared<ILParams>(n, bigModulus, bigRootOfUnity);

  auto zero_alloc = DCRTPoly::Allocator(ilDCRTParams, Format::COEFFICIENT);

  DCRTPoly::DggType dgg = DCRTPoly::DggType(SIGMA);
  DCRTPoly::DugType dug = DCRTPoly::DugType();

  // Trapdoor Generation
  std::pair<Matrix<DCRTPoly>, RLWETrapdoorPair<DCRTPoly>> trapdoorA =
      RLWETrapdoorUtility<DCRTPoly>::TrapdoorGen(
          ilDCRTParams, SIGMA, base);  // A.first is the public element

  DCRTPoly pubElemBeta(dug, ilDCRTParams, Format::EVALUATION);

  Matrix<DCRTPoly> publicElementB(zero_alloc, ell + 1, m);
  Matrix<DCRTPoly> ctCin(zero_alloc, ell + 2, m);
  DCRTPoly c1(dug, ilDCRTParams, Format::EVALUATION);

  KPABE<DCRTPoly, Poly> pkg(SUBGAUSSIAN), sender(SUBGAUSSIAN),
      receiver(SUBGAUSSIAN);

  pkg.Setup(ilDCRTParams, base, ell, dug, &publicElementB);
  sender.Setup(ilDCRTParams, base, ell);
  receiver.Setup(ilDCRTParams, base, ell);

  // Attribute values all are set to 1 for NAND gate evaluation
  std::vector<usint> x(ell + 1);
  x[0] = x[1] = x[2] = 1;
  usint y;

  Poly ptext1(ilParamsConsolidated, Format::COEFFICIENT, true);
  ptext1.SetValues(bug.GenerateVector(ringDimension, bigModulus),
                   Format::COEFFICIENT);

  DCRTPoly ptext(ptext1, ilDCRTParams);

  // circuit outputs
  // evaluated Bs
  Matrix<DCRTPoly> pubElemBf(
      DCRTPoly::Allocator(ilDCRTParams, Format::EVALUATION), 1, m);
  // evaluated Cs
  Matrix<DCRTPoly> ctCf(DCRTPoly::Allocator(ilDCRTParams, Format::EVALUATION),
                        1, m);
  // CA
  Matrix<DCRTPoly> ctCA(DCRTPoly::Allocator(ilDCRTParams, Format::EVALUATION),
                        1, m);

  // secret key corresponding to the circuit output
  Matrix<DCRTPoly> sk(zero_alloc, 2, m);

  // decrypted text
  DCRTPoly dtext(ilDCRTParams, Format::EVALUATION, true);

  ptext.SwitchFormat();

  sender.Encrypt(ilDCRTParams, trapdoorA.first, publicElementB, pubElemBeta,
                 &x[0], ptext, dgg, dug, bug, &ctCin, &c1);

  ctCA = ctCin.ExtractRow(0);

  receiver.NANDGateEvalPKDCRT(ilDCRTParams, publicElementB.ExtractRow(0),
                              publicElementB.ExtractRows(1, 2), &pubElemBf,
                              ilParamsConsolidated, 10);

  receiver.NANDGateEvalCTDCRT(ilDCRTParams, ctCin.ExtractRow(1), &x[1],
                              publicElementB.ExtractRows(1, 2),
                              ctCin.ExtractRows(2, 3), &y, &ctCf,
                              ilParamsConsolidated, 10);

  pkg.KeyGen(ilDCRTParams, trapdoorA.first, pubElemBf, pubElemBeta,
             trapdoorA.second, dgg, &sk);

  receiver.Decrypt(ilDCRTParams, sk, ctCA, ctCf, c1, &dtext);

  Poly dtextPoly(dtext.CRTInterpolate());

  receiver.Decode(&dtextPoly);

  ptext.SwitchFormat();
  if (ptext1.GetValues() == dtextPoly.GetValues()) {
    std::cout << "Success" << std::endl;
  } else {
    std::cout << "Failure" << std::endl;
  }
}

void KPABEANDGate(int32_t base, usint k, usint ringDimension) {
  usint n = ringDimension * 2;
  usint ell = 4;  // No of attributes for AND gate

  BigInteger q = BigInteger(1) << (k - 1);
  q = lbcrypto::FirstPrime<BigInteger>(k, n);
  BigInteger rootOfUnity(RootOfUnity(n, q));

  double val = q.ConvertToDouble();
  double logTwo = log(val - 1.0) / log(base) + 1.0;
  size_t k_ = (usint)floor(logTwo) + 1;
  usint m = k_ + 2;

  auto ilParams = std::make_shared<ILParams>(n, q, rootOfUnity);

  auto zero_alloc = Poly::Allocator(ilParams, Format::COEFFICIENT);

  DiscreteGaussianGenerator dgg = DiscreteGaussianGenerator(SIGMA);
  Poly::DugType dug = Poly::DugType();
  dug.SetModulus(q);
  BinaryUniformGenerator bug = BinaryUniformGenerator();

  // Precompuations for FTT
  ChineseRemainderTransformFTT<BigVector>::PreCompute(rootOfUnity, n, q);

  // Trapdoor Generation
  std::pair<Matrix<Poly>, RLWETrapdoorPair<Poly>> A =
      RLWETrapdoorUtility<Poly>::TrapdoorGen(ilParams, SIGMA, base, true);

  Poly pubElemBeta(dug, ilParams, Format::EVALUATION);

  Matrix<Poly> publicElementB(zero_alloc, ell + 1, m);
  Matrix<Poly> ctCin(zero_alloc, ell + 2, m);
  Poly c1(dug, ilParams, Format::EVALUATION);

  KPABE<Poly, Poly> pkg, sender, receiver;

  pkg.Setup(ilParams, base, ell, dug, &publicElementB);
  sender.Setup(ilParams, base, ell);
  receiver.Setup(ilParams, base, ell);

  // Attribute values all are set to 1 for NAND gate evaluation
  std::vector<usint> x(ell);
  x[0] = x[1] = x[2] = 0;
  usint y;

  // plain text in $R_2$
  Poly ptext(ilParams, Format::COEFFICIENT, true);

  // circuit outputs
  Matrix<Poly> pubElemBf(Poly::Allocator(ilParams, Format::EVALUATION), 1, m);
  Matrix<Poly> ctCf(Poly::Allocator(ilParams, Format::EVALUATION), 1, m);
  Matrix<Poly> ctCA(Poly::Allocator(ilParams, Format::EVALUATION), 1, m);

  // Secret key for the output of the circuit
  Matrix<Poly> sk(zero_alloc, 2, m);

  // text after the decryption
  Poly dtext(ilParams, Format::EVALUATION, true);

  // Encrypt a uniformly randomly selected message ptext (in ptext in $R_2$)
  ptext.SetValues(bug.GenerateVector(ringDimension, q), Format::COEFFICIENT);
  ptext.SwitchFormat();
  sender.Encrypt(ilParams, A.first, publicElementB, pubElemBeta, &x[0], ptext,
                 dgg, dug, bug, &ctCin, &c1);

  ctCA = ctCin.ExtractRow(0);

  receiver.ANDGateEvalPK(ilParams, publicElementB.ExtractRows(1, 2),
                         &pubElemBf);
  receiver.ANDGateEvalCT(ilParams, &x[1], publicElementB.ExtractRows(1, 2),
                         ctCin.ExtractRows(2, 3), &y, &ctCf);

  pkg.KeyGen(ilParams, A.first, pubElemBf, pubElemBeta, A.second, dgg, &sk);

  receiver.Decrypt(ilParams, sk, ctCA, ctCf, c1, &dtext);
  receiver.Decode(&dtext);

  ptext.SwitchFormat();
  if (ptext.GetValues() == dtext.GetValues()) {
    std::cout << "Success" << std::endl;
  }
}

void KPABEANDGateDCRT(int32_t base, usint k, usint ringDimension) {
  usint n = ringDimension * 2;  // cyclotomic order

  //  usint k = 21;
  usint ell = 4;  // No of attributes

  //  NativeInteger q = NativeInteger::ONE << (k - 1);
  //  q = lbcrypto::FirstPrime<NativeInteger>(k, n);

  NativeInteger q("2101249");

  NativeInteger rootOfUnity(RootOfUnity(n, q));

  //  NativeInteger rootOfUnity("794438271477401");

  double val = q.ConvertToDouble();
  double logTwo = log(val - 1.0) / log(base) + 1.0;
  size_t k_ = (usint)floor(logTwo) + 1;  //  (+1) is For NAF
  std::cout << "q: " << q << std::endl;
  std::cout << "modulus length: " << k_ << std::endl;
  std::cout << "root of unity: " << rootOfUnity << std::endl;
  std::cout << "Standard deviation: " << SIGMA << std::endl;

  //  NativeInteger nextQ = NativeInteger::ONE << (k-1);
  //  nextQ = lbcrypto::NextPrime<NativeInteger>(q, n);
  //  std::cout << "nextQ: " << nextQ << std::endl;

  NativeInteger nextQ("2236417");

  NativeInteger nextRootOfUnity(RootOfUnity<NativeInteger>(n, nextQ));

  //  NativeInteger nextQ2 = NativeInteger::ONE << (k-1);
  //  nextQ2 = lbcrypto::NextPrime<NativeInteger>(nextQ, n);

  NativeInteger nextQ2("2277377");
  NativeInteger nextRootOfUnity2(RootOfUnity<NativeInteger>(n, nextQ2));

  usint m = 3 * k_ + 2;

  std::vector<NativeInteger> moduli;
  std::vector<NativeInteger> roots_Of_Unity;
  moduli.reserve(3);
  roots_Of_Unity.reserve(3);

  moduli.push_back(q);
  moduli.push_back(nextQ);
  moduli.push_back(nextQ2);

  roots_Of_Unity.push_back(rootOfUnity);
  roots_Of_Unity.push_back(nextRootOfUnity);
  roots_Of_Unity.push_back(nextRootOfUnity2);

  BigInteger bigModulus =
      BigInteger("2101249") * BigInteger("2236417") * BigInteger("2277377");

  BigInteger bigRootOfUnity(RootOfUnity(n, bigModulus));

  BinaryUniformGenerator bug = BinaryUniformGenerator();

  auto ilDCRTParams =
      std::make_shared<ILDCRTParams<BigInteger>>(n, moduli, roots_Of_Unity);
  auto ilParamsConsolidated =
      std::make_shared<ILParams>(n, bigModulus, bigRootOfUnity);

  auto zero_alloc = DCRTPoly::Allocator(ilDCRTParams, Format::COEFFICIENT);

  DCRTPoly::DggType dgg = DCRTPoly::DggType(SIGMA);
  DCRTPoly::DugType dug = DCRTPoly::DugType();

  // Trapdoor Generation
  // A.first is the public element
  std::pair<Matrix<DCRTPoly>, RLWETrapdoorPair<DCRTPoly>> trapdoorA =
      RLWETrapdoorUtility<DCRTPoly>::TrapdoorGen(ilDCRTParams, SIGMA, base,
                                                 true);

  DCRTPoly pubElemBeta(dug, ilDCRTParams, Format::EVALUATION);

  Matrix<DCRTPoly> publicElementB(zero_alloc, ell + 1, m);
  Matrix<DCRTPoly> ctCin(zero_alloc, ell + 2, m);
  DCRTPoly c1(dug, ilDCRTParams, Format::EVALUATION);

  KPABE<DCRTPoly, Poly> pkg, sender, receiver;

  pkg.Setup(ilDCRTParams, base, ell, dug, &publicElementB);
  sender.Setup(ilDCRTParams, base, ell);
  receiver.Setup(ilDCRTParams, base, ell);

  // Attribute values all are set to 1 for NAND gate evaluation
  std::vector<usint> x(ell);
  x[0] = x[1] = x[2] = 0;
  usint y;

  Poly ptext1(ilParamsConsolidated, Format::COEFFICIENT, true);
  ptext1.SetValues(bug.GenerateVector(ringDimension, bigModulus),
                   Format::COEFFICIENT);

  DCRTPoly ptext(ptext1, ilDCRTParams);

  // circuit outputs
  // evaluated Bs
  Matrix<DCRTPoly> pubElemBf(
      DCRTPoly::Allocator(ilDCRTParams, Format::EVALUATION), 1, m);
  // evaluated Cs
  Matrix<DCRTPoly> ctCf(DCRTPoly::Allocator(ilDCRTParams, Format::EVALUATION),
                        1, m);
  // CA
  Matrix<DCRTPoly> ctCA(DCRTPoly::Allocator(ilDCRTParams, Format::EVALUATION),
                        1, m);

  // secret key corresponding to the circuit output
  Matrix<DCRTPoly> sk(zero_alloc, 2, m);

  // decrypted text
  DCRTPoly dtext(ilDCRTParams, Format::EVALUATION, true);

  ptext.SwitchFormat();

  sender.Encrypt(ilDCRTParams, trapdoorA.first, publicElementB, pubElemBeta,
                 &x[0], ptext, dgg, dug, bug, &ctCin, &c1);

  ctCA = ctCin.ExtractRow(0);

  receiver.ANDGateEvalPKDCRT(ilDCRTParams, publicElementB.ExtractRows(1, 2),
                             &pubElemBf, ilParamsConsolidated);
  receiver.ANDGateEvalCTDCRT(
      ilDCRTParams, &x[1], publicElementB.ExtractRows(1, 2),
      ctCin.ExtractRows(2, 3), &y, &ctCf, ilParamsConsolidated);

  pkg.KeyGen(ilDCRTParams, trapdoorA.first, pubElemBf, pubElemBeta,
             trapdoorA.second, dgg, &sk);

  receiver.Decrypt(ilDCRTParams, sk, ctCA, ctCf, c1, &dtext);

  Poly dtextPoly(dtext.CRTInterpolate());

  receiver.Decode(&dtextPoly);

  ptext.SwitchFormat();
  if (ptext1.GetValues() == dtextPoly.GetValues()) {
    std::cout << "Success" << std::endl;
  }
}
