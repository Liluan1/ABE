// @file UnitTestKPABErns.cpp - Unit tests for the KPABE scheme
//
// @copyright Copyright (c) 2019, Duality Technologies Inc.
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

#include "gtest/gtest.h"
#include <iostream>
#include <vector>
#include <list>
#include <ctime>

#include "utils/debug.h"

#include "utils/parallel.h"

#include "palisade.h"
#include "abe/kp_abe_rns.h"

using namespace std;
using namespace lbcrypto;

class UTKPABErns : public ::testing::Test {
 public:
  UTKPABErns() {}
  ~UTKPABErns() {}

 protected:
  void SetUp() {}

  void TearDown() {}
};

static usint EvalNANDTree(usint* x, usint ell) {
  usint y;

  if (ell == 2) {
    y = 1 - x[0] * x[1];
  } else {
    ell >>= 1;
    y = 1 - (EvalNANDTree(&x[0], ell) * EvalNANDTree(&x[ell], ell));
  }
  return y;
}

static void checkEquality(const NativeVector& a, const NativeVector& b,
                          const string& failmsg) {
  int vectorSize = a.GetLength();
  std::vector<usint> allTrue(vectorSize);
  std::vector<usint> tmp(vectorSize);
  for (int i = 0; i < vectorSize; i++) {
    allTrue[i] = 1;
    tmp[i] = (a[i] == b[i]);
  }
  EXPECT_TRUE(tmp == allTrue) << failmsg;
}

static void UnitTest_Eval_Nand_Tree(usint n, usint bitsize, usint attributeNum,
                                    usint primeNum, int32_t base) {
  std::vector<NativeInteger> moduli;
  std::vector<NativeInteger> roots_Of_Unity;

  NativeInteger firstInteger = FirstPrime<NativeInteger>(bitsize, 2 * n);

  NativeInteger q = PreviousPrime<NativeInteger>(firstInteger, 2 * n);
  moduli.push_back(q);
  roots_Of_Unity.push_back(RootOfUnity<NativeInteger>(2 * n, moduli[0]));

  NativeInteger prevQ = q;
  for (size_t i = 1; i < primeNum; i++) {
    prevQ = lbcrypto::PreviousPrime<NativeInteger>(prevQ, 2 * n);
    NativeInteger nextRootOfUnity(RootOfUnity<NativeInteger>(2 * n, prevQ));
    moduli.push_back(prevQ);
    roots_Of_Unity.push_back(nextRootOfUnity);
  }

  auto ilDCRTParams =
      std::make_shared<ILDCRTParams<BigInteger>>(2 * n, moduli, roots_Of_Unity);

  ChineseRemainderTransformFTT<NativeVector>::PreCompute(roots_Of_Unity, 2 * n,
                                                         moduli);

  size_t digitCount = static_cast<long>(
      ceil(log2(ilDCRTParams->GetParams()[0]->GetModulus().ConvertToDouble()) /
           log2(base)));
  size_t k = digitCount * ilDCRTParams->GetParams().size();

  size_t m = k + 2;

  auto zero_alloc = DCRTPoly::Allocator(ilDCRTParams, Format::COEFFICIENT);

  DCRTPoly::DggType dgg(SIGMA);
  DCRTPoly::DugType dug;
  DCRTPoly::BugType bug;

  // Trapdoor Generation
  std::pair<Matrix<DCRTPoly>, RLWETrapdoorPair<DCRTPoly>> trapdoorA =
      RLWETrapdoorUtility<DCRTPoly>::TrapdoorGen(
          ilDCRTParams, SIGMA, base);  // A.first is the public element

  DCRTPoly pubElemBeta(dug, ilDCRTParams, Format::EVALUATION);

  Matrix<DCRTPoly> publicElementB(zero_alloc, attributeNum + 1, m);
  Matrix<DCRTPoly> ctCin(zero_alloc, attributeNum + 2, m);
  DCRTPoly c1(dug, ilDCRTParams, Format::EVALUATION);

  KPABErns pkg, sender, receiver;

  pkg.Setup(ilDCRTParams, base, attributeNum, dug, &publicElementB);
  sender.Setup(ilDCRTParams, base, attributeNum);
  receiver.Setup(ilDCRTParams, base, attributeNum);

  std::vector<usint> x(attributeNum + 1);
  x[0] = 1;

  do {
    for (usint i = 1; i < attributeNum + 1; i++)
      x[i] = bug.GenerateInteger().ConvertToInt();
  } while (EvalNANDTree(&x[1], attributeNum) != 0);

  usint y;

  for (usint i = 0; i < 5; i++) {
    NativePoly ptext(bug, ilDCRTParams->GetParams()[0], Format::COEFFICIENT);

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
    NativePoly dtext;

    // Switches to Format::EVALUATION representation
    // ptext.SwitchFormat();
    sender.Encrypt(ilDCRTParams, trapdoorA.first, publicElementB, pubElemBeta,
                   &x[0], ptext, dgg, dug, bug, &ctCin,
                   &c1);  // Cin and c1 are the ciphertext

    ctCA = ctCin.ExtractRow(0);  // CA is A^T * s + e 0,A

    receiver.EvalCT(ilDCRTParams, publicElementB, &x[0],
                    ctCin.ExtractRows(1, attributeNum + 1), &y, &evalCf);

    pkg.EvalPK(ilDCRTParams, publicElementB, &evalBf);

    pkg.KeyGen(ilDCRTParams, trapdoorA.first, evalBf, pubElemBeta,
               trapdoorA.second, dgg, &sk);

    receiver.Decrypt(ilDCRTParams, sk, ctCA, evalCf, c1, &dtext);

    NativeVector ptext2 = ptext.GetValues();

    ptext2.SetModulus(NativeInteger(2));
    checkEquality(ptext2, dtext.GetValues(), "EvalNand fails");
  }
}

TEST(UTKPABErns, kpabe_nandtree_attr_2_primes_2) {
  UnitTest_Eval_Nand_Tree(1 << 12, 50, 2, 2, 1 << 20);
}
TEST(UTKPABErns, kpabe_nandtree_attr_4_primes_2) {
  UnitTest_Eval_Nand_Tree(1 << 12, 50, 4, 2, 1 << 20);
}
TEST(UTKPABErns, kpabe_nandtree_attr_8_primes_3) {
  UnitTest_Eval_Nand_Tree(1 << 12, 50, 8, 3, 1 << 20);
}
