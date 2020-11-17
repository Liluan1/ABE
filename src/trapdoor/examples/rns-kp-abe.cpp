// @file rns-kp-abe.cpp - Example of Key-Policy Attribute-Based Encryption
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

#include <fstream>
#include <iostream>

#include "utils/debug.h"

#include "utils/parallel.h"

#define PROFIILE

#include "abe/kp_abe_rns.h"

using namespace lbcrypto;

int KPABE_BenchmarkCircuitTestDCRT(usint iter, int32_t base);
usint EvalNANDTree(usint *x, usint ell);

int main() {
  PalisadeParallelControls.Enable();
  PseudoRandomNumberGenerator::InitPRNG();
  KPABE_BenchmarkCircuitTestDCRT(5, 1 << 20);

  return 0;
}

int KPABE_BenchmarkCircuitTestDCRT(usint iter, int32_t base) {
  usint n = 1 << 12;  // cyclotomic order
  size_t kRes = 50;   // CRT modulus size
  usint ell = 4;      // No of attributes
  size_t size = 2;    // Number of CRT moduli

  std::cout << "Number of attributes: " << ell << std::endl;

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

  auto zero_alloc = DCRTPoly::Allocator(ilDCRTParams, Format::COEFFICIENT);

  DCRTPoly::DggType dgg = DCRTPoly::DggType(SIGMA);
  DCRTPoly::DugType dug = DCRTPoly::DugType();
  DCRTPoly::BugType bug = DCRTPoly::BugType();

  // Trapdoor Generation
  std::pair<Matrix<DCRTPoly>, RLWETrapdoorPair<DCRTPoly>> trapdoorA =
      RLWETrapdoorUtility<DCRTPoly>::TrapdoorGen(
          ilDCRTParams, SIGMA, base);  // A.first is the public element

  DCRTPoly pubElemBeta(dug, ilDCRTParams, Format::EVALUATION);

  Matrix<DCRTPoly> publicElementB(zero_alloc, ell + 1, m);
  Matrix<DCRTPoly> ctCin(zero_alloc, ell + 2, m);
  DCRTPoly c1(dug, ilDCRTParams, Format::EVALUATION);

  KPABErns pkg, sender, receiver;

  pkg.Setup(ilDCRTParams, base, ell, dug, &publicElementB);
  sender.Setup(ilDCRTParams, base, ell);
  receiver.Setup(ilDCRTParams, base, ell);

  // Attribute values all are set to 1 for NAND gate Format::EVALUATION
  std::vector<usint> x(ell + 1);
  x[0] = 1;

  usint found = 0;
  while (found == 0) {
    for (usint i = 1; i < ell + 1; i++)
      // x[i] = rand() & 0x1;
      x[i] = bug.GenerateInteger().ConvertToInt();
    if (EvalNANDTree(&x[1], ell) == 0) found = 1;
  }

  usint y;

  TimeVar t1;
  double avg_keygen(0.0), avg_evalct(0.0), avg_evalpk(0.0), avg_enc(0.0),
      avg_dec(0.0);

  // plaintext
  for (usint i = 0; i < iter; i++) {
    std::cout << "running iter " << i + 1 << std::endl;

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
    TIC(t1);
    sender.Encrypt(ilDCRTParams, trapdoorA.first, publicElementB, pubElemBeta,
                   &x[0], ptext, dgg, dug, bug, &ctCin,
                   &c1);  // Cin and c1 are the ciphertext
    avg_enc += TOC(t1);

    ctCA = ctCin.ExtractRow(0);  // CA is A^T * s + e 0,A

    TIC(t1);
    receiver.EvalCT(ilDCRTParams, publicElementB, &x[0],
                    ctCin.ExtractRows(1, ell + 1), &y, &evalCf);
    avg_evalct += TOC(t1);

    TIC(t1);
    pkg.EvalPK(ilDCRTParams, publicElementB, &evalBf);
    avg_evalpk += TOC(t1);

    TIC(t1);
    pkg.KeyGen(ilDCRTParams, trapdoorA.first, evalBf, pubElemBeta,
               trapdoorA.second, dgg, &sk);
    avg_keygen += TOC(t1);
    //  CheckSecretKeyKPDCRT(m, trapdoorA.first, evalBf, sk, pubElemBeta);

    TIC(t1);
    receiver.Decrypt(ilDCRTParams, sk, ctCA, evalCf, c1, &dtext);
    avg_dec += TOC_US(t1);

    NativeVector ptext2 = ptext.GetValues();
    ptext2.SetModulus(NativeInteger(2));

    if (ptext2 != dtext.GetValues()) {
      std::cout << "Decryption fails at iteration: " << i << std::endl;
      // std::cerr << ptext << std::endl;
      // std::cerr << dtext << std::endl;
      return 0;
    }

    // std::cerr << ptext << std::endl;
    // std::cerr << dtext << std::endl;
  }

  std::cout << "Encryption is successful after " << iter << " iterations!\n";
  std::cout << "Average key generation time : "
            << "\t" << (avg_keygen) / iter << " ms" << std::endl;
  std::cout << "Average ciphertext Format::EVALUATION time : "
            << "\t" << (avg_evalct) / iter << " ms" << std::endl;
  std::cout << "Average public key Format::EVALUATION time : "
            << "\t" << (avg_evalpk) / iter << " ms" << std::endl;
  std::cout << "Average encryption time : "
            << "\t" << (avg_enc) / iter << " ms" << std::endl;
  std::cout << "Average decryption time : "
            << "\t" << (avg_dec) / (iter * 1000) << " ms" << std::endl;

  return 0;
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
