// @file cpa-abe.cpp - Example for ciphertext-policy attribute based encryption
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

#include "abecontext.h"
#include "lattice/backend.h"

using namespace lbcrypto;

ABEContext<NativePoly> context;
CPABEMasterPublicKey<NativePoly> mpk;
CPABEMasterSecretKey<NativePoly> msk;
usint ringsize = 1024;
usint base = 64;
TimeVar t1;

void iterTest(double durations[][4], usint numAttributes,
              std::vector<int64_t> &vectorOfInts) {
  // Create context under security level and number of attributes
  // std::cout << "This is a demo file of the CPABE scheme" << std::endl
  //           << std::endl;
  std::cout << "=======================================" << std::endl;
  // usint numAttributes = 16;
  // std::cout << "Used parameters:" << std::endl;
  // std::cout << "Ring size: " << ringsize << std::endl;
  std::cout << "Number of attributes: " << numAttributes << std::endl;
  // std::cout << "Base: " << base << std::endl << std::endl;

  // std::cout << "Generating a context under these parameters" << std::endl
  //           << std::endl;

  // std::cout << "Generating master secret key and master public key"
  //           << std::endl;
  // Generate master keys
  // TIC(t1);

  // double duration = TOC_US(t1);
  // durations[numAttributes / 2 - 1][0] += duration;
  // std::cout << "Setup: " << duration << " us" << std::endl;

  // Create a random access policy and user attribute set
  // std::cout << " Creating access policy and user attribute sets" <<
  // std::endl;
  std::vector<usint> s(numAttributes);
  std::vector<int> w(numAttributes);

  for (usint j = 0; j < numAttributes; j++) s[j] = rand() % 2;
  // for (usint j = 0; j < numAttributes; j++) w[j] = rand() % 3 - 1;

  for (usint j = 0; j < numAttributes; j++) {
    if (rand() % 2 == 0) {
      w[j] = s[j];
    } else {
      w[j] = s[j] - 1;
    }
  }

  // for (usint j = 0; j < 6; j++) w[j] = s[j];

  // for (usint j = 0; j < 6; j++)
  //   if (w[j] == 1) {
  //     w[j] = 0;
  //     break;
  //   }
  // for (usint j = 0; j < 6; j++)
  //   if (s[j] == 0) {
  //     w[j] = -1;
  //     break;
  //   }
  std::cout << "User attribute set: " << s << std::endl;
  std::cout << "Access policy defined:" << w << std::endl;
  CPABEUserAccess<NativePoly> ua(s);
  CPABEAccessPolicy<NativePoly> ap(w);

  // Create the key corresponding to the access policy
  CPABESecretKey<NativePoly> sk;
  // std::cout << "Creating secret key for the attribute set" << std::endl;
  TIC(t1);
  context.KeyGen(msk, mpk, ua, &sk);
  double duration = TOC_US(t1);
  durations[numAttributes / 2 - 1][0] += duration;
  std::cout << "KeyGen: " << duration << " us" << std::endl;

  // Create a plaintext
  // std::vector<int64_t> vectorOfInts = {1, 0, 0, 1, 1, 0, 1, 0, 1, 0};
  Plaintext pt = context.MakeCoefPackedPlaintext(vectorOfInts);
  std::cout << "Plaintext vector of bits: " << vectorOfInts.size() << std::endl;
  //           << std::endl;

  // Encrypt the plaintext
  // std::cout << "Encrypting the plaintext under the access policy" <<
  // std::endl;
  TIC(t1);
  CPABECiphertext<NativePoly> ct;
  NativePoly secret;
  context.Encrypt(mpk, ap, pt, &ct, &secret);
  duration = TOC_US(t1);
  durations[numAttributes / 2 - 1][1] += duration;
  std::cout << "Encryption: " << duration << " us" << std::endl;

  // Decrypt the ciphertext
  // std::cout << "Decrpyting the ciphertext" << std::endl;
  TIC(t1);
  Plaintext dt = context.Decrypt(ap, ua, sk, ct);
  duration = TOC_US(t1);
  durations[numAttributes / 2 - 1][2] += duration;
  std::cout << "Decryption: " << duration << " us" << std::endl;

  // std::cout << "Checking if the plaintext & decrypted text match" <<
  // std::endl; Check if original plaintext and decrypted plaintext match
  if (pt->GetElement<NativePoly>() == dt->GetElement<NativePoly>()) {
    std::cout << "Encryption & decryption successful" << std::endl;
  } else {
    std::cout << "Encryption & decryption failed" << std::endl;
  }
  std::cout << "Updating the ciphertext" << std::endl;
  std::vector<usint> new_s(numAttributes);
  std::vector<int> new_w(numAttributes);

  for (usint j = 0; j < numAttributes; j++) new_s[j] = rand() % 2;

  for (usint j = 0; j < numAttributes; j++) new_w[j] = new_s[j];

  for (usint j = 0; j < numAttributes; j++)
    if (new_w[j] == 1) {
      new_w[j] = 0;
      break;
    }
  for (usint j = 0; j < numAttributes; j++)
    if (new_s[j] == 0) {
      new_w[j] = -1;
      break;
    }
  std::cout << "User attribute set: " << new_s << std::endl;
  std::cout << "Access policy defined:" << new_w << std::endl;
  CPABEUserAccess<NativePoly> new_ua(new_s);
  CPABEAccessPolicy<NativePoly> new_ap(new_w);
  TIC(t1);
  context.Update(mpk, new_ap, &ct, &secret);
  duration = TOC_US(t1);
  durations[numAttributes / 2 - 1][3] += duration;
  std::cout << "Update: " << duration << " us" << std::endl << std::endl;
}

int main() {
  double durations[8][4] = {{0}};
  const int vectorSize = 16 * 8;
  int arr[vectorSize] = {0};
  // std::cout << vectorOfInts << std::endl;
  for (int i = 1; i <= 8; ++i) {
    context.GenerateCPABEContext(i * 2, ringsize, base);
    context.Setup(&mpk, &msk);
    for (int j = 1; j <= 100; ++j) {
      for (int k = 0; k < vectorSize; ++k) {
        arr[k] = rand() % 2;
      }
      std::vector<int64_t> vectorOfInts(arr, arr + vectorSize);
      iterTest(durations, i * 2, vectorOfInts);
    }
  }
  std::cout << "=======================================" << std::endl;
  std::cout.width(8);
  std::cout << "keygen "
            << "encrypt "
            << "decrypt "
            << "update" << std::endl;
  for (int i = 0; i < 8; ++i) {
    for (int j = 0; j < 4; ++j) {
      std::cout << durations[i][j] / 100 << ' ';
    }
    std::cout << std::endl;
  }
  std::cout << "=======================================" << std::endl;
  return 0;
}