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
#include "utils/serial.h"

using namespace lbcrypto;

const std::string DATAFOLDER = "serial-data";

int main() {
  // Create context under security level and number of attributes
  std::cout << "This is a demo file of the CPABE scheme" << std::endl
            << std::endl;
  usint ringsize = 1024;
  usint numAttributes = 6;
  usint base = 64;
  TimeVar t1;
  std::cout << "Used parameters:" << std::endl;
  std::cout << "Ring size: " << ringsize << std::endl;
  std::cout << "Number of attributes: " << numAttributes << std::endl;
  std::cout << "Base: " << base << std::endl << std::endl;

  ABEContext<NativePoly> context;
  std::cout << "Generating a context under these parameters" << std::endl
            << std::endl;
  context.GenerateCPABEContext(numAttributes, ringsize, base);

  std::cout << "Generating master secret key and master public key"
            << std::endl;

  // Generate master keys
  TIC(t1);
  CPABEMasterPublicKey<NativePoly> mpk;
  CPABEMasterSecretKey<NativePoly> msk;
  context.Setup(&mpk, &msk);
  double duration = TOC(t1);
  std::cout << "Setup: " << duration << " ms" << std::endl << std::endl;

  // Serialize master public keys
  // if (!Serial::SerializeToFile(DATAFOLDER + "/mpk.bin", mpk,
  // SerType::BINARY)) {
  //   std::cerr << "Error writing serialization of the master public key to "
  //                "master-public-key.bin"
  //             << std::endl;
  //   return 1;
  // }

  // Serialize master secret keys
  // if (!Serial::SerializeToFile(DATAFOLDER + "/master-secret-key.bin", msk,
  //                              SerType::BINARY)) {
  //   std::cerr << "Error writing serialization of the master secret key to "
  //                "master-secret-key.bin"
  //             << std::endl;
  //   return 1;
  // }

  // Deserialize master public keys
  // CPABEMasterPublicKey<NativePoly> s_mpk;
  // if (!Serial::DeserializeFromFile(DATAFOLDER + "/master-public-key.bin",
  // s_mpk,
  //                                  SerType::BINARY)) {
  //   std::cerr << "Error reading serialization of the master public key from "
  //                "master-public-key.bin"
  //             << std::endl;
  //   return 1;
  // }

  // Deserialize master secret keys
  // CPABEMasterSecretKey<NativePoly> s_msk;
  // if (!Serial::DeserializeFromFile(DATAFOLDER + "/master-secret-key.bin",
  // s_msk,
  //                                  SerType::BINARY)) {
  //   std::cerr << "Error reading serialization of the master secret key from "
  //                "master-secret-key.bin"
  //             << std::endl;
  //   return 1;
  // }

  // Test master public key serialization
  // if (mpk.GetBNeg() == s_mpk.GetBNeg() && mpk.GetBPos() == s_mpk.GetBPos() &&
  //     mpk.GetPubElemD() == s_mpk.GetPubElemD() && mpk.GetA() == s_mpk.GetA())
  //     {
  //   std::cout << "Master public key serialization test passed" << std::endl
  //             << std::endl;
  // } else {
  //   std::cout << "Master public key serialization test failed" << std::endl
  //             << std::endl;
  //   return 1;
  // }

  // Test master secret key serialization
  // if (msk.GetTA().m_r == s_msk.GetTA().m_r &&
  //     msk.GetTA().m_e == s_msk.GetTA().m_e) {
  //   std::cout << "Master secret key serialization test passed" << std::endl
  //             << std::endl;
  // } else {
  //   std::cout << "Master secret key serialization test failed" << std::endl
  //             << std::endl;
  //   return 1;
  // }

  // Create a random access policy and user attribute set
  std::cout << "Creating access policy and user attribute sets" << std::endl;
  std::vector<usint> s(6);
  std::vector<int> w(6);

  for (usint j = 0; j < 6; j++) s[j] = rand() % 2;

  for (usint j = 0; j < 6; j++) w[j] = s[j];

  for (usint j = 0; j < 6; j++)
    if (w[j] == 1) {
      w[j] = 0;
      break;
    }
  for (usint j = 0; j < 6; j++)
    if (s[j] == 0) {
      w[j] = -1;
      break;
    }
  std::cout << "User attribute set: " << s << std::endl;
  std::cout << "Access policy defined:" << w << std::endl << std::endl;
  CPABEUserAccess<NativePoly> ua(s);
  CPABEAccessPolicy<NativePoly> ap(w);

  // Create the key corresponding to the access policy
  CPABESecretKey<NativePoly> sk;
  std::cout << "Creating secret key for the attribute set" << std::endl;
  TIC(t1);
  context.KeyGen(msk, mpk, ua, &sk);
  duration = TOC(t1);
  std::cout << "KeyGen: " << duration << " ms" << std::endl << std::endl;

  // Serialize secret key
  // if (!Serial::SerializeToFile(DATAFOLDER + "/secret-key.bin", sk,
  //                              SerType::BINARY)) {
  //   std::cerr
  //       << "Error writing serialization of the secret key to secret-key.bin"
  //       << std::endl;
  //   return 1;
  // }

  // Deserialize secret key
  // CPABESecretKey<NativePoly> s_sk;
  // if (!Serial::DeserializeFromFile(DATAFOLDER + "/secret-key.bin", s_sk,
  //                                  SerType::BINARY)) {
  //   std::cerr
  //       << "Error reading serialization of the secret key from
  //       secret-key.bin"
  //       << std::endl;
  //   return 1;
  // }

  // Test secret key serialization
  // if (s_sk.GetSK() == sk.GetSK()) {
  //   std::cout << "Secret key serialization test passed" << std::endl
  //             << std::endl;
  // } else {
  //   std::cout << "Secret key serialization test failed" << std::endl
  //             << std::endl;
  //   return 1;
  // }

  // Create a plaintext
  std::string str = "Hello World";
  std::vector<int64_t> vectorOfInts;
  for (usint i = 0; i < str.length(); ++i) {
    std::bitset<8> bs4(str[i]);
    std::string binStr = bs4.to_string();
    for (usint i = 0; i < binStr.length(); ++i) {
      vectorOfInts.push_back(binStr[i] - '0');
    }
  }
  // std::vector<int64_t> vectorOfInts = {1, 0, 0, 1, 1, 0, 1, 0, 1, 0};
  Plaintext pt = context.MakeCoefPackedPlaintext(vectorOfInts);
  std::cout << "Plaintext vector of bits: " << vectorOfInts << std::endl
            << std::endl;

  // Encrypt the plaintext
  std::cout << "Encrypting the plaintext under the access policy" << std::endl;
  TIC(t1);
  CPABECiphertext<NativePoly> ct;
  NativePoly secret;
  context.Encrypt(mpk, ap, pt, &ct, &secret);
  std::cout << "secret: " << secret << std::endl;
  duration = TOC(t1);
  std::cout << "Encryption: " << duration << " ms" << std::endl << std::endl;

  // Serialize ciphertext
  // if (!Serial::SerializeToFile(DATAFOLDER + "/ciphertext.bin", ct,
  //                              SerType::BINARY)) {
  //   std::cerr
  //       << "Error writing serialization of the ciphertext to ciphertext.bin"
  //       << std::endl;
  //   return 1;
  // }

  // Deserialize ciphertext
  // CPABECiphertext<NativePoly> s_ct;
  // if (!Serial::DeserializeFromFile(DATAFOLDER + "/ciphertext.bin", s_ct,
  //                                  SerType::BINARY)) {
  //   std::cerr
  //       << "Error reading serialization of the ciphertext from
  //       ciphertext.bin"
  //       << std::endl;
  //   return 1;
  // }

  // Test ciphertext serialization
  // if (s_ct.GetCNeg() == ct.GetCNeg() && s_ct.GetCPos() == ct.GetCPos() &&
  //     s_ct.GetCW() == ct.GetCW() && s_ct.GetC1() == ct.GetC1()) {
  //   std::cout << "Ciphertext serialization test passed" << std::endl
  //             << std::endl;
  // } else {
  //   std::cout << "Ciphertext serialization test failed" << std::endl
  //             << std::endl;
  //   return 1;
  // }

  // Decrypt the ciphertext
  std::cout << "Decrpyting the ciphertext" << std::endl;
  TIC(t1);
  Plaintext dt = context.Decrypt(ap, ua, sk, ct);
  duration = TOC(t1);
  std::cout << "Decryption: " << duration << " ms" << std::endl << std::endl;

  std::cout << "Checking if the plaintext & decrypted text match" << std::endl;
  // Check if original plaintext and decrypted plaintext match
  if (pt->GetElement<NativePoly>() == dt->GetElement<NativePoly>()) {
    std::cout << "Encryption & decryption successful" << std::endl;
    dt->Decode();
    std::vector<int64_t> dtv = dt->GetCoefPackedValue();
    size_t i = dt->GetLength();
    while (--i > 0)
      if (dtv[i] != 0) break;
    dtv.resize((i + 1) % 8 == 0 ? i + 1 : ((i + 1) / 8 + 1) * 8);
    std::string dts;
    for (size_t i = 0; i < dtv.size(); i += 8) {
      std::bitset<8> bs;
      for (size_t j = 0; j < 8; j++) {
        bs[j] = dtv[i - j + 7];
      }
      // std::cout << "Decrypted char: " << (char)bs.to_ulong() << " ";
      dts += (char)bs.to_ulong();
    }
    std::cout << std::endl;
    std::cout << "Decrypted text: " << dts << std::endl;
  } else {
    std::cout << "Encryption & decryption failed" << std::endl;
  }

  std::cout << "Creating access policy and user attribute sets" << std::endl;
  std::vector<usint> new_s(6);
  std::vector<int> new_w(6);

  for (usint j = 0; j < 6; j++) new_s[j] = rand() % 2;

  for (usint j = 0; j < 6; j++) new_w[j] = new_s[j];

  for (usint j = 0; j < 6; j++)
    if (new_w[j] == 1) {
      new_w[j] = 0;
      break;
    }
  for (usint j = 0; j < 6; j++)
    if (new_s[j] == 0) {
      new_w[j] = -1;
      break;
    }
  std::cout << "User attribute set: " << new_s << std::endl;
  std::cout << "Access policy defined:" << new_w << std::endl << std::endl;
  CPABEUserAccess<NativePoly> new_ua(new_s);
  CPABEAccessPolicy<NativePoly> new_ap(new_w);
  std::cout << "Updating the ciphertext" << std::endl;
  TIC(t1);
  context.Update(mpk, new_ap, &ct, &secret);
  duration = TOC(t1);
  std::cout << "Update: " << duration << " ms" << std::endl << std::endl;
  CPABESecretKey<NativePoly> new_sk;
  context.KeyGen(msk, mpk, new_ua, &new_sk);
  Plaintext new_dt = context.Decrypt(new_ap, new_ua, new_sk, ct);

  if (pt->GetElement<NativePoly>() == new_dt->GetElement<NativePoly>()) {
    std::cout << "Encryption & decryption successful" << std::endl;
    new_dt->Decode();
    std::vector<int64_t> new_dtv = new_dt->GetCoefPackedValue();
    size_t i = new_dt->GetLength();
    while (--i > 0)
      if (new_dtv[i] != 0) break;
    new_dtv.resize((i + 1) % 8 == 0 ? i + 1 : ((i + 1) / 8 + 1) * 8);
    std::string new_dts;
    for (size_t i = 0; i < new_dtv.size(); i += 8) {
      std::bitset<8> bs;
      for (size_t j = 0; j < 8; j++) {
        bs[j] = new_dtv[i - j + 7];
      }
      // std::cout << "Decrypted char: " << (char)bs.to_ulong() << " ";
      new_dts += (char)bs.to_ulong();
    }
    std::cout << std::endl;
    std::cout << "Decrypted text: " << new_dts << std::endl;
  } else {
    std::cout << "Encryption & decryption failed" << std::endl;
  }
  return 0;
}
