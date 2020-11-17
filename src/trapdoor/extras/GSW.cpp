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

template <class Integer, class Vector>
void RunHETests(uint32_t n, uint32_t base, const Integer &q, double std,
                GINVERSEMODE mode);

int main() {
  TimeVar t1;

  double timeExec(0);

  uint32_t n = 16;
  uint32_t l = 30;
  uint32_t base = 2;

  double std = 3.19;

  BigInteger q(uint64_t(1) << l);

  std::cout << " ==== BigInteger Deterministic =====" << std::endl;

  TIC(t1);

  RunHETests<BigInteger, BigVector>(n, base, q, std, DETERMINISTIC);

  timeExec = TOC(t1);

  std::cerr << "Time of BigInteger: " << timeExec << "ms" << std::endl;

  std::cout << " ==== BigInteger Randomized =====" << std::endl;

  TIC(t1);

  RunHETests<BigInteger, BigVector>(n, base, q, std, RANDOM);

  timeExec = TOC(t1);

  std::cerr << "Time of BigInteger: " << timeExec << "ms" << std::endl;

  // NativeInteger qNative(uint64_t(1)<<l);
  NativeInteger qNative(1073741827);

  std::cout << " \n==== NativeInteger Deterministic =====" << std::endl;

  TIC(t1);

  RunHETests<NativeInteger, NativeVector>(n, base, qNative, std, DETERMINISTIC);

  timeExec = TOC(t1);

  std::cerr << "Time of NativeInteger: " << timeExec << "ms" << std::endl;

  std::cout << " ==== NativeInteger Randomized =====" << std::endl;

  TIC(t1);

  RunHETests<NativeInteger, NativeVector>(n, base, qNative, std, RANDOM);

  timeExec = TOC(t1);

  std::cerr << "Time of NativeInteger Randomized: " << timeExec << "ms"
            << std::endl;

  return 0;
}

template <class Integer, class Vector>
void RunHETests(uint32_t n, uint32_t base, const Integer &q, double std,
                GINVERSEMODE mode) {
  GSWScheme<Integer, Vector> scheme;

  scheme.Setup(n, base, q, std, mode);

  auto sk = scheme.SecretKeyGen();

  auto c = scheme.Encrypt(Integer(1), sk);

  auto cPlus = scheme.EvalAdd(c, c);

  auto cPlus2 = scheme.EvalAdd(cPlus, c);

  auto cMult = scheme.EvalMult(c, c);

  auto cMultByZero = scheme.EvalMult(c, cPlus);

  auto cMult3 = scheme.EvalMult(cMult, cMult);

  auto cMult3Alt = scheme.EvalMult(cMult, cMultByZero);

  auto p = scheme.Decrypt(c, sk);

  std::cout << "plaintext = " << p << std::endl;

  auto pPlus = scheme.Decrypt(cPlus, sk);

  std::cout << "1 + 1 mod 2 = " << pPlus << std::endl;

  auto pPlus2 = scheme.Decrypt(cPlus2, sk);

  std::cout << "1 + 1 + 1 mod 2 = " << pPlus2 << std::endl;

  auto pMult = scheme.Decrypt(cMult, sk);

  std::cout << "1 * 1 mod 2 = " << pMult << std::endl;

  auto pMultByZero = scheme.Decrypt(cMultByZero, sk);

  std::cout << "1 * 0 mod 2 = " << pMultByZero << std::endl;

  auto pMult3 = scheme.Decrypt(cMult3, sk);

  std::cout << "1 * 1 * 1 * 1 mod 2 = " << pMult3 << std::endl;

  auto pMult3Alt = scheme.Decrypt(cMult3Alt, sk);

  std::cout << "1 * 1 * 1 * 0 mod 2 = " << pMult3Alt << std::endl;
}
