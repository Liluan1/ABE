// @file TODO
// @author TPOC: contact@palisade-crypto.org
//
// @copyright Copyright (c) TODO

// define this to enable PROFILELOG and TIC/TOC
#define PROFILE

// native libs
// #include <sys/resource.h>
#include <math.h>
#include <iostream>
#include <ctime>
#include <cstdlib>
#include <fstream>

#include "utils/debug.h"

#include "subgaussian.h"

using namespace std;
using namespace lbcrypto;

void NativeIntegerTest() {
  TimeVar t1;  // for TIC TOC

  double timeEval;

  long b = 1 << 5;
  long q = 1073741827;
  long k = (long)ceil(log2(q) / log2(b));

  long u = pow(3, 5);

  const size_t count = 10000;

  LatticeSubgaussianUtility<NativeInteger> sampler(b, q, k);

  vector<int64_t> nativeOutput(k);

  TIC(t1);  // start timer for total time
  for (size_t i = 0; i < count; i++)
    sampler.InverseG(u, PseudoRandomNumberGenerator::GetPRNG(), &nativeOutput);
  timeEval = TOC_US(t1);

  cout << "PALISADE impl sampling time: " << timeEval / count << " microseconds"
       << std::endl;
  cout << "input = " << u << endl;
  // test the output
  int64_t test1 = 0;
  int64_t b_i1 = 1;
  for (int i = 0; i < k; i++) {
    test1 += nativeOutput[i] * b_i1;
    b_i1 = b_i1 * b;
    //    std::cout << nativeOutput[i] << std::endl;
  }

  cout << "g^t * output 2 = " << test1 << endl;
}

void BigIntegerTest() {
  TimeVar t1;  // for TIC TOC

  double timeEval;

  long b = 1 << 5;

  const size_t count = 10000;
  BigInteger bigModulus = BigInteger("3079705401285115676503558331198");
  long kBig = (long)ceil(log2(bigModulus.ConvertToDouble()) / log2(b));

  LatticeSubgaussianUtility<BigInteger> samplerBig(b, bigModulus, kBig);

  vector<int64_t> nativeOutputBig(kBig);

  TIC(t1);  // start timer for total time
  for (size_t i = 0; i < count; i++)
    samplerBig.InverseG(bigModulus >> 2, PseudoRandomNumberGenerator::GetPRNG(),
                        &nativeOutputBig);
  timeEval = TOC_US(t1);

  cout << "PALISADE impl sampling time: " << timeEval / count << " microseconds"
       << std::endl;
  cout << "input = " << (bigModulus >> 2) << endl;
  // test the output
  BigInteger testBig1 = 0;
  BigInteger bBig_i1 = 1;
  for (int i = 0; i < kBig; i++) {
    if (nativeOutputBig[i] >= 0) {
      testBig1 += BigInteger(nativeOutputBig[i]) * bBig_i1;
      testBig1.ModEq(bigModulus);
    } else {
      testBig1 +=
          BigInteger(bigModulus - BigInteger(-nativeOutputBig[i])) * bBig_i1;
      testBig1.ModEq(bigModulus);
    }
    bBig_i1 = bBig_i1 * BigInteger(b);
    //    std::cout << nativeOutputBig[i] << std::endl;
  }

  cout << "g^t * output 2 = " << testBig1 << endl;
}

int main() {
  NativeIntegerTest();
  BigIntegerTest();
  return 0;
}

// 3/23/2018: We must change the code in BcB to update the target each step!
// 4/1/2018: The code works! Now we must change the modular arithmetic to use
// the residues {-q/2, ... , 0, 1, ... , q/2}. (I'm not sure if this matters.)
