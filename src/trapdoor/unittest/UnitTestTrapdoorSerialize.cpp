// @file
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
/*
 This code exercises the (de)serialization of objects used in the PALISADE
 trapdoor components.
*/

#include <iostream>
#include "gtest/gtest.h"

#include "lattice/backend.h"
#include "lattice/elemparamfactory.h"
#include "lattice/trapdoor.h"
#include "math/backend.h"
#include "math/distrgen.h"
#include "math/nbtheory.h"
#include "utils/inttypes.h"
#include "utils/parmfactory.h"
#include "utils/serialize-binary.h"
#include "utils/utilities.h"

using namespace std;
using namespace lbcrypto;

class UnitTestSerialize : public ::testing::Test {
 protected:
  virtual void SetUp() {}

  virtual void TearDown() {
    // Code here will be called immediately after each test
    // (right before the destructor).
  }
};

TEST(UTTDSer, serialize_vector_RLWETrapdoorPair) {
  // Serialize/DeserializeVectorOfRLWETrapdoorPair is a helper function to test
  // note the object has to be created outside of the function.

  DEBUG_FLAG(false);
  const int vecsize = 4;

  DEBUG("step 0");
  // build test vector (note needs allocator for Matrix<>

  vector<RLWETrapdoorPair<BigInteger>> testvec(vecsize,
                                               RLWETrapdoorPair<BigInteger>());
  vector<RLWETrapdoorPair<BigInteger>> newvec(vecsize,
                                              RLWETrapdoorPair<BigInteger>());

  // zero matricies
  Matrix<BigInteger> zeromat(BigInteger::Allocator, 0, 0);
  // Matrix<BigInteger> *tm_p; //pointer to a M<I>

  DEBUG("step 1");
  // build test input matricies
  usint nrows(3);
  usint ncols(5);

  DEBUG("step 2");

  DEBUG("step 3");

  for (usint i = 0; i < vecsize; i++) {
    // point to zero matricies
    auto tm_p = make_shared<Matrix<BigInteger>>(BigInteger::Allocator, 0, 0);
    tm_p->SetSize(nrows + i, ncols + i);
    for (usint row = 0; row < nrows + i; row++) {
      for (usint col = 0; col < ncols + i; col++) {
        // write a unique value
        (*tm_p)(row, col) = BigInteger(100 * i + 10 * row + col);
      }
    }
    testvec[i].m_r = *tm_p;

    for (usint row = 0; row < nrows + i; row++) {
      for (usint col = 0; col < ncols + i; col++) {
        // write a unique value
        (*tm_p)(row, col) = BigInteger(1000 * i + 10 * row + col);
      }
    }
    testvec[i].m_e = *tm_p;

    newvec[i].m_r = zeromat;
    newvec[i].m_e = zeromat;
  }

  DEBUG("step 4");
  stringstream ss;

  // serialize the vector
  Serial::Serialize(testvec, ss, SerType::BINARY);

#if !defined(NDEBUG)
  if (dbg_flag) {
    // write the result to cout for debug
    std::cout << Serial::SerializeToString(testvec) << std::endl;
  }
#endif

  DEBUG("step 5");
  Serial::Deserialize(newvec, ss, SerType::BINARY);

  DEBUG("step 6");

  // loop over vector  and dereference matricies and compare
  auto it_1 = testvec.begin();
  auto it_2 = newvec.begin();
  auto i = 0;
  for (; (it_1 != testvec.end()) && (it_2 != newvec.end());
       it_1++, it_2++, i++) {
    DEBUG("testing iteration " << i);
    // compare dereferenced matricies
    EXPECT_EQ(it_1->m_r, it_2->m_r)
        << "Mismatch after ser/deser in entry " << i << " m_r ";

    EXPECT_EQ(it_1->m_e, it_2->m_e)
        << "Mismatch after ser/deser in entry " << i << " m_e ";
  }
}
