// @file obfuscatelp.h Abstract classes for conjunction obfuscator
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

#ifndef LBCRYPTO_OBFUSCATE_OBFUSCATELP_H
#define LBCRYPTO_OBFUSCATE_OBFUSCATELP_H

#include <string>
#include <vector>

#include "encoding/plaintext.h"
#include "lattice/elemparams.h"
#include "lattice/ildcrtparams.h"
#include "lattice/ilelement.h"
#include "lattice/ilparams.h"
#include "math/distrgen.h"
#include "utils/inttypes.h"
#include "utils/serial.h"
#include "utils/serializable.h"

/**
 * @namespace lbcrypto
 * The namespace of lbcrypto
 */
namespace lbcrypto {

/**
 * @brief Abstract interface class for patterns (common for cleartext and
 * conjunction patterns)
 * @tparam Element a ring element.
 */
template <class Element>
class Pattern {};

/**
 * @brief Abstract interface class for conjunction patterns (common for
 * cleartext and conjunction patterns)
 * @tparam Element a ring element.
 */
template <class Element>
class ConjunctionPattern : public Pattern<Element> {
 public:
  virtual ~ConjunctionPattern() {}

  /**
   * Method to return the length of the pattern.
   *
   * @return the length of the pattern.
   */
  virtual usint GetLength() const = 0;

  template <class Archive>
  void save(Archive& ar, std::uint32_t const version) const {}

  template <class Archive>
  void load(Archive& ar, std::uint32_t const version) {}

  std::string SerializedObjectName() const { return "ConjunctionPattern"; }
};

/**
 * @brief Class for cleartext patterns; includes methods specific to clear
 * patterns
 * @tparam Element a ring element.
 */
template <class Element>
class ClearPattern {};

/**
 * @brief Class for obfuscated patterns; includes methods specific to obfuscated
 * patterns
 * @tparam Element a ring element.
 */
template <class Element>
class ObfuscatedPattern {
 public:
  virtual ~ObfuscatedPattern() {}

  template <class Archive>
  void save(Archive& ar, std::uint32_t const version) const {}

  template <class Archive>
  void load(Archive& ar, std::uint32_t const version) {}

  std::string SerializedObjectName() const { return "ObfuscatedPattern"; }
};

}  // namespace lbcrypto
#endif

CEREAL_REGISTER_TYPE(lbcrypto::ConjunctionPattern<lbcrypto::Poly>);
CEREAL_REGISTER_TYPE(lbcrypto::ConjunctionPattern<lbcrypto::NativePoly>);
CEREAL_REGISTER_TYPE(lbcrypto::ConjunctionPattern<lbcrypto::DCRTPoly>);

CEREAL_REGISTER_TYPE(lbcrypto::ObfuscatedPattern<lbcrypto::Poly>);
CEREAL_REGISTER_TYPE(lbcrypto::ObfuscatedPattern<lbcrypto::NativePoly>);
CEREAL_REGISTER_TYPE(lbcrypto::ObfuscatedPattern<lbcrypto::DCRTPoly>);
