/*
 * Copyright (c) Facebook, Inc. and its affiliates.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#pragma once

#include <stdexcept>
#include <folly/Range.h>
#include <folly/io/IOBuf.h>
#include <folly/portability/OpenSSL.h>
#include <folly/ssl/OpenSSLPtrTypes.h>

namespace folly {
namespace ssl {

/// Warning:
/// These functions are not thread-safe unless you initialize OpenSSL.
class OpenSSLHash {
 public:
  class Digest {
   public:
    Digest() noexcept = default;

    Digest(const Digest& that) { copy_impl(that); }

    Digest(Digest&& that) noexcept { move_impl(std::move(that)); }

    Digest& operator=(const Digest& that) {
      if (this != &that) {
        copy_impl(that);
      }
      return *this;
    }

    Digest& operator=(Digest&& that) noexcept {
      if (this != &that) {
        move_impl(std::move(that));
        that.hash_reset();
      }
      return *this;
    }

    void hash_init(const EVP_MD* md) {
      if (nullptr == ctx_) {
        ctx_.reset(EVP_MD_CTX_new());
        if (nullptr == ctx_) {
          throw_exception<std::runtime_error>(
              "EVP_MD_CTX_new() returned nullptr");
        }
      }
      check_libssl_result(1, EVP_DigestInit_ex(ctx_.get(), md, nullptr));
      md_ = md;
    }

    void hash_update(ByteRange data) {
      if (nullptr == ctx_) {
        throw_exception<std::runtime_error>(
            "hash_update() called without hash_init()");
      }
      check_libssl_result(
          1, EVP_DigestUpdate(ctx_.get(), data.data(), data.size()));
    }

    void hash_update(const IOBuf& data) {
      for (auto r : data) {
        hash_update(r);
      }
    }

    void hash_final(MutableByteRange out) {
      if (nullptr == ctx_) {
        throw_exception<std::runtime_error>(
            "hash_final() called without hash_init()");
      }
      const auto size = EVP_MD_size(md_);
      check_out_size(size_t(size), out);
      unsigned int len = 0;
      check_libssl_result(1, EVP_DigestFinal_ex(ctx_.get(), out.data(), &len));
      check_libssl_result(size, int(len));
      hash_reset();
    }

   private:
    const EVP_MD* md_{nullptr};
    EvpMdCtxUniquePtr ctx_{nullptr};

    void hash_reset() noexcept {
      ctx_.reset(nullptr);
      md_ = nullptr;
    }

    void copy_impl(const Digest& that) {
      if (that.md_ != nullptr && that.ctx_ != nullptr) {
        hash_init(that.md_);
        check_libssl_result(1, EVP_MD_CTX_copy_ex(ctx_.get(), that.ctx_.get()));
      } else {
        this->hash_reset();
      }
    }

    void move_impl(Digest&& that) noexcept {
      std::swap(this->md_, that.md_);
      std::swap(this->ctx_, that.ctx_);
    }
  };

  static void hash(MutableByteRange out, const EVP_MD* md, ByteRange data) {
    Digest hash;
    hash.hash_init(md);
    hash.hash_update(data);
    hash.hash_final(out);
  }
  static void hash(MutableByteRange out, const EVP_MD* md, const IOBuf& data) {
    Digest hash;
    hash.hash_init(md);
    hash.hash_update(data);
    hash.hash_final(out);
  }
  static void sha1(MutableByteRange out, ByteRange data) {
    hash(out, EVP_sha1(), data);
  }
  static void sha1(MutableByteRange out, const IOBuf& data) {
    hash(out, EVP_sha1(), data);
  }
  static void sha256(MutableByteRange out, ByteRange data) {
    hash(out, EVP_sha256(), data);
  }
  static void sha256(MutableByteRange out, const IOBuf& data) {
    hash(out, EVP_sha256(), data);
  }

  class Hmac {
   public:
    Hmac() : ctx_(HMAC_CTX_new()) {}

    void hash_init(const EVP_MD* md, ByteRange key) {
      md_ = md;
      check_libssl_result(
          1,
          HMAC_Init_ex(ctx_.get(), key.data(), int(key.size()), md_, nullptr));
    }
    void hash_update(ByteRange data) {
      check_libssl_result(1, HMAC_Update(ctx_.get(), data.data(), data.size()));
    }
    void hash_update(const IOBuf& data) {
      for (auto r : data) {
        hash_update(r);
      }
    }
    void hash_final(MutableByteRange out) {
      const auto size = EVP_MD_size(md_);
      check_out_size(size_t(size), out);
      unsigned int len = 0;
      check_libssl_result(1, HMAC_Final(ctx_.get(), out.data(), &len));
      check_libssl_result(size, int(len));
      md_ = nullptr;
    }

   private:
    const EVP_MD* md_ = nullptr;
    HmacCtxUniquePtr ctx_{nullptr};
  };

  static void hmac(
      MutableByteRange out, const EVP_MD* md, ByteRange key, ByteRange data) {
    Hmac hmac;
    hmac.hash_init(md, key);
    hmac.hash_update(data);
    hmac.hash_final(out);
  }
  static void hmac(
      MutableByteRange out,
      const EVP_MD* md,
      ByteRange key,
      const IOBuf& data) {
    Hmac hmac;
    hmac.hash_init(md, key);
    hmac.hash_update(data);
    hmac.hash_final(out);
  }
  static void hmac_sha1(MutableByteRange out, ByteRange key, ByteRange data) {
    hmac(out, EVP_sha1(), key, data);
  }
  static void hmac_sha1(
      MutableByteRange out, ByteRange key, const IOBuf& data) {
    hmac(out, EVP_sha1(), key, data);
  }
  static void hmac_sha256(MutableByteRange out, ByteRange key, ByteRange data) {
    hmac(out, EVP_sha256(), key, data);
  }
  static void hmac_sha256(
      MutableByteRange out, ByteRange key, const IOBuf& data) {
    hmac(out, EVP_sha256(), key, data);
  }

 private:
  static inline void check_out_size(size_t size, MutableByteRange out) {
    if (LIKELY(size == out.size())) {
      return;
    }
    check_out_size_throw(size, out);
  }
  [[noreturn]] static void check_out_size_throw(
      size_t size, MutableByteRange out);

  static inline void check_libssl_result(int expected, int result) {
    if (LIKELY(result == expected)) {
      return;
    }
    throw_exception<std::runtime_error>("openssl crypto function failed");
  }
};

} // namespace ssl
} // namespace folly
