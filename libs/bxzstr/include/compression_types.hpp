/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *
 * This file is a part of bxzstr (https://github.com/tmaklin/bxzstr)
 * Written by Tommi Mäklin (tommi@maklin.fi) */

#ifndef BXZSTR_COMPRESSION_TYPES_HPP
#define BXZSTR_COMPRESSION_TYPES_HPP

#include <exception>

#include "stream_wrapper.hpp"
#include "bz_stream_wrapper.hpp"
#include "lzma_stream_wrapper.hpp"
#include "z_stream_wrapper.hpp"
#include "zstd_stream_wrapper.hpp"

namespace bxz {
    enum Compression { z, bz2, lzma, zstd, plaintext };
inline Compression detect_type(char* in_buff_start, char* in_buff_end) {
    unsigned char b0 = *reinterpret_cast< unsigned char * >(in_buff_start);
    unsigned char b1 = *reinterpret_cast< unsigned char * >(in_buff_start + 1);
#if defined(BXZSTR_Z_STREAM_WRAPPER_HPP) || defined(BXZSTR_BZ_STREAM_WRAPPER_HPP) || defined(BXZSTR_LZMA_STREAM_WRAPPER_HPP) || defined(BXZSTR_ZSTD_STREAM_WRAPPER_HPP)
    bool gzip_header = (b0 == 0x1F && b1 == 0x8B);
    bool zlib_header = (b0 == 0x78 && (b1 == 0x01 || b1 == 0x9C || b1 == 0xDA));
    if (in_buff_start + 1 <= in_buff_end && (gzip_header || zlib_header)) return z;
#endif
#if defined(BXZSTR_BZ_STREAM_WRAPPER_HPP) || defined(BXZSTR_LZMA_STREAM_WRAPPER_HPP) || defined(BXZSTR_ZSTD_STREAM_WRAPPER_HPP)
    unsigned char b2 = *reinterpret_cast< unsigned char * >(in_buff_start + 2);
    bool bz2_header = (b0 == 0x42 && b1 == 0x5a && b2 == 0x68);
    if (in_buff_start + 2 <= in_buff_end && bz2_header) return bz2;
#endif
#if defined(BXZSTR_LZMA_STREAM_WRAPPER_HPP) || defined(BXZSTR_ZSTD_STREAM_WRAPPER_HPP)
    unsigned char b3 = *reinterpret_cast< unsigned char * >(in_buff_start + 3);
    unsigned char b4 = *reinterpret_cast< unsigned char * >(in_buff_start + 4);
    unsigned char b5 = *reinterpret_cast< unsigned char * >(in_buff_start + 5);
    bool lzma_header = (b0 == 0xFD && b1 == 0x37 && b2 == 0x7A
			&& b3 == 0x58 && b4 == 0x5A && b5 == 0x00);
    if (in_buff_start + 5 <= in_buff_end && lzma_header) return lzma;
#endif
#if defined(BXZSTR_ZSTD_STREAM_WRAPPER_HPP)
    bool zstd_header = (b0 == 0x28 && b1 == 0xB5 && b2 == 0x2F && b3 == 0xFD);
    if (in_buff_start + 3 <= in_buff_end && zstd_header) return zstd;
#endif
    return plaintext;
}

inline void init_stream(const Compression &type, const bool is_input, const int level,
			std::unique_ptr<detail::stream_wrapper> *strm_p) {
    switch (type) {
#ifdef BXZSTR_LZMA_STREAM_WRAPPER_HPP
        case lzma : strm_p->reset(new detail::lzma_stream_wrapper(is_input, level));
	break;
#endif
#ifdef BXZSTR_BZ_STREAM_WRAPPER_HPP
        case bz2 : strm_p->reset(new detail::bz_stream_wrapper(is_input, level));
	break;
#endif
#ifdef BXZSTR_Z_STREAM_WRAPPER_HPP
        case z : strm_p->reset(new detail::z_stream_wrapper(is_input, level));
	break;
#endif
#ifdef BXZSTR_ZSTD_STREAM_WRAPPER_HPP
        case zstd : strm_p->reset(new detail::zstd_stream_wrapper(is_input, level));
	break;
#endif
	default : throw std::runtime_error("Unrecognized compression type.");
    }
}
inline void init_stream(const Compression &type, const bool is_input,
			std::unique_ptr<detail::stream_wrapper> *strm_p) {
    init_stream(type, is_input, 6, strm_p);
}

inline int bxz_run(const Compression &type) {
    switch(type){
#ifdef BXZSTR_LZMA_STREAM_WRAPPER_HPP
	case lzma: return 0;
	break; // LZMA_RUN
#endif
#ifdef BXZSTR_BZ_STREAM_WRAPPER_HPP
        case bz2: return 0;
	break; // BZ_RUN
#endif
#ifdef BXZSTR_Z_STREAM_WRAPPER_HPP
        case z: return 0;
	break;// Z_NO_FLUSH
#endif
#ifdef BXZSTR_ZSTD_STREAM_WRAPPER_HPP
        case zstd: return 0;
	break;// ZSTD_NO_FLUSH
#endif
	default: throw std::runtime_error("Unrecognized compression type.");
    }
}
inline int bxz_finish(const Compression &type) {
    switch(type){
#ifdef BXZSTR_LZMA_STREAM_WRAPPER_HPP
	case lzma: return 3;
	break; // LZMA_FINISH
#endif
#ifdef BXZSTR_BZ_STREAM_WRAPPER_HPP
        case bz2: return 2;
	break; // BZ_FINISH
#endif
#ifdef BXZSTR_Z_STREAM_WRAPPER_HPP
        case z: return 4;
	break; // Z_FINISH
#endif
#ifdef BXZSTR_ZSTD_STREAM_WRAPPER_HPP
        case zstd: return 1;
	break; // endStream == true
#endif
	default: throw std::runtime_error("Unrecognized compression type.");
    }
}
}

#endif
