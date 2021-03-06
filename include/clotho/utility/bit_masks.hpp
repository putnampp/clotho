//   Copyright 2015 Patrick Putnam
//
//   Licensed under the Apache License, Version 2.0 (the "License");
//   you may not use this file except in compliance with the License.
//   You may obtain a copy of the License at
//
//       http://www.apache.org/licenses/LICENSE-2.0
//
//   Unless required by applicable law or agreed to in writing, software
//   distributed under the License is distributed on an "AS IS" BASIS,
//   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//   See the License for the specific language governing permissions and
//   limitations under the License.
#ifndef BIT_MASK_HPP_
#define BIT_MASK_HPP_

#ifndef POWER2S
#define POWER2S

#define POW2_0 0x0000000000000001
#define POW2_1 0x0000000000000002
#define POW2_2 0x0000000000000004
#define POW2_3 0x0000000000000008
#define POW2_4 0x0000000000000010
#define POW2_5 0x0000000000000020
#define POW2_6 0x0000000000000040
#define POW2_7 0x0000000000000080
#define POW2_8 0x0000000000000100
#define POW2_9 0x0000000000000200
#define POW2_10 0x0000000000000400
#define POW2_11 0x0000000000000800
#define POW2_12 0x0000000000001000
#define POW2_13 0x0000000000002000
#define POW2_14 0x0000000000004000
#define POW2_15 0x0000000000008000
#define POW2_16 0x0000000000010000
#define POW2_17 0x0000000000020000
#define POW2_18 0x0000000000040000
#define POW2_19 0x0000000000080000
#define POW2_20 0x0000000000100000
#define POW2_21 0x0000000000200000
#define POW2_22 0x0000000000400000
#define POW2_23 0x0000000000800000
#define POW2_24 0x0000000001000000
#define POW2_25 0x0000000002000000
#define POW2_26 0x0000000004000000
#define POW2_27 0x0000000008000000
#define POW2_28 0x0000000010000000
#define POW2_29 0x0000000020000000
#define POW2_30 0x0000000040000000
#define POW2_31 0x0000000080000000
#define POW2_32 0x0000000100000000
#define POW2_33 0x0000000200000000
#define POW2_34 0x0000000400000000
#define POW2_35 0x0000000800000000
#define POW2_36 0x0000001000000000
#define POW2_37 0x0000002000000000
#define POW2_38 0x0000004000000000
#define POW2_39 0x0000008000000000
#define POW2_40 0x0000010000000000
#define POW2_41 0x0000020000000000
#define POW2_42 0x0000040000000000
#define POW2_43 0x0000080000000000
#define POW2_44 0x0000100000000000
#define POW2_45 0x0000200000000000
#define POW2_46 0x0000400000000000
#define POW2_47 0x0000800000000000
#define POW2_48 0x0001000000000000
#define POW2_49 0x0002000000000000
#define POW2_50 0x0004000000000000
#define POW2_51 0x0008000000000000
#define POW2_52 0x0010000000000000
#define POW2_53 0x0020000000000000
#define POW2_54 0x0040000000000000
#define POW2_55 0x0080000000000000
#define POW2_56 0x0100000000000000
#define POW2_57 0x0200000000000000
#define POW2_58 0x0400000000000000
#define POW2_59 0x0800000000000000
#define POW2_60 0x1000000000000000
#define POW2_61 0x2000000000000000
#define POW2_62 0x4000000000000000
#define POW2_63 0x8000000000000000

#define LOW_POW2_0 0x0000000000000001
#define LOW_POW2_1 0x0000000000000003
#define LOW_POW2_2 0x0000000000000007
#define LOW_POW2_3 0x000000000000000F
#define LOW_POW2_4 0x000000000000001F
#define LOW_POW2_5 0x000000000000003F
#define LOW_POW2_6 0x000000000000007F
#define LOW_POW2_7 0x00000000000000FF
#define LOW_POW2_8 0x00000000000001FF
#define LOW_POW2_9 0x00000000000003FF
#define LOW_POW2_10 0x00000000000007FF
#define LOW_POW2_11 0x0000000000000FFF
#define LOW_POW2_12 0x0000000000001FFF
#define LOW_POW2_13 0x0000000000003FFF
#define LOW_POW2_14 0x0000000000007FFF
#define LOW_POW2_15 0x000000000000FFFF
#define LOW_POW2_16 0x000000000001FFFF
#define LOW_POW2_17 0x000000000003FFFF
#define LOW_POW2_18 0x000000000007FFFF
#define LOW_POW2_19 0x00000000000FFFFF
#define LOW_POW2_20 0x00000000001FFFFF
#define LOW_POW2_21 0x00000000003FFFFF
#define LOW_POW2_22 0x00000000007FFFFF
#define LOW_POW2_23 0x0000000000FFFFFF
#define LOW_POW2_24 0x0000000001FFFFFF
#define LOW_POW2_25 0x0000000003FFFFFF
#define LOW_POW2_26 0x0000000007FFFFFF
#define LOW_POW2_27 0x000000000FFFFFFF
#define LOW_POW2_28 0x000000001FFFFFFF
#define LOW_POW2_29 0x000000003FFFFFFF
#define LOW_POW2_30 0x000000007FFFFFFF
#define LOW_POW2_31 0x00000000FFFFFFFF
#define LOW_POW2_32 0x00000001FFFFFFFF
#define LOW_POW2_33 0x00000003FFFFFFFF
#define LOW_POW2_34 0x00000007FFFFFFFF
#define LOW_POW2_35 0x0000000FFFFFFFFF
#define LOW_POW2_36 0x0000001FFFFFFFFF
#define LOW_POW2_37 0x0000003FFFFFFFFF
#define LOW_POW2_38 0x0000007FFFFFFFFF
#define LOW_POW2_39 0x000000FFFFFFFFFF
#define LOW_POW2_40 0x000001FFFFFFFFFF
#define LOW_POW2_41 0x000003FFFFFFFFFF
#define LOW_POW2_42 0x000007FFFFFFFFFF
#define LOW_POW2_43 0x00000FFFFFFFFFFF
#define LOW_POW2_44 0x00001FFFFFFFFFFF
#define LOW_POW2_45 0x00003FFFFFFFFFFF
#define LOW_POW2_46 0x00007FFFFFFFFFFF
#define LOW_POW2_47 0x0000FFFFFFFFFFFF
#define LOW_POW2_48 0x0001FFFFFFFFFFFF
#define LOW_POW2_49 0x0003FFFFFFFFFFFF
#define LOW_POW2_50 0x0007FFFFFFFFFFFF
#define LOW_POW2_51 0x000FFFFFFFFFFFFF
#define LOW_POW2_52 0x001FFFFFFFFFFFFF
#define LOW_POW2_53 0x003FFFFFFFFFFFFF
#define LOW_POW2_54 0x007FFFFFFFFFFFFF
#define LOW_POW2_55 0x00FFFFFFFFFFFFFF
#define LOW_POW2_56 0x01FFFFFFFFFFFFFF
#define LOW_POW2_57 0x03FFFFFFFFFFFFFF
#define LOW_POW2_58 0x07FFFFFFFFFFFFFF
#define LOW_POW2_59 0x0FFFFFFFFFFFFFFF
#define LOW_POW2_60 0x1FFFFFFFFFFFFFFF
#define LOW_POW2_61 0x3FFFFFFFFFFFFFFF
#define LOW_POW2_62 0x7FFFFFFFFFFFFFFF
#define LOW_POW2_63 0xFFFFFFFFFFFFFFFF

#define CAT( x, y ) x ## y
#define OFFSET( x ) CAT(POW2_, x )

#define DEBRUIJN32 0x077CB531U
// hashed bit offsets for bits [0,31]
static const unsigned int deBruijn_bit_position[ 32 ] = {
    0, 1, 28, 2, 29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4, 8,
    31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6, 11, 5, 10, 9
};

// hashed bit offsets for bits [32,63]
static const unsigned int deBruijn_bit_position_hi[ 32 ] = {
    32, 33, 60, 34, 61, 46, 56, 35, 62, 54, 52, 47, 57, 49, 36, 40,
    63, 59, 45, 55, 53, 51, 48, 39, 58, 44, 50, 38, 43, 37, 42, 41
};

#define LEAST_SIG_BIT( x ) (x & (-x))

#define DEBRUIJNBIT_HASH( LSB ) ((LSB * DEBRUIJN32) >> 27)

#define DEBRUIJNBIT_HASH_LOOKUP( LSB ) (deBruijn_bit_position[ DEBRUIJNBIT_HASH( LSB ) ])

#define DEBRUIJNBIT_HASH_LOOKUP_HI( LSB ) (deBruijn_bit_position_hi[ DEBRUIJNBIT_HASH( LSB ) ])

#endif  // POWER2S

namespace clotho {
namespace utility {

// 64 * 8 = 512 bytes
static const unsigned long bit_masks[ 64 ] = {
    POW2_0, POW2_1, POW2_2, POW2_3,
    POW2_4, POW2_5, POW2_6, POW2_7,
    POW2_8, POW2_9, POW2_10, POW2_11,
    POW2_12, POW2_13, POW2_14, POW2_15,
    POW2_16, POW2_17, POW2_18, POW2_19,
    POW2_20, POW2_21, POW2_22, POW2_23,
    POW2_24, POW2_25, POW2_26, POW2_27,
    POW2_28, POW2_29, POW2_30, POW2_31,
    POW2_32, POW2_33, POW2_34, POW2_35,
    POW2_36, POW2_37, POW2_38, POW2_39,
    POW2_40, POW2_41, POW2_42, POW2_43,
    POW2_44, POW2_45, POW2_46, POW2_47,
    POW2_48, POW2_49, POW2_50, POW2_51,
    POW2_52, POW2_53, POW2_54, POW2_55,
    POW2_56, POW2_57, POW2_58, POW2_59,
    POW2_60, POW2_61, POW2_62, POW2_63
};

static const unsigned long low_bit_masks[ 64 ] = {
    LOW_POW2_0, LOW_POW2_1, LOW_POW2_2, LOW_POW2_3,
    LOW_POW2_4, LOW_POW2_5, LOW_POW2_6, LOW_POW2_7,
    LOW_POW2_8, LOW_POW2_9, LOW_POW2_10, LOW_POW2_11,
    LOW_POW2_12, LOW_POW2_13, LOW_POW2_14, LOW_POW2_15,
    LOW_POW2_16, LOW_POW2_17, LOW_POW2_18, LOW_POW2_19,
    LOW_POW2_20, LOW_POW2_21, LOW_POW2_22, LOW_POW2_23,
    LOW_POW2_24, LOW_POW2_25, LOW_POW2_26, LOW_POW2_27,
    LOW_POW2_28, LOW_POW2_29, LOW_POW2_30, LOW_POW2_31,
    LOW_POW2_32, LOW_POW2_33, LOW_POW2_34, LOW_POW2_35,
    LOW_POW2_36, LOW_POW2_37, LOW_POW2_38, LOW_POW2_39,
    LOW_POW2_40, LOW_POW2_41, LOW_POW2_42, LOW_POW2_43,
    LOW_POW2_44, LOW_POW2_45, LOW_POW2_46, LOW_POW2_47,
    LOW_POW2_48, LOW_POW2_49, LOW_POW2_50, LOW_POW2_51,
    LOW_POW2_52, LOW_POW2_53, LOW_POW2_54, LOW_POW2_55,
    LOW_POW2_56, LOW_POW2_57, LOW_POW2_58, LOW_POW2_59,
    LOW_POW2_60, LOW_POW2_61, LOW_POW2_62, LOW_POW2_63
};

unsigned int scan_set_bits( unsigned int b, unsigned int * indices );
unsigned int scan_set_bits( unsigned long b, unsigned int * indices );

unsigned int hash_set_bits( unsigned int b, unsigned int * indices );
unsigned int hash_set_bits( unsigned long b, unsigned int * indices );

}   // namespace utility
}   // namespace clotho

#endif  // BIT_MASK_HPP_
