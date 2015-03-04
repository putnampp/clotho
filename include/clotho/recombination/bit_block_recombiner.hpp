#ifndef BIT_BLOCK_RECOMBINER_HPP_
#define BIT_BLOCK_RECOMBINER_HPP_

#include "clotho/utility/bit_block_iterator.hpp"
#include "clotho/utility/bit_fill.hpp"

namespace clotho {
namespace recombine {

struct copy_matching_classify_mismatch {};

template < class Classifier, class Tag = copy_matching_classify_mismatch >
class bit_block_recombiner {
public:
    typedef Classifier classifier_type;

    bit_block_recombiner( const classifier_type & cfier ) : m_cfier( cfier ) {}

    template < class Block, class ElementIterator >
    Block operator()( Block b0, Block b1, ElementIterator first ) {
        typedef Block block_type;
//        typedef clotho::utility::bit_block_iterator< block_type > iterator;
/*
        // b0 = 01001
        // b1 = 10101
        block_type matching = (b0 & b1); // = 00001
        block_type mismatch = (b0 ^ b1); // = 11100

        block_type b0_mask = (block_type)0;
        if( mismatch ) {
            iterator bit_it( mismatch ), bit_end;
            do { 
                unsigned int bit_idx = (*bit_it);
                if( m_cfier(first, bit_idx) ) {
                    b0_mask |= ((block_type)1 << (bit_idx));
                }
            } while( ++bit_it != bit_end);
            // assume:
            // classify( bit_idx = 2 ) true
            // classify( bit_idx = 3 ) true
            // classify( bit_idx = 4 ) false
            // => b0_mask = 01100
        }

        // ~b0_mask & mismatch = ~01100 & 11100 = 10011 & 11100 = 10000
        // b1_mismatch = (10000 & 10101) = 10000
        // res = (matching | b1_mismatch) = (00001 | 10000) = 10001
        // res = res | (b0_mask & b0) = (10001 | (01100 & 01001)) = (10001 | (01000)) = 11001
        block_type b1_mismatch = ((~b0_mask & mismatch) & b1);
        block_type res = (matching | b1_mismatch);
        res |= (b0_mask & b0);
        return  res;
*/
        // Simplified Logic:
        // b0 = 01001
        // b1 = 10101
        block_type hets = (b0 ^ b1);    // = 11100
        block_type b0_mask = (block_type)0;
        if( hets ) {
//            iterator bit_it( hets ), bit_end;
//            do {
//               unsigned int bit_idx = (*bit_it);
//                if( m_cfier(first, bit_idx ) ) {
//                    b0_mask |= ((block_type)1 << bit_idx);
//                }
//            } while( ++bit_it != bit_end );

//            unsigned int indices[ sizeof( block_type ) * 8];

            // preprocess set bits 
            //
//            unsigned int max_index = fill_set_bits( hets, indices );
//
//            for( unsigned int offset = 0; offset < max_index; ++offset ) {
//                unsigned int idx = indices[offset];
//                if( m_cfier( *(first + idx ) ) ) {
//                    b0_mask |= ((block_type) 1 << idx);
//                }
//            }
            // assume:
            // classify( bit_idx = 2 ) true
            // classify( bit_idx = 3 ) true
            // classify( bit_idx = 4 ) false
            // => b0_mask = 01100

            b0_mask = block_logic256( hets, first );
        }
        // b0 & b0_mask = (01001  &  01100) = 01000
        // b1 & ~b0_mask = (10101 & ~01100) = (10101 & 10011) = 10001
        // res = (01000 | 10001) = 11001
        block_type res = ((b0 & b0_mask) | (b1 & ~b0_mask));

        return res;
    }

    virtual ~bit_block_recombiner() {}

protected:
    classifier_type m_cfier;

#define BIT_SET_EVEN_1( z )                     \
    first += z;                                 \
    if( m_cfier( *first ) ) { res |= (m << z); }\
    first += (8 - z);                           \
    break;

#define BIT_SET_ODD_2( z )                          \
    if( m_cfier( *first ) ) { res |= m; }           \
    first += z;                                     \
    if( m_cfier( *first ) ) { res |= (m << z); }    \
    first += (8 - z);                               \
    break;

#define BIT_SET_EVEN_2( y, z )                      \
    first += y;                                     \
    if( m_cfier( *first ) ) { res |= (m << y); }    \
    first += (z - y);                               \
    if( m_cfier( *first ) ) { res |= (m << z); }    \
    first += (8 - z);                               \
    break;

#define BIT_SET_ODD_3( y, z )                       \
    if( m_cfier( *first ) ) { res |= m; }           \
    first += y;                                     \
    if( m_cfier( *first ) ) { res |= (m << y); }    \
    first += (z - y);                               \
    if( m_cfier( *first ) ) { res |= (m << z); }    \
    first += (8 - z);                               \
    break;

#define BIT_SET_EVEN_3( x, y, z )                   \
    first += x;                                     \
    if( m_cfier( *first ) ) { res |= (m << x); }    \
    first += (y - x);                               \
    if( m_cfier( *first ) ) { res |= (m << y); }    \
    first += (z - y);                               \
    if( m_cfier( *first ) ) { res |= (m << z); }    \
    first += (8 - z);                               \
    break;

#define BIT_SET_ODD_4( x, y, z )                    \
    if( m_cfier( *first ) ) { res |= m; }           \
    first += x;                                     \
    if( m_cfier( *first ) ) { res |= (m << x); }    \
    first += (y - x);                               \
    if( m_cfier( *first ) ) { res |= (m << y); }    \
    first += (z - y);                               \
    if( m_cfier( *first ) ) { res |= (m << z); }    \
    first += (8 - z);                               \
    break;


#define BIT_SET_EVEN_4( w, x, y, z )                \
    first += w;                                     \
    if( m_cfier( *first ) ) { res |= (m << w); }    \
    first += (x - w);                               \
    if( m_cfier( *first ) ) { res |= (m << x); }    \
    first += (y - x);                               \
    if( m_cfier( *first ) ) { res |= (m << y); }    \
    first += (z - y);                               \
    if( m_cfier( *first ) ) { res |= (m << z); }    \
    first += (8 - z);                               \
    break;

#define BIT_SET_ODD_5( w, x, y, z )                 \
    if( m_cfier( *first ) ) { res |= m; }           \
    first += w;                                     \
    if( m_cfier( *first ) ) { res |= (m << w); }    \
    first += (x - w);                               \
    if( m_cfier( *first ) ) { res |= (m << x); }    \
    first += (y - x);                               \
    if( m_cfier( *first ) ) { res |= (m << y); }    \
    first += (z - y);                               \
    if( m_cfier( *first ) ) { res |= (m << z); }    \
    first += (8 - z);                               \
    break;

#define BIT_SET_EVEN_5(v, w, x, y, z )              \
    first += v;                                     \
    if( m_cfier( *first ) ) { res |= (m << v); }    \
    first += (w - v);                               \
    if( m_cfier( *first ) ) { res |= (m << w); }    \
    first += (x - w);                               \
    if( m_cfier( *first ) ) { res |= (m << x); }    \
    first += (y - x);                               \
    if( m_cfier( *first ) ) { res |= (m << y); }    \
    first += (z - y);                               \
    if( m_cfier( *first ) ) { res |= (m << z); }    \
    first += (8 - z);                               \
    break;

#define BIT_SET_ODD_6( v, w, x, y, z )              \
    if( m_cfier( *first ) ) { res |= m; }           \
    first += v;                                     \
    if( m_cfier( *first ) ) { res |= (m << v); }    \
    first += (w - v);                               \
    if( m_cfier( *first ) ) { res |= (m << w); }    \
    first += (x - w);                               \
    if( m_cfier( *first ) ) { res |= (m << x); }    \
    first += (y - x);                               \
    if( m_cfier( *first ) ) { res |= (m << y); }    \
    first += (z - y);                               \
    if( m_cfier( *first ) ) { res |= (m << z); }    \
    first += (8 - z);                               \
    break;

#define BIT_SET_EVEN_6(u, v, w, x, y, z )           \
    first += u;                                     \
    if( m_cfier( *first ) ) { res |= (m << u); }    \
    first += (v - u);                               \
    if( m_cfier( *first ) ) { res |= (m << v); }    \
    first += (w - v);                               \
    if( m_cfier( *first ) ) { res |= (m << w); }    \
    first += (x - w);                               \
    if( m_cfier( *first ) ) { res |= (m << x); }    \
    first += (y - x);                               \
    if( m_cfier( *first ) ) { res |= (m << y); }    \
    first += (z - y);                               \
    if( m_cfier( *first ) ) { res |= (m << z); }    \
    first += (8 - z);                               \
    break;

#define BIT_SET_ODD_7( u, v, w, x, y, z )           \
    if( m_cfier( *first ) ) { res |= m; }           \
    first += u;                                     \
    if( m_cfier( *first ) ) { res |= (m << u); }    \
    first += (v - u);                               \
    if( m_cfier( *first ) ) { res |= (m << v); }    \
    first += (w - v);                               \
    if( m_cfier( *first ) ) { res |= (m << w); }    \
    first += (x - w);                               \
    if( m_cfier( *first ) ) { res |= (m << x); }    \
    first += (y - x);                               \
    if( m_cfier( *first ) ) { res |= (m << y); }    \
    first += (z - y);                               \
    if( m_cfier( *first ) ) { res |= (m << z); }    \
    first += (8 - z);                               \
    break;

#define BIT_SET_EVEN_7(t, u, v, w, x, y, z )        \
    first += t;                                     \
    if( m_cfier( *first ) ) { res |= (m << t); }    \
    first += (u - t);                               \
    if( m_cfier( *first ) ) { res |= (m << u); }    \
    first += (v - u);                               \
    if( m_cfier( *first ) ) { res |= (m << v); }    \
    first += (w - v);                               \
    if( m_cfier( *first ) ) { res |= (m << w); }    \
    first += (x - w);                               \
    if( m_cfier( *first ) ) { res |= (m << x); }    \
    first += (y - x);                               \
    if( m_cfier( *first ) ) { res |= (m << y); }    \
    first += (z - y);                               \
    if( m_cfier( *first ) ) { res |= (m << z); }    \
    first += (8 - z);                               \
    break;

    template < class Block, class ElementIterator >
    inline Block block_logic256( Block b, ElementIterator first ) {
        Block res = (Block)0, m = (Block)1;
        while( b ) {
            switch( b & 255 ) {
            case 0: 
                first += 8;
                break;
            case 1:
                if( m_cfier( *first ) ) { res |= m; }
                first += 8;
                break;

            case 2:  BIT_SET_EVEN_1(1)
            case 3:  BIT_SET_ODD_2(1)
            case 4:  BIT_SET_EVEN_1(2)
            case 5:  BIT_SET_ODD_2(2)
            case 6:  BIT_SET_EVEN_2(1,2)
            case 7:  BIT_SET_ODD_3(1,2)
            case 8:  BIT_SET_EVEN_1(3)
            case 9:  BIT_SET_ODD_2(3)
            case 10:  BIT_SET_EVEN_2(1,3)
            case 11:  BIT_SET_ODD_3(1,3)
            case 12:  BIT_SET_EVEN_2(2,3)
            case 13:  BIT_SET_ODD_3(2,3)
            case 14:  BIT_SET_EVEN_3(1,2,3)
            case 15:  BIT_SET_ODD_4(1,2,3)
            case 16:  BIT_SET_EVEN_1(4)
            case 17:  BIT_SET_ODD_2(4)
            case 18:  BIT_SET_EVEN_2(1,4)
            case 19:  BIT_SET_ODD_3(1,4)
            case 20:  BIT_SET_EVEN_2(2,4)
            case 21:  BIT_SET_ODD_3(2,4)
            case 22:  BIT_SET_EVEN_3(1,2,4)
            case 23:  BIT_SET_ODD_4(1,2,4)
            case 24:  BIT_SET_EVEN_2(3,4)
            case 25:  BIT_SET_ODD_3(3,4)
            case 26:  BIT_SET_EVEN_3(1,3,4)
            case 27:  BIT_SET_ODD_4(1,3,4)
            case 28:  BIT_SET_EVEN_3(2,3,4)
            case 29:  BIT_SET_ODD_4(2,3,4)
            case 30:  BIT_SET_EVEN_4(1,2,3,4)
            case 31:  BIT_SET_ODD_5(1,2,3,4)
            case 32:  BIT_SET_EVEN_1(5)
            case 33:  BIT_SET_ODD_2(5)
            case 34:  BIT_SET_EVEN_2(1,5)
            case 35:  BIT_SET_ODD_3(1,5)
            case 36:  BIT_SET_EVEN_2(2,5)
            case 37:  BIT_SET_ODD_3(2,5)
            case 38:  BIT_SET_EVEN_3(1,2,5)
            case 39:  BIT_SET_ODD_4(1,2,5)
            case 40:  BIT_SET_EVEN_2(3,5)
            case 41:  BIT_SET_ODD_3(3,5)
            case 42:  BIT_SET_EVEN_3(1,3,5)
            case 43:  BIT_SET_ODD_4(1,3,5)
            case 44:  BIT_SET_EVEN_3(2,3,5)
            case 45:  BIT_SET_ODD_4(2,3,5)
            case 46:  BIT_SET_EVEN_4(1,2,3,5)
            case 47:  BIT_SET_ODD_5(1,2,3,5)
            case 48:  BIT_SET_EVEN_2(4,5)
            case 49:  BIT_SET_ODD_3(4,5)
            case 50:  BIT_SET_EVEN_3(1,4,5)
            case 51:  BIT_SET_ODD_4(1,4,5)
            case 52:  BIT_SET_EVEN_3(2,4,5)
            case 53:  BIT_SET_ODD_4(2,4,5)
            case 54:  BIT_SET_EVEN_4(1,2,4,5)
            case 55:  BIT_SET_ODD_5(1,2,4,5)
            case 56:  BIT_SET_EVEN_3(3,4,5)
            case 57:  BIT_SET_ODD_4(3,4,5)
            case 58:  BIT_SET_EVEN_4(1,3,4,5)
            case 59:  BIT_SET_ODD_5(1,3,4,5)
            case 60:  BIT_SET_EVEN_4(2,3,4,5)
            case 61:  BIT_SET_ODD_5(2,3,4,5)
            case 62:  BIT_SET_EVEN_5(1,2,3,4,5)
            case 63:  BIT_SET_ODD_6(1,2,3,4,5)
            case 64:  BIT_SET_EVEN_1(6)
            case 65:  BIT_SET_ODD_2(6)
            case 66:  BIT_SET_EVEN_2(1,6)
            case 67:  BIT_SET_ODD_3(1,6)
            case 68:  BIT_SET_EVEN_2(2,6)
            case 69:  BIT_SET_ODD_3(2,6)
            case 70:  BIT_SET_EVEN_3(1,2,6)
            case 71:  BIT_SET_ODD_4(1,2,6)
            case 72:  BIT_SET_EVEN_2(3,6)
            case 73:  BIT_SET_ODD_3(3,6)
            case 74:  BIT_SET_EVEN_3(1,3,6)
            case 75:  BIT_SET_ODD_4(1,3,6)
            case 76:  BIT_SET_EVEN_3(2,3,6)
            case 77:  BIT_SET_ODD_4(2,3,6)
            case 78:  BIT_SET_EVEN_4(1,2,3,6)
            case 79:  BIT_SET_ODD_5(1,2,3,6)
            case 80:  BIT_SET_EVEN_2(4,6)
            case 81:  BIT_SET_ODD_3(4,6)
            case 82:  BIT_SET_EVEN_3(1,4,6)
            case 83:  BIT_SET_ODD_4(1,4,6)
            case 84:  BIT_SET_EVEN_3(2,4,6)
            case 85:  BIT_SET_ODD_4(2,4,6)
            case 86:  BIT_SET_EVEN_4(1,2,4,6)
            case 87:  BIT_SET_ODD_5(1,2,4,6)
            case 88:  BIT_SET_EVEN_3(3,4,6)
            case 89:  BIT_SET_ODD_4(3,4,6)
            case 90:  BIT_SET_EVEN_4(1,3,4,6)
            case 91:  BIT_SET_ODD_5(1,3,4,6)
            case 92:  BIT_SET_EVEN_4(2,3,4,6)
            case 93:  BIT_SET_ODD_5(2,3,4,6)
            case 94:  BIT_SET_EVEN_5(1,2,3,4,6)
            case 95:  BIT_SET_ODD_6(1,2,3,4,6)
            case 96:  BIT_SET_EVEN_2(5,6)
            case 97:  BIT_SET_ODD_3(5,6)
            case 98:  BIT_SET_EVEN_3(1,5,6)
            case 99:  BIT_SET_ODD_4(1,5,6)
            case 100:  BIT_SET_EVEN_3(2,5,6)
            case 101:  BIT_SET_ODD_4(2,5,6)
            case 102:  BIT_SET_EVEN_4(1,2,5,6)
            case 103:  BIT_SET_ODD_5(1,2,5,6)
            case 104:  BIT_SET_EVEN_3(3,5,6)
            case 105:  BIT_SET_ODD_4(3,5,6)
            case 106:  BIT_SET_EVEN_4(1,3,5,6)
            case 107:  BIT_SET_ODD_5(1,3,5,6)
            case 108:  BIT_SET_EVEN_4(2,3,5,6)
            case 109:  BIT_SET_ODD_5(2,3,5,6)
            case 110:  BIT_SET_EVEN_5(1,2,3,5,6)
            case 111:  BIT_SET_ODD_6(1,2,3,5,6)
            case 112:  BIT_SET_EVEN_3(4,5,6)
            case 113:  BIT_SET_ODD_4(4,5,6)
            case 114:  BIT_SET_EVEN_4(1,4,5,6)
            case 115:  BIT_SET_ODD_5(1,4,5,6)
            case 116:  BIT_SET_EVEN_4(2,4,5,6)
            case 117:  BIT_SET_ODD_5(2,4,5,6)
            case 118:  BIT_SET_EVEN_5(1,2,4,5,6)
            case 119:  BIT_SET_ODD_6(1,2,4,5,6)
            case 120:  BIT_SET_EVEN_4(3,4,5,6)
            case 121:  BIT_SET_ODD_5(3,4,5,6)
            case 122:  BIT_SET_EVEN_5(1,3,4,5,6)
            case 123:  BIT_SET_ODD_6(1,3,4,5,6)
            case 124:  BIT_SET_EVEN_5(2,3,4,5,6)
            case 125:  BIT_SET_ODD_6(2,3,4,5,6)
            case 126:  BIT_SET_EVEN_6(1,2,3,4,5,6)
            case 127:  BIT_SET_ODD_7(1,2,3,4,5,6)
            case 128:  BIT_SET_EVEN_1(7)
            case 129:  BIT_SET_ODD_2(7)
            case 130:  BIT_SET_EVEN_2(1,7)
            case 131:  BIT_SET_ODD_3(1,7)
            case 132:  BIT_SET_EVEN_2(2,7)
            case 133:  BIT_SET_ODD_3(2,7)
            case 134:  BIT_SET_EVEN_3(1,2,7)
            case 135:  BIT_SET_ODD_4(1,2,7)
            case 136:  BIT_SET_EVEN_2(3,7)
            case 137:  BIT_SET_ODD_3(3,7)
            case 138:  BIT_SET_EVEN_3(1,3,7)
            case 139:  BIT_SET_ODD_4(1,3,7)
            case 140:  BIT_SET_EVEN_3(2,3,7)
            case 141:  BIT_SET_ODD_4(2,3,7)
            case 142:  BIT_SET_EVEN_4(1,2,3,7)
            case 143:  BIT_SET_ODD_5(1,2,3,7)
            case 144:  BIT_SET_EVEN_2(4,7)
            case 145:  BIT_SET_ODD_3(4,7)
            case 146:  BIT_SET_EVEN_3(1,4,7)
            case 147:  BIT_SET_ODD_4(1,4,7)
            case 148:  BIT_SET_EVEN_3(2,4,7)
            case 149:  BIT_SET_ODD_4(2,4,7)
            case 150:  BIT_SET_EVEN_4(1,2,4,7)
            case 151:  BIT_SET_ODD_5(1,2,4,7)
            case 152:  BIT_SET_EVEN_3(3,4,7)
            case 153:  BIT_SET_ODD_4(3,4,7)
            case 154:  BIT_SET_EVEN_4(1,3,4,7)
            case 155:  BIT_SET_ODD_5(1,3,4,7)
            case 156:  BIT_SET_EVEN_4(2,3,4,7)
            case 157:  BIT_SET_ODD_5(2,3,4,7)
            case 158:  BIT_SET_EVEN_5(1,2,3,4,7)
            case 159:  BIT_SET_ODD_6(1,2,3,4,7)
            case 160:  BIT_SET_EVEN_2(5,7)
            case 161:  BIT_SET_ODD_3(5,7)
            case 162:  BIT_SET_EVEN_3(1,5,7)
            case 163:  BIT_SET_ODD_4(1,5,7)
            case 164:  BIT_SET_EVEN_3(2,5,7)
            case 165:  BIT_SET_ODD_4(2,5,7)
            case 166:  BIT_SET_EVEN_4(1,2,5,7)
            case 167:  BIT_SET_ODD_5(1,2,5,7)
            case 168:  BIT_SET_EVEN_3(3,5,7)
            case 169:  BIT_SET_ODD_4(3,5,7)
            case 170:  BIT_SET_EVEN_4(1,3,5,7)
            case 171:  BIT_SET_ODD_5(1,3,5,7)
            case 172:  BIT_SET_EVEN_4(2,3,5,7)
            case 173:  BIT_SET_ODD_5(2,3,5,7)
            case 174:  BIT_SET_EVEN_5(1,2,3,5,7)
            case 175:  BIT_SET_ODD_6(1,2,3,5,7)
            case 176:  BIT_SET_EVEN_3(4,5,7)
            case 177:  BIT_SET_ODD_4(4,5,7)
            case 178:  BIT_SET_EVEN_4(1,4,5,7)
            case 179:  BIT_SET_ODD_5(1,4,5,7)
            case 180:  BIT_SET_EVEN_4(2,4,5,7)
            case 181:  BIT_SET_ODD_5(2,4,5,7)
            case 182:  BIT_SET_EVEN_5(1,2,4,5,7)
            case 183:  BIT_SET_ODD_6(1,2,4,5,7)
            case 184:  BIT_SET_EVEN_4(3,4,5,7)
            case 185:  BIT_SET_ODD_5(3,4,5,7)
            case 186:  BIT_SET_EVEN_5(1,3,4,5,7)
            case 187:  BIT_SET_ODD_6(1,3,4,5,7)
            case 188:  BIT_SET_EVEN_5(2,3,4,5,7)
            case 189:  BIT_SET_ODD_6(2,3,4,5,7)
            case 190:  BIT_SET_EVEN_6(1,2,3,4,5,7)
            case 191:  BIT_SET_ODD_7(1,2,3,4,5,7)
            case 192:  BIT_SET_EVEN_2(6,7)
            case 193:  BIT_SET_ODD_3(6,7)
            case 194:  BIT_SET_EVEN_3(1,6,7)
            case 195:  BIT_SET_ODD_4(1,6,7)
            case 196:  BIT_SET_EVEN_3(2,6,7)
            case 197:  BIT_SET_ODD_4(2,6,7)
            case 198:  BIT_SET_EVEN_4(1,2,6,7)
            case 199:  BIT_SET_ODD_5(1,2,6,7)
            case 200:  BIT_SET_EVEN_3(3,6,7)
            case 201:  BIT_SET_ODD_4(3,6,7)
            case 202:  BIT_SET_EVEN_4(1,3,6,7)
            case 203:  BIT_SET_ODD_5(1,3,6,7)
            case 204:  BIT_SET_EVEN_4(2,3,6,7)
            case 205:  BIT_SET_ODD_5(2,3,6,7)
            case 206:  BIT_SET_EVEN_5(1,2,3,6,7)
            case 207:  BIT_SET_ODD_6(1,2,3,6,7)
            case 208:  BIT_SET_EVEN_3(4,6,7)
            case 209:  BIT_SET_ODD_4(4,6,7)
            case 210:  BIT_SET_EVEN_4(1,4,6,7)
            case 211:  BIT_SET_ODD_5(1,4,6,7)
            case 212:  BIT_SET_EVEN_4(2,4,6,7)
            case 213:  BIT_SET_ODD_5(2,4,6,7)
            case 214:  BIT_SET_EVEN_5(1,2,4,6,7)
            case 215:  BIT_SET_ODD_6(1,2,4,6,7)
            case 216:  BIT_SET_EVEN_4(3,4,6,7)
            case 217:  BIT_SET_ODD_5(3,4,6,7)
            case 218:  BIT_SET_EVEN_5(1,3,4,6,7)
            case 219:  BIT_SET_ODD_6(1,3,4,6,7)
            case 220:  BIT_SET_EVEN_5(2,3,4,6,7)
            case 221:  BIT_SET_ODD_6(2,3,4,6,7)
            case 222:  BIT_SET_EVEN_6(1,2,3,4,6,7)
            case 223:  BIT_SET_ODD_7(1,2,3,4,6,7)
            case 224:  BIT_SET_EVEN_3(5,6,7)
            case 225:  BIT_SET_ODD_4(5,6,7)
            case 226:  BIT_SET_EVEN_4(1,5,6,7)
            case 227:  BIT_SET_ODD_5(1,5,6,7)
            case 228:  BIT_SET_EVEN_4(2,5,6,7)
            case 229:  BIT_SET_ODD_5(2,5,6,7)
            case 230:  BIT_SET_EVEN_5(1,2,5,6,7)
            case 231:  BIT_SET_ODD_6(1,2,5,6,7)
            case 232:  BIT_SET_EVEN_4(3,5,6,7)
            case 233:  BIT_SET_ODD_5(3,5,6,7)
            case 234:  BIT_SET_EVEN_5(1,3,5,6,7)
            case 235:  BIT_SET_ODD_6(1,3,5,6,7)
            case 236:  BIT_SET_EVEN_5(2,3,5,6,7)
            case 237:  BIT_SET_ODD_6(2,3,5,6,7)
            case 238:  BIT_SET_EVEN_6(1,2,3,5,6,7)
            case 239:  BIT_SET_ODD_7(1,2,3,5,6,7)
            case 240:  BIT_SET_EVEN_4(4,5,6,7)
            case 241:  BIT_SET_ODD_5(4,5,6,7)
            case 242:  BIT_SET_EVEN_5(1,4,5,6,7)
            case 243:  BIT_SET_ODD_6(1,4,5,6,7)
            case 244:  BIT_SET_EVEN_5(2,4,5,6,7)
            case 245:  BIT_SET_ODD_6(2,4,5,6,7)
            case 246:  BIT_SET_EVEN_6(1,2,4,5,6,7)
            case 247:  BIT_SET_ODD_7(1,2,4,5,6,7)
            case 248:  BIT_SET_EVEN_5(3,4,5,6,7)
            case 249:  BIT_SET_ODD_6(3,4,5,6,7)
            case 250:  BIT_SET_EVEN_6(1,3,4,5,6,7)
            case 251:  BIT_SET_ODD_7(1,3,4,5,6,7)
            case 252:  BIT_SET_EVEN_6(2,3,4,5,6,7)
            case 253:  BIT_SET_ODD_7(2,3,4,5,6,7)
            case 254:  BIT_SET_EVEN_7(1,2,3,4,5,6,7)

            case 255:
                for( unsigned int i = 0; i < 8; ++i) {
                    if( m_cfier( *first ) ) {
                        res |= (m << i);
                    }
                    ++first;
                }
                break;
            }
            m <<= 8;
            b >>= 8;
        }

        return res;
    }
};

}   // namespace recombinations {
}   // namespace clotho {

#endif  // BLOCK_RECOMBINER_HPP_
