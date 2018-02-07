/* Copyright (C) 2017 Daniel Page <csdsp@bristol.ac.uk>
 *
 * Use of this source code is restricted per the CC BY-NC-ND license, a copy of
 * which can be found via http://creativecommons.org (and should be included as
 * LICENSE.txt within the associated archive or repository).
 */

#ifndef __MODMUL_H
#define __MODMUL_H

#include  <stdio.h>
#include <stdlib.h>

#include <string.h>
#include    <gmp.h>

#include <unistd.h>
#include <fcntl.h>
#include <math.h>
#include <stdint.h>

#endif

#define SIGN_POS 0
#define SIGN_NEG 1

typedef uint32_t limb_t;

typedef struct my_num_t {
    limb_t* d; //  Conatins the actual numbers

    int l; // The number of limbs in use
    int s; // The sign of the number
} my_num_t;

typedef struct {
    my_num_t N, e;
    my_num_t m;
} RSA_public_key;

typedef struct {
    mpz_t N, d;
    mpz_t p, q;
    mpz_t d_p, d_q;
    mpz_t i_p, i_q;
    mpz_t c;
} RSA_private_key;

typedef struct {
    mpz_t p, q;
    mpz_t g;
    mpz_t h;
    mpz_t k;
    mpz_t m;
} ElGamal_public_key;

typedef struct {
    mpz_t p, q;
    mpz_t g;
    mpz_t x;
    mpz_t c1, c2;
} ElGamal_private_key;

int min(int x, int y) {
    if (x > y) {
        return y;
    }

    return x;
}

int max(int x, int y) {
    if (x >= y) {
        return x;
    }

    return y;
}

#define LIMB_ADD0(r_1, r_0, e, g) { \
  asm( "addl %2 ,%0 ; adcl $0 ,%1 ;    \
        addl %3 ,%0 ; adcl $0 ,%1 ;"   \
                                       \
        : "+&r" (r_0), "+&r" (r_1)     \
        :   "r" (e),     "r" (g)       \
        : "cc");                       \
}

#define LIMB_ADD1(r_1, r_0, e, f, g) { \
  asm( "movl %2 ,%0 ; movl $0 ,%1 ;    \
        addl %3 ,%0 ; adcl $0 ,%1 ;    \
        addl %4 ,%0 ; adcl $0 ,%1 ;"   \
                                       \
        : "+&r" (r_0), "+&r" (r_1)     \
        :   "r" (e),     "r" (f),      \
            "r" (g)                    \
        : "cc");                       \
}

#define LIMB_SUB1(r_1, r_0, e, f, g) { \
  asm( "movl %2 ,%0 ; movl $0 ,%1 ;    \
        subl %3 ,%0 ; adcl $0 ,%1 ;    \
        subl %4 ,%0 ; adcl $0 ,%1 ;"   \
                                       \
        : "+&r" (r_0), "+&r" (r_1)     \
        :   "r" (e),     "r" (f),      \
            "r" (g)                    \
        : "cc");                       \
}

// -1 = x < y ;; 0 = x == y ;; 1 = x > y
int limbs_cmp( const limb_t* x, int l_x, const limb_t* y, int l_y ) {
    int i = max(l_x, l_y);

    if (i > l_x) {
        return -1;
    } else if ( i > l_y ) {
        return 1;
    }

    while( i >= 0 ) {
        if ( x[i] > y[i] ) {
            return 1;
        }
        if ( x[i] < y[i] ) {
            return -1;
        }
        i--;
    }

    return 0;
}

//TODO: Fix lIMB_ADD0 and LIMB_SUB0
limb_t limb_add( limb_t* r, const limb_t* x, int l_x , const limb_t* y, int l_y ) {
    int l_r = min( l_x , l_y ), i = 0;

    limb_t c = 0;

    while( i < l_r ) {
        limb_t d_x = x[ i ];
        limb_t d_y = y[ i ];

        LIMB_ADD1 ( c, r[ i ], d_x , d_y , c ); i++;
    }
    while( i < l_x ) {
        limb_t d_x = x[ i ];

        LIMB_ADD1 ( c, r[ i ], d_x, 0, c ); i++;
    }
    while( i < l_y ) {
        limb_t d_y = y[ i ];

        LIMB_ADD1 ( c, r[ i ], d_y, 0, c ); i++;
    }

    return c;
}

int limb_sub( limb_t* r, const limb_t* x, int l_x , const limb_t* y, int l_y ) {
    limb_t c = 0;
    int s = 0;
    int l_r = max(l_x, l_y);

    // Second oporand is bigger
    if ( limbs_cmp(x, l_x, y, l_y) < 0 ) {
        s = 1;
        for (int i=0; i<l_r; i++) {
            limb_t d_x = x[ i ];
            limb_t d_y = y[ i ];

            LIMB_SUB1 ( c, r[ i ], d_y , d_x , c );
        }

    } else {
        for (int i=0; i<l_r; i++) {
            limb_t d_x = x[ i ];
            limb_t d_y = y[ i ];

            LIMB_SUB1 ( c, r[ i ], d_x , d_y , c );

            if (c != 0) {
                int d = 0;
                LIMB_SUB1 ( d, r[ i ], r[ i ] , c , d );
            }
        }
    }

    return s;
}

#define MYN_ADD(r, x, y, sign) {                              \
    int l_r = max(x->l, y->l) + 1;                            \
    r->d[l_r - 1] = limb_add( r->d, x->d, x->l, y->d, y->l ); \
    r->l  = myn_lop(r->d, l_r);                               \
    r->s = sign;                                              \
}

#define MYN_SUB(r, x, y, sign) {                      \
    int l_r = max(x->l, y->l);                        \
    int s = limb_sub( r->d, x->d, x->l, y->d, y->l ); \
    r->l = myn_lop(r->d, l_r);                        \
    r->s = (sign + s) % 2;                            \
}

int myn_lop(const limb_t* x, int l_x ) {
    while( ( l_x > 1 ) && ( x[ l_x - 1 ] == 0 ) ) {
        l_x --;
    }

    return l_x;
}

void myn_init(my_num_t* const n) {
	n->d =  malloc(sizeof(limb_t)*32);
    for (int i=0; i<32; i++) {
        n->d[i] = 0;
    }
    n->l = 1;
    n->s = 0;
}

void myn_add( my_num_t* r, const my_num_t* x, const my_num_t* y ) {
    if ( ( x->s == SIGN_NEG ) && ( y->s == SIGN_POS ) ) {
        if( limbs_cmp( x->d, x->l, y->d, y->l ) >= 0 ) {
            MYN_SUB( r, x, y, SIGN_NEG ); // r = -( abs(x) - abs(y) )
        }
        else {
            MYN_SUB( r, y, x, SIGN_POS ); // r = +( abs(y) - abs(x) )
        }
    }
    else if( ( x->s == SIGN_POS ) && ( y->s == SIGN_NEG ) ) {
        if( limbs_cmp( x->d, x->l, y->d, y->l ) >= 0 ) {
            MYN_SUB( r, x, y, SIGN_POS ); // r = +( abs(x) - abs(y) )
        }
        else {
            MYN_SUB( r, y, x, SIGN_NEG ); // r = -( abs(y) - abs(x) )
        }
    }
    else if( ( x->s == SIGN_POS ) && ( y->s == SIGN_POS ) ) {
        MYN_ADD( r, x, y, SIGN_POS ); // r = +( abs(x) + abs(y) )
    }
    else if( ( x->s == SIGN_NEG ) && ( y->s == SIGN_NEG ) ) {
        MYN_ADD( r, x, y, SIGN_NEG ); // r = -( abs(x) + abs(y) )
    }
}

void myn_sub( my_num_t* r, const my_num_t* x, const my_num_t* y ) {
    if ( ( x->s == SIGN_POS ) && ( y->s == SIGN_NEG ) ) {
        MYN_ADD( r, x, y, SIGN_POS ); // r = +( abs(x) + abs(y) )
    }
    else if( ( x->s == SIGN_NEG ) && ( y->s == SIGN_POS ) ) {
        MYN_ADD( r, x, y, SIGN_NEG ); // r = -( abs(x) + abs(y) )
    }
    else if( ( x->s == SIGN_POS ) && ( y->s == SIGN_POS ) ) {
        MYN_SUB( r, x, y, SIGN_POS ); // r = +( abs(x) - abs(y) )
    }
    else if( ( x->s == SIGN_NEG ) && ( y->s == SIGN_NEG ) ) {
        MYN_ADD( r, x, y, SIGN_NEG ); // r = -( abs(x) - abs(y) )
    }
}

void myn_set_ui(my_num_t* x, const unsigned long int y) {
    x->d[0] = y;
    x->l = 1;
}

void myn_set(my_num_t* x, my_num_t* y) {
    myn_init(x);
    for (int i=0; i<32; i++) {
        x->d[i] = y->d[i];
    }

    x->l = y->l;
    x->s = y->s;
}

void myn_mul_ui(my_num_t* r, my_num_t* x, const unsigned long int y) {
    my_num_t *ans = malloc(sizeof(my_num_t));
    myn_init(ans);


    for (int i=0; i<y; i++) {
        myn_add(ans, ans, x);
    }

    (*r) = (*ans);
}

void myn_ui_pow_ui(my_num_t* r, const unsigned long int base, const unsigned long int expoent) {
    myn_init(r);

    if (expoent == 0) {
        myn_set_ui(r, 1);
        return;
    }

    myn_set_ui(r, base);
    for (int i=1; i<expoent; i++) {
        myn_mul_ui(r, r, base);
    }
}

int myn_cmp(my_num_t* x, my_num_t* y) {
    return limbs_cmp(x->d, x->l, y->d, y->l);
}

int myn_mod_quote(my_num_t* x, my_num_t* mod) {
    int i = 0;

    while ( myn_cmp(x, mod) >= 0 ) {
        myn_sub(x, x, mod);
        i++;
    }

    return i;
}
