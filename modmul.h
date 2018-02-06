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

typedef struct {
    mpz_t N, e;
    mpz_t m;
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

typedef uint32_t limb_t;

typedef struct my_num {
    limb_t d[ 32 ]; //  Conatins the actual numbers

    int l; // The number of limbs in use
    int s; // The sign of the number
} my_num;

int min(int x, int y) {
    if (x > y) {
        return y;
    }

    return x;
}

int max(int x, int y) {
    if (x > y) {
        return x;
    }

    return y;
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

int myn_lop(const limb_t* x, int l_x ) {
    while( ( l_x > 1 ) && ( x[ l_x - 1 ] == 0 ) ) {
        l_x --;
    }

    return l_x;
}

void myn_add(my_num *r, my_num *x, my_num *y, int sign) {
    int l_r = max(x->l, y->l) + 1;

    r->d[l_r - 1] = limb_add( r->d, x->d, x->l, y->d, y->l);

    r->l  = myn_lop(r->d, l_r);
    r->s = sign;
}

void myn_init(my_num **n) {
	*n = (my_num*) malloc(sizeof(my_num));
    for (int i=0; i<32; i++) {
        (*n)->d[i] = 0;
    }
    (*n)->l = 1;
    (*n)->s = 0;
}
