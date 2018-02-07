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
    n->d = malloc(sizeof(limb_t)*32);
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

void myn_sub_ui(my_num_t* r, my_num_t* x, const unsigned long int y) {
    my_num_t *op2 = malloc(sizeof(my_num_t));
    myn_init(op2);
    myn_set_ui(op2, y);

    myn_sub(r, x, op2);
}


void myn_mul(my_num_t* r, my_num_t* x, my_num_t* y) {
    my_num_t *ans = malloc(sizeof(my_num_t));
    my_num_t *times = malloc(sizeof(my_num_t));
    my_num_t *one = malloc(sizeof(my_num_t));
    myn_init(ans);
    myn_init(times);
    myn_init(one);
    myn_set(times, y);
    myn_set_ui(one, 1);

    int i =0;

    while(times->l > 1 || times->d[0] > 0) {
        myn_add(ans, ans, x);
        myn_sub(times, times, one);
        i++;
    }

    (*r) = (*ans);
    //for (int i=31; i>=0; i--) {
    //    fprintf(stderr, ">%d\n", r->d[i]);
    //}
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

void myn_mod(my_num_t* x, my_num_t* mod) {
    while ( myn_cmp(x, mod) >= 0 ) {
        myn_sub(x, x, mod);
    }
}

int myn_mod_quote(my_num_t* x, my_num_t* mod) {
    int i = 0;

    while ( myn_cmp(x, mod) >= 0 ) {
        myn_sub(x, x, mod);
        i++;
    }

    return i;
}

int myn_digit_at(my_num_t* y, int n) {
    int k = n / 10;
    int p = n % 10;
    //y->d[k];
    //pth pos in k limb

    int r, t1, t2;
    t1 = pow(10, p+1);
    r = y->d[k] % t1;

    if (n > 0) {
        t2 = pow(10, p);
        r = r /t2;
    }

    return r;
}

void myn_xmuly(my_num_t* r, my_num_t* x, my_num_t* y) {
    myn_init(r);

    int limit = floor (log10 (abs (y->d[y->l-1]))) + 1;
    if (y->l > 2) {
        limit += (y->l -1) * 10;
    }

    for(int i=limit; i>=0; i--) {
        my_num_t tmp1;
        my_num_t tmp2;
        myn_init(&tmp1);
        myn_init(&tmp2);
        myn_mul_ui(&tmp1, r, 10);
        myn_mul_ui(&tmp2, x, myn_digit_at(y, i));
        myn_add(r, &tmp1, &tmp2);
    }
}

//TODO: Fix this
void myn_div10(my_num_t* r, my_num_t* x) {
    myn_init(r);
    int64_t invDivisor = 0x1999999A;
    for (int i=0; i<x->l; i++) {
        r->d[i] = (int32_t) ((invDivisor * x->d[i]) >> 32);
    }
}

void myn_xmulydiv10(my_num_t* r, my_num_t* x, my_num_t* y) {
    myn_init(r);

    int limit = floor (log10 (abs (y->d[y->l-1]))) + 1;
    if (y->l > 2) {
        limit += (y->l -1) * 10;
    }

    for (int i=0; i<limit; i++) {
        my_num_t tmp1;
        my_num_t tmp2;
        myn_init(&tmp1);
        myn_init(&tmp2);
        myn_mul_ui(&tmp1, x, myn_digit_at(y, i));
        myn_add(&tmp2, r, &tmp1);
        myn_div10(r, &tmp2);
    }
}

void myn_xmulymodN(my_num_t* r, my_num_t* x, my_num_t* y, my_num_t* N) {
    myn_init(r);

    int limit = floor (log10 (abs (y->d[y->l-1]))) + 1;
    if (y->l > 2) {
        limit += (y->l -1) * 10;
    }

    for (int i=0; i<limit; i++) {

    }

}

//void myn_mul_mod(my_num_t* r, my_num_t* x, my_num_t* y, my_num_t* N) {
//}

void print_myn(my_num_t *x) {
    for (int i=x->l-1; i>=0; i--) {
        fprintf(stderr, ">%u\n", x->d[i]);
    }
}

void myn_powm(my_num_t* r, my_num_t* base, my_num_t* expoent, my_num_t* modulo) {
    myn_init(r);
    myn_set_ui(r, 1);
    myn_mod(base, modulo);

    int i=0;

    while(expoent->l > 1 || expoent->d[0] > 0) {
        myn_mul(r, r, base);
        myn_mod(r, modulo);
        myn_sub_ui(expoent, expoent, 1);
        i++;
    }
}

void results(my_num_t *n, char* name) {
    fprintf(stderr, "Results: %s->l: %d, %s->d[1]: %u, %s->d[0]: %u\n", name, n->l, name, n->d[1], name, n->d[0]);
}

int test() {
    int fail = 0;

    my_num_t *x = malloc(sizeof(my_num_t));
    my_num_t *y = malloc(sizeof(my_num_t));
    myn_init(x);
    myn_init(y);

    myn_set_ui(x, 10);
    if (x->d[0] != 10) {
        fail += 1;
    }

   myn_add(x, x, y);
   if (x->d[0] != 10 || y->d[0] != 0) {
       results(x, "x");
       results(y, "y");
       fail += 1;
   }

   myn_set_ui(y, 10);
   if (y->d[0] != 10) {
       results(y, "y");
       fail += 1;
   }

   myn_add(x, x, y);
   if (x->d[0] != 20 || y->d[0] != 10) {
       results(x, "x");
       results(y, "y");
       fail += 1;
   }
   for (int i=x->l-1; i>0; i--) {
       if (x->d[i] != 0) {
           results(x, "x");
           fail += 1;
       }
   }

   myn_set_ui(x, 20);
   myn_mul_ui(x, x, 10);
   if (x->d[0] != 200) {
       results(x, "x");
       fail += 1;
   }
   for (int i=31; i>0; i--) {
       if (x->d[i] != 0) {
           results(x, "x");
           fail += 1;
       }
   }

   x->d[1] = 200;
   x->l = 2;
   myn_mul_ui(x, x, 10);
   if (x->d[0] != 2000 || x->d[1] != 2000) {
       results(x, "x");
       fail += 1;
   }
   for (int i=31; i>1; i--) {
       if (x->d[i] != 0) {
           results(x, "x");
           fail += 1;
       }
   }

   myn_ui_pow_ui(x, 2, 2);
   if (x->d[0] != 4) {
       results(x, "x");
       fail += 1;
   }
   for (int i=31; i>0; i--) {
       if (x->d[i] != 0) {
           results(x, "x");
           fail += 1;
       }
   }

   myn_ui_pow_ui(x, 2, 31);
   if (x->l != 1 ) {
       results(x, "x");
       fail += 1;
   }
   myn_ui_pow_ui(x, 2, 32);
   if (x->l != 2 || x->d[0] != 0 || x->d[1] != 1) {
       results(x, "x");
       fail += 1;
   }

   myn_init(x);
   myn_init(y);
   myn_set_ui(x, 10);
   myn_set_ui(y, 7);
   int quote = myn_mod_quote(x, y);
   if (x->l != 1 || x->d[0] != 3 || quote != 1 ) {
       results(x, "x");
       fail += 1;
   }

   myn_init(x);
   myn_init(y);
   myn_ui_pow_ui(x, 2, 32);

   myn_set_ui(y, 1);

   myn_add(x, x, y);
   if (x->l != 2 || x->d[0] != 1 || x->d[1] != 1 ) {
       results(x, "x");
       results(y, "y");
       fail += 1;
   }
   myn_ui_pow_ui(y, 2, 32);
   quote = myn_mod_quote(x, y);
   if (x->l != 1 || x->d[0] != 1 || quote != 1 ) {
       results(x, "x");
       fail += 1;
   }

   myn_init(x);
   myn_init(y);
   myn_set_ui(x, 5);
   myn_set_ui(y, 4);
   if( myn_cmp(x, y) != 1) {
       fail += 1;
   }
   myn_set_ui(x, 3);
   if( myn_cmp(x, y) != -1) {
       fail += 1;
   }
   myn_set_ui(x, 4);
   if( myn_cmp(x, y) != 0) {
       fail += 1;
   }

   myn_init(x);
   myn_init(y);
   myn_ui_pow_ui(x, 2, 32);
   myn_set_ui(y, 100);
   myn_sub(x, x, y);
   if(x->l != 1 || x->d[0] != (UINT32_MAX-100) || y->l != 1 || y->d[0] != 100) {
       results(x, "x");
       results(y, "y");
       fail += 1;
   }

   myn_init(x);
   myn_init(y);
   x->l = 3;
   x->d[2] = 1;
   myn_set_ui(y, 1);
   myn_sub(x, x, y);
   if (x->l != 2) fail += 1;
   if (x->d[2] != 0) fail += 1;
   for (int i=1; i>=0; i--) {
       if (x->d[i] != UINT32_MAX-1) {
           fail += 1;
           break;
       }
   }

   my_num_t *z = malloc(sizeof(my_num_t));
   my_num_t *a = malloc(sizeof(my_num_t));
   myn_init(x);
   myn_init(y);
   myn_init(z);
   myn_init(a);
   myn_set_ui(x, 3);
   myn_set_ui(y, 3);
   myn_set_ui(z, 7);
   //myn_mul(x, z, y);
   myn_powm(a, x, y, z);
    if (a->d[0] != 6) {
        results(a, "a");
        fail++;
    }

    myn_init(x);
    x->l = 2;
    x->d[1] = 200;
    if (myn_digit_at(x, 12) != 2) {
        results(x, "x");
        fail++;
    }

    myn_init(x);
    myn_init(y);
    myn_init(z);
    myn_set_ui(x, 123);
    myn_set_ui(y, 456);
    myn_xmuly(z, x, y);
    if (z->d[0] != 56088) {
        results(z, "z");
        fail++;
    }

    myn_init(x);
    myn_init(y);
    myn_set_ui(x, 100);
    myn_div10(y, x);
    if (y->d[0] != 10) {
        results(y, "y");
        fail++;
    }

    myn_init(x);
    myn_init(y);
    myn_init(z);
    myn_set_ui(x, 123000);
    myn_set_ui(y, 456000);
    myn_xmulydiv10(z, x, y);
    if (z->d[0] != 56088) {
        results(z, "z");
        fail++;
    }


    return fail;
}
