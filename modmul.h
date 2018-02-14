
///////////////////////////////////////////////////////////
//                                                       //
//                 Joshua Van Leeuwen                    //
//                                                       //
//                University of Bristol                  //
//                                                       //
///////////////////////////////////////////////////////////

#ifndef __MODMUL_H
#define __MODMUL_H

#include  <stdio.h>
#include <stdlib.h>

#include <string.h>
#include    <gmp.h>

#include <unistd.h>
#include <fcntl.h>

#endif

typedef struct {
    mpz_t N, e;
    mpz_t m;
    mpz_t ro2;
    mp_limb_t *o;
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
    mpz_t g, h;
    mpz_t k;
    mpz_t m;
    mpz_t ro2;
    mp_limb_t *o;
} ElGamal_public_key;

typedef struct {
    mpz_t p, q;
    mpz_t g;
    mpz_t x;
    mpz_t c1, c2;
    mpz_t ro2;
    mp_limb_t *o;
} ElGamal_private_key;

typedef struct {
    mpz_t p, q;
    mpz_t N;
    mpz_t s;
    mpz_t ro2, two;
    mp_limb_t *o;
} BBS;
