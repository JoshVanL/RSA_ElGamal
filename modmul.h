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
