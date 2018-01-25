/* Copyright (C) 2017 Daniel Page <csdsp@bristol.ac.uk>
 *
 * Use of this source code is restricted per the CC BY-NC-ND license, a copy of
 * which can be found via http://creativecommons.org (and should be included as
 * LICENSE.txt within the associated archive or repository).
 */

#include "modmul.h"

#define BASE 16
#define WORD_LENGTH 256

/* Perform stage 1:
 *
 * - read each 3-tuple of N, e and m from stdin,
 * - compute the RSA encryption c, then
 * - write the ciphertext c to stdout.
 */

typedef struct {
    mpz_t N;
    mpz_t e;
    mpz_t m;
} RSA_public_key;

typedef struct {
    mpz_t N;
    mpz_t d;
    mpz_t p, q;
    mpz_t d_p, d_q;
    mpz_t i_p, i_q;
    mpz_t c;
} RSA_private_key;

typedef struct {
    mpz_t p;
    mpz_t q;
    mpz_t g;
    mpz_t h;
    mpz_t m;
} ElGamal_public_key;

typedef struct {
    mpz_t p;
    mpz_t q;
    mpz_t g;
    mpz_t x;
    mpz_t c1, c2;
} ElGamal_private_key;

int hexToInt(char ch) {
    if (ch >= '0' && ch <= '9')
        return ch - '0';
    if (ch >= 'A' && ch <= 'F')
        return ch - 'A' + 10;
    if (ch >= 'a' && ch <= 'f')
        return ch - 'a' + 10;
    return -1;
}

char intToHex(int n) {
    if (n > 9) {
        return 'A' + (n - 10);
    }
    return n + '0';
}


char* intToStr(mpz_t num) {
    int loc =0;
    char* str = malloc((WORD_LENGTH)*sizeof(char));
    mpz_t pow;
    mpz_t quot;
    mpz_init(pow);
    mpz_init(quot);

    for(int i=(WORD_LENGTH - 1); i>= 0; i--) {
        mpz_ui_pow_ui(pow, BASE, i);

        if (mpz_cmp(num, pow) < 0) {
            str[loc] = '0';

        } else {
            mpz_tdiv_q(quot, num, pow);
            mpz_mod(num, num, pow);

            int lquote = mpz_get_ui(quot);

            char c = intToHex(lquote);
            str[loc] = c;
        }
        loc++;
    }

    char * out = str;
    while (*out=='0') out++;

    return out;
}

void strToInt(mpz_t num, char* str) {
    int pow = 0;
    mpz_t tmp;
    mpz_init(num);
    mpz_init(tmp);
    int n = strcspn(str, "\n\r");
    str[n] = 0;

    for (int i=(n - 1); i >= 0; i--) {
        mpz_ui_pow_ui(tmp, BASE, pow);
        mpz_mul_si(tmp, tmp, hexToInt(str[i]));
        mpz_add(num, num, tmp);
        mpz_init(tmp);
        pow++;
    }
}

int readGroup(int size, mpz_t field[size]) {

    for (int i=0; i < size; i++) {
        char* line = NULL;
        size_t n = 1024;

        if (getline(&line, &n, stdin) == -1) {
            if (i != 0) {
                fprintf( stderr, "ERROR: expected %d further fields (total %d)\n", size - i, size);
                exit(1);
            }
            return -1;
        }

        strToInt(field[i], line);
    }

    return 0;
}

// N, e, m
void RSAEncrypt(RSA_public_key *pk, mpz_t cipher) {
    mpz_init(cipher);
    mpz_powm(cipher, pk->m, pk->e, pk->N);
}

//Uses CRT with Garner's formula
void RSADecrypt(RSA_private_key *sk, mpz_t message) {
    mpz_t m1, m2;
    mpz_init(m1);
    mpz_init(m2);
    mpz_init(message);

    mpz_powm(m1, sk->c, sk->d_p, sk->p);
    mpz_powm(m2, sk->c, sk->d_q, sk->q);
    mpz_sub(m1, m1, m2);
    mpz_mul(m1, sk->i_q, m1);
    mpz_mod(m1, m1, sk->p);
    mpz_mul(m1, m1, sk->q);
    mpz_add(message, m1, m2);
}

void ElGamalEncrypt(ElGamal_public_key *pk, mpz_t c1, mpz_t c2) {
    // Make random k
    int k = 1;
    mpz_init(c1);
    mpz_init(c2);

    mpz_powm_ui(c1, pk->g, k, pk->p);
    mpz_pow_ui(c2, pk->h, k);
    mpz_mul(c2, c2, pk->m);
    mpz_mod(c2, c2, pk->p);
}

void ElGamalDecryption(ElGamal_private_key *sk, mpz_t message) {
    mpz_init(message);

    mpz_powm(sk->c1, sk->c1, sk->x, sk->p);
    mpz_invert(sk->c1, sk->c1, sk->p);
    mpz_mul(message, sk->c2, sk->c1);
    mpz_mod(message, message, sk->p);
}

// m^e (mod N)
void stage1() {
    const int exp_size = 3;

    mpz_t fields[exp_size];
    for (int i=0; i < exp_size; i++) {
        mpz_init(fields[i]);
    }

    while(readGroup(exp_size, fields) != -1) {
        mpz_t cipher;
        RSA_public_key pk;
        mpz_init(pk.N);
        mpz_init(pk.e);
        mpz_init(pk.m);

        mpz_set(pk.N, fields[0]);
        mpz_set(pk.e, fields[1]);
        mpz_set(pk.m, fields[2]);

        RSAEncrypt(&pk, cipher);

        char* out = intToStr(cipher);
        fprintf( stdout, "%s\n", out);
    }
}

void stage2() {
    const int exp_size = 9;

    mpz_t fields[exp_size];
    for (int i=0; i < exp_size; i++) {
        mpz_init(fields[i]);
    }

    while(readGroup(exp_size, fields) != -1) {
        mpz_t message;
        RSA_private_key sk;
        mpz_init(sk.N);
        mpz_init(sk.d);
        mpz_init(sk.p);
        mpz_init(sk.q);
        mpz_init(sk.d_p);
        mpz_init(sk.d_q);
        mpz_init(sk.i_p);
        mpz_init(sk.i_q);
        mpz_init(sk.c);

        mpz_set(sk.N, fields[0]);
        mpz_set(sk.d, fields[1]);
        mpz_set(sk.p, fields[2]);
        mpz_set(sk.q, fields[3]);
        mpz_set(sk.d_p, fields[4]);
        mpz_set(sk.d_q, fields[5]);
        mpz_set(sk.i_p, fields[6]);
        mpz_set(sk.i_q, fields[7]);
        mpz_set(sk.c, fields[8]);

        RSADecrypt(&sk, message);

        char* out = intToStr(message);

        fprintf( stdout, "%s\n", out);
    }
}

void stage3() {
    const int exp_size = 5;

    mpz_t fields[exp_size];
    for (int i=0; i < exp_size; i++) {
        mpz_init(fields[i]);
    }

    while(readGroup(exp_size, fields) != -1) {
        mpz_t c1;
        mpz_t c2;
        ElGamal_public_key pk;
        mpz_init(pk.p);
        mpz_init(pk.q);
        mpz_init(pk.g);
        mpz_init(pk.h);
        mpz_init(pk.m);

        mpz_set(pk.p, fields[0]);
        mpz_set(pk.q, fields[1]);
        mpz_set(pk.g, fields[2]);
        mpz_set(pk.h, fields[3]);
        mpz_set(pk.m, fields[4]);

        ElGamalEncrypt(&pk, c1, c2);

        char* out1 = intToStr(c1);
        char* out2 = intToStr(c2);
        fprintf( stdout, "%s\n", out1);
        fprintf( stdout, "%s\n", out2);
    }
}

void stage4() {
    const int exp_size = 6;

    mpz_t fields[exp_size];
    for (int i=0; i < exp_size; i++) {
        mpz_init(fields[i]);
    }

    while(readGroup(exp_size, fields) != -1) {
        mpz_t message;
        ElGamal_private_key sk;
        mpz_init(sk.p);
        mpz_init(sk.q);
        mpz_init(sk.g);
        mpz_init(sk.x);
        mpz_init(sk.c1);
        mpz_init(sk.c2);

        mpz_set(sk.p, fields[0]);
        mpz_set(sk.q, fields[1]);
        mpz_set(sk.g, fields[2]);
        mpz_set(sk.x, fields[3]);
        mpz_set(sk.c1, fields[4]);
        mpz_set(sk.c2, fields[5]);

        ElGamalDecryption(&sk, message);

        char* out = intToStr(message);
        fprintf( stdout, "%s\n", out);
    }
}

/* The main function acts as a driver for the assignment by simply invoking the
 * correct function for the requested stage.
 */

int main( int argc, char* argv[] ) {
  if( 2 != argc ) {
    abort();
  }

  if     ( !strcmp( argv[ 1 ], "stage1" ) ) {
    stage1();
  }
  else if( !strcmp( argv[ 1 ], "stage2" ) ) {
    stage2();
  }
  else if( !strcmp( argv[ 1 ], "stage3" ) ) {
    stage3();
  }
  else if( !strcmp( argv[ 1 ], "stage4" ) ) {
    stage4();
  }
  else {
    abort();
  }

  return 0;
}
