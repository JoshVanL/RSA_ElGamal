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
    str[strcspn(str, "\n\r")] = 0;

    for (int i=(WORD_LENGTH - 1); i >= 0; i--) {
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
        unsigned long n = 100;

        if (getline(&line, &n, stdin) == -1) {
            if (i != 0) {
                fprintf( stderr, "ERROR: epected %d further fields\n", size - i);
                exit(1);
            }
            return -1;
        }

        strToInt(field[i], line);
    }

    return 0;
}

// N, e, m
void RSAEncrypt(mpz_t N, mpz_t e, mpz_t message, mpz_t cipher) {
    mpz_init(cipher);
    mpz_powm(cipher, message, e, N);
}

//TODO: Put all input into this function
// N, d, p, q, d_p, d_q, i_p, i_q, c
//void RSAencrypt(mpz_t N, mpz_t d, mpz_t p, mpz_t d_p, mpz_t d_q, mpz_t i_p, mpz_t i_q, mpz_t c) {
void RSADecrypt(mpz_t N, mpz_t d, mpz_t cipher, mpz_t message) {
    mpz_init(message);
    mpz_powm(message, cipher, d, N);
}

// g, p q, h, message, c1, c2
void ElGamalEncrypt(mpz_t p, mpz_t q, mpz_t g, mpz_t h, mpz_t message, mpz_t c1, mpz_t c2) {
    // Make random k
    int k = 1;
    mpz_init(c1);
    mpz_init(c2);

    mpz_powm_ui(c1, g, k, p);
    mpz_pow_ui(c2, h, k);
    mpz_mul(c2, c2, message);
    mpz_mod(c2, c2, p);
}

// p, q, g, x, c1, c2, message
// * - read each 5-tuple of p, q, g, x and c = (c_1,c_2) from stdin,
void ElGamalDecryption(mpz_t p, mpz_t q, mpz_t g, mpz_t x, mpz_t c1, mpz_t c2, mpz_t message) {
    mpz_init(message);

    mpz_powm_sec(c1, c1, x, p);
    mpz_invert(c1, c1, p);
    mpz_mul(message, c1, c2);
    mpz_mod(message, message, p);
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
        RSAEncrypt(fields[0], fields[1], fields[2], cipher);

        char* out = intToStr(cipher);

        fprintf( stdout, "%s\n", out);
    }
}

/* Perform stage 2:
 *
 * - read each 9-tuple of N, d, p, q, d_p, d_q, i_p, i_q and c from stdin,
 * - compute the RSA decryption m, then
 * - write the plaintext m to stdout.
 */

void stage2() {
    const int exp_size = 9;

    mpz_t fields[exp_size];
    for (int i=0; i < exp_size; i++) {
        mpz_init(fields[i]);
    }

    while(readGroup(exp_size, fields) != -1) {
        mpz_t message;

        RSADecrypt(fields[0], fields[1], fields[8], message);

        char* out = intToStr(message);

        fprintf( stdout, "%s\n", out);
    }
}

/* Perform stage 3:
 *
 * - read each 5-tuple of p, q, g, h and m from stdin,
 * - compute the ElGamal encryption c = (c_1,c_2), then
 * - write the ciphertext c to stdout.
 */

void stage3() {
    const int exp_size = 5;

    mpz_t fields[exp_size];
    for (int i=0; i < exp_size; i++) {
        mpz_init(fields[i]);
    }

    while(readGroup(exp_size, fields) != -1) {
        mpz_t c1;
        mpz_t c2;
        ElGamalEncrypt(fields[0], fields[1], fields[2], fields[3], fields[4], c1, c2);

        char* out1 = intToStr(c1);
        char* out2 = intToStr(c2);
        fprintf( stdout, "%s\n", out1);
        fprintf( stdout, "%s\n", out2);
    }
}

/* Perform stage 4:
 * 
 * - read each 5-tuple of p, q, g, x and c = (c_1,c_2) from stdin,
 * - compute the ElGamal decryption m, then
 * - write the plaintext m to stdout.
 */
void stage4() {
    const int exp_size = 6;

    mpz_t fields[exp_size];
    for (int i=0; i < exp_size; i++) {
        mpz_init(fields[i]);
    }

    while(readGroup(exp_size, fields) != -1) {
        mpz_t message;
        ElGamalDecryption(fields[0], fields[1], fields[2], fields[3], fields[4], fields[5], message);

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
