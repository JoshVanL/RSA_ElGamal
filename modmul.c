/* Copyright (C) 2017 Daniel Page <csdsp@bristol.ac.uk>
 *
 * Use of this source code is restricted per the CC BY-NC-ND license, a copy of
 * which can be found via http://creativecommons.org (and should be included as
 * LICENSE.txt within the associated archive or repository).
 */

#include "modmul.h"

#define INPUT_BASE 16

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


void intToStr(char* str, mpz_t num) {
    int loc =0;
    mpz_t curr;
    mpz_t pow;
    mpz_t quot;
    mpz_inits(curr, pow, quot);

    str = malloc(256*sizeof(char));
    for(int i=255; i>= 0; i--) {
        mpz_ui_pow_ui(pow, 16, i);

        if (mpz_cmp(curr, pow) < 0) {
            str[loc] = '0';

        } else {
            mpz_tdiv_q(quot, curr, pow);
            mpz_mod(curr, curr, pow);
            sprintf(str, "%x", num);

        }
        str++;
    }
}

void strToInt(mpz_t num, char* str) {
    int pow = 0;
    mpz_t tmp;
    mpz_inits(num, tmp);
    str[strcspn(str, "\n\r")] = 0;

    for (int i=strlen(str); i >= 0; i--) {
        mpz_ui_pow_ui(tmp, hexToInt(str[i]), pow);
        mpz_add(num, num, tmp);
        pow++;
    }
}

int readGroup(mpz_t *field, int size) {

    for (int i=0; i < size; i++) {
        char* line;
        size_t n = 1024;

        if (getline(&line, &n, stdin) == -1) {
            fprintf( stderr, "ERROR: failed to read line [%d]", i);
            return -1;
        }

        strToInt(field[i], line);
        //fprintf( stdout, "%s", line);
        //gmp_printf("%Zd\n", field[i]);
        //fprintf( stdout, "\n");
    }

    return 0;
}

// N, e, m
void RSAencrypt(mpz_t N, mpz_t e, mpz_t message, mpz_t cipher) {
    mpz_powm(cipher, message, e, N);
}

// m^e (mod N)
void stage1() {
    const int exp_size = 3;

    mpz_t fields[exp_size];
    for (int i=0; i < exp_size; i++) {
        mpz_init(fields[i]);
    }

    if (readGroup(fields, exp_size) == -1) {
        exit(1);
    }


//gmp_printf ("%s is an mpz %Zd\n", "here", z);

    //for (int i=0; i < exp_size; i++) {
    //    gmp_printf(">%Zd\n", fields[i]);
    //}

    mpz_t cipher;
    char* out = NULL;

    RSAencrypt(fields[0], fields[1], fields[2], cipher);
    intToStr(out, cipher);

    fprintf( stdout, "%s\n", out);
}

/* Perform stage 2:
 *
 * - read each 9-tuple of N, d, p, q, d_p, d_q, i_p, i_q and c from stdin,
 * - compute the RSA decryption m, then
 * - write the plaintext m to stdout.
 */

void stage2() {

  // fill in this function with solution

}

/* Perform stage 3:
 *
 * - read each 5-tuple of p, q, g, h and m from stdin,
 * - compute the ElGamal encryption c = (c_1,c_2), then
 * - write the ciphertext c to stdout.
 */

void stage3() {

  // fill in this function with solution

}

/* Perform stage 4:
 * 
 * - read each 5-tuple of p, q, g, x and c = (c_1,c_2) from stdin,
 * - compute the ElGamal decryption m, then
 * - write the plaintext m to stdout.
 */

void stage4() {

  // fill in this function with solution

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
