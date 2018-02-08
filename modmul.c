/* Copyright (C) 2017 Daniel Page <csdsp@bristol.ac.uk>
 *
 * Use of this source code is restricted per the CC BY-NC-ND license, a copy of
 * which can be found via http://creativecommons.org (and should be included as
 * LICENSE.txt within the associated archive or repository).
 */

#include "modmul.h"

#define BASE 16
#define WORD_LENGTH 256
#define P_MOD_SIZE 1024

void init_RSA_pk(RSA_public_key* const pk) {
	//*pk = (RSA_public_key*) malloc(sizeof(RSA_public_key));
    myn_init(&(pk->N));
    myn_init(&(pk->e));
    myn_init(&(pk->m));
}

void init_RSA_sk(RSA_private_key **sk) {
	*sk = (RSA_private_key*) malloc(sizeof(RSA_private_key));
    mpz_init((*sk)->N);
    mpz_init((*sk)->d);
    mpz_init((*sk)->p);
    mpz_init((*sk)->q);
    mpz_init((*sk)->d_p);
    mpz_init((*sk)->d_q);
    mpz_init((*sk)->i_p);
    mpz_init((*sk)->i_q);
    mpz_init((*sk)->c);
}

void init_ElGamal_pk(ElGamal_public_key **pk) {
	*pk = (ElGamal_public_key*) malloc(sizeof(ElGamal_public_key));
    mpz_init((*pk)->p);
    mpz_init((*pk)->q);
    mpz_init((*pk)->g);
    mpz_init((*pk)->h);
    mpz_init((*pk)->m);
}

void init_ElGamal_sk(ElGamal_private_key **sk) {
	*sk = (ElGamal_private_key*) malloc(sizeof(ElGamal_private_key));
    mpz_init((*sk)->p);
    mpz_init((*sk)->q);
    mpz_init((*sk)->g);
    mpz_init((*sk)->x);
    mpz_init((*sk)->c1);
    mpz_init((*sk)->c2);
}

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


char* zToStr(my_num_t* num) {
    int loc =0;
    char* str = malloc((WORD_LENGTH)*sizeof(char));
    my_num_t *power = malloc(sizeof(my_num_t));

    for(int i=(WORD_LENGTH - 1); i>= 0; i--) {
        myn_init(power);
        myn_ui_pow_ui(power, BASE, i);

        if (myn_cmp(num, power) < 0) {
            str[loc] = '0';

        } else {
            int lquote = myn_mod_quote(num, power);

            char c = intToHex(lquote);
            str[loc] = c;
        }
        loc++;
    }

    char * out = str;
    while (*out=='0') out++;
    return out;
}

void strToZ(my_num_t* num, char* str) {
    int n = strcspn(str, "\n\r");
    str[n] = 0;

    myn_init(num);
    my_num_t *tmp = malloc (sizeof(my_num_t));
    myn_init(tmp);
    int power = 0;

    for (int i=(n - 1); i >= 0; i--) {
        myn_ui_pow_ui(tmp, BASE, power);
        myn_mul_ui(tmp, tmp, hexToInt(str[i]));
        myn_add(num, num, tmp);
        power++;
    }
}

int readGroup(int size, my_num_t field[size]) {
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

        strToZ(&field[i], line);
        //for (int j=field[i].l-1; j>=0; j--) {
        //    fprintf(stderr, ">%d\n", field[i].d[j]);
        //}
    }

    return 0;
}

int genRandomKey(mpz_t k, const size_t size) {
	char* dat = (char*) malloc(sizeof(char)*size);

    int file = open("/dev/random", O_RDONLY);
    if (file < 0) {
        fprintf(stderr, "ERROR: something went wrong opening /dev/random");
        return -1;
    }

	if (read(file, dat, size) < 0) {
        fprintf(stderr, "ERROR: something went wrong reading /dev/random");
        return -1;
    }

    if (close(file) < 0) {
        fprintf(stderr, "ERROR: something went wrong closing /dev/random");
        return -1;
    }

	mpz_import(k, size, 1, sizeof(char), 0, 0, dat);

    return 0;
}

// N, e, m
void RSAEncrypt(RSA_public_key *pk, my_num_t *cipher) {
    myn_init(cipher);
    my_num_t* rr = malloc(sizeof(my_num_t));
    uint32_t o;

    myn_mont_init(rr, &o, &(pk->N));
    myn_mont_mul(&pk->m, &pk->m, rr, &pk->N, o);
    myn_exp_mod_mont(cipher, &pk->m, &pk->e, &pk->N, rr, o);
    myn_mont_red(cipher, &pk->N, o);

    //myn_powm(cipher, &pk->m, &pk->e, &pk->N);
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
    mpz_init(pk->k);
    //Enable for random generation
    //if (genRandomKey(pk->k, P_MOD_SIZE) < 0) {
    //    return;
    //}
    //mpz_mod(pk->k, pk->k, pk->p);

    //Enable for testing
    mpz_set_ui(pk->k, 1);

    mpz_init(c1);
    mpz_init(c2);

    mpz_powm_sec(c1, pk->g, pk->k, pk->p);
    mpz_powm_sec(c2, pk->h, pk->k, pk->p);
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

    my_num_t fields[exp_size];
    for (int i=0; i < exp_size; i++) {
        myn_init(&fields[i]);
    }

    while(readGroup(exp_size, fields) != -1) {
        my_num_t cipher;
        RSA_public_key *pk = malloc(sizeof(RSA_public_key));
        init_RSA_pk(pk);

        myn_set(&pk->N, &fields[0]);
        myn_set(&pk->e, &fields[1]);
        myn_set(&pk->m, &fields[2]);

        RSAEncrypt(pk, &cipher);

        //char out[256];
        //for (int i=0; i<256; i++)
        //    out[i] = 0;



        char* out = zToStr(&cipher);
        fprintf( stdout, "%s\n", out);
    }
}

//void stage2() {
//    const int exp_size = 9;
//
//    mpz_t fields[exp_size];
//    for (int i=0; i < exp_size; i++) {
//        mpz_init(fields[i]);
//    }
//
//
//    while(readGroup(exp_size, fields) != -1) {
//        mpz_t message;
//        RSA_private_key *sk;
//        init_RSA_sk(&sk);
//
//        mpz_set(sk->N, fields[0]);
//        mpz_set(sk->d, fields[1]);
//        mpz_set(sk->p, fields[2]);
//        mpz_set(sk->q, fields[3]);
//        mpz_set(sk->d_p, fields[4]);
//        mpz_set(sk->d_q, fields[5]);
//        mpz_set(sk->i_p, fields[6]);
//        mpz_set(sk->i_q, fields[7]);
//        mpz_set(sk->c, fields[8]);
//
//        RSADecrypt(sk, message);
//
//        char* out = zToStr(message);
//
//        fprintf( stdout, "%s\n", out);
//    }
//}
//
//void stage3() {
//    const int exp_size = 5;
//
//    mpz_t fields[exp_size];
//    for (int i=0; i < exp_size; i++) {
//        mpz_init(fields[i]);
//    }
//
//    while(readGroup(exp_size, fields) != -1) {
//        mpz_t c1, c2;
//        ElGamal_public_key *pk;
//        init_ElGamal_pk(&pk);
//
//        mpz_set(pk->p, fields[0]);
//        mpz_set(pk->q, fields[1]);
//        mpz_set(pk->g, fields[2]);
//        mpz_set(pk->h, fields[3]);
//        mpz_set(pk->m, fields[4]);
//
//        ElGamalEncrypt(pk, c1, c2);
//
//        char* out1 = zToStr(c1);
//        char* out2 = zToStr(c2);
//        fprintf( stdout, "%s\n", out1);
//        fprintf( stdout, "%s\n", out2);
//    }
//}
//
//void stage4() {
//    const int exp_size = 6;
//
//    mpz_t fields[exp_size];
//    for (int i=0; i < exp_size; i++) {
//        mpz_init(fields[i]);
//    }
//
//    while(readGroup(exp_size, fields) != -1) {
//        mpz_t message;
//        ElGamal_private_key *sk;
//        init_ElGamal_sk(&sk);
//
//        mpz_set(sk->p, fields[0]);
//        mpz_set(sk->q, fields[1]);
//        mpz_set(sk->g, fields[2]);
//        mpz_set(sk->x, fields[3]);
//        mpz_set(sk->c1, fields[4]);
//        mpz_set(sk->c2, fields[5]);
//
//        ElGamalDecryption(sk, message);
//
//        char* out = zToStr(message);
//        fprintf( stdout, "%s\n", out);
//    }
//}

/* The main function acts as a driver for the assignment by simply invoking the
 * correct function for the requested stage.
 */

int main( int argc, char* argv[] ) {


    if ( !strcmp( argv[ 1 ], "test" ) ) {
        int fail = test();
        if (fail > 0) {
            fprintf(stderr, "[%d] TESTS FAILED.\n", fail);
            return 1;
        }
        fprintf(stdout, "ALL TESTS PASSED.\n");
        return 0;
    }

    if( 2 != argc ) {
        fprintf(stderr, "Expected 2 args; got %d\n", argc);
        return 1;
    }

    if ( !strcmp( argv[ 1 ], "stage1" ) ) {
        stage1();
    }

    //else if( !strcmp( argv[ 1 ], "stage2" ) ) {
    //  stage2();
    //}
    //else if( !strcmp( argv[ 1 ], "stage3" ) ) {
    //  stage3();
    //}
    //else if( !strcmp( argv[ 1 ], "stage4" ) ) {
    //  stage4();
    //}
    else {
        fprintf(stderr, "ERROR: Unknown stage '%s'\n", argv[1]);
        return 1;
    }

    return 0;
}
