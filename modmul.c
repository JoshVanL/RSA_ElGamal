// Joshua Van Leeuwen
// University of Bristol

#include "modmul.h"

#define BASE 16
#define WORD_LENGTH 256
#define P_MOD_SIZE 1024

void init_RSA_pk(RSA_public_key **pk) {
    *pk = (RSA_public_key*) malloc(sizeof(RSA_public_key));

    mpz_inits((*pk)->N, (*pk)->e,  NULL);
    mpz_init((*pk)->m);
}

void init_RSA_sk(RSA_private_key **sk) {
    *sk = (RSA_private_key*) malloc(sizeof(RSA_private_key));

    mpz_inits((*sk)->N,   (*sk)->d,   (*sk)->p,   (*sk)->q,   NULL);
    mpz_inits((*sk)->d_p, (*sk)->d_q, (*sk)->i_p, (*sk)->i_q, NULL);
    mpz_init( (*sk)->c);
}

void init_ElGamal_pk(ElGamal_public_key **pk) {
    *pk = (ElGamal_public_key*) malloc(sizeof(ElGamal_public_key));

    mpz_inits((*pk)->p, (*pk)->q, (*pk)->g, (*pk)->h, NULL);
    mpz_init( (*pk)->m);
}

void init_ElGamal_sk(ElGamal_private_key **sk) {
    *sk = (ElGamal_private_key*) malloc(sizeof(ElGamal_private_key));

    mpz_inits((*sk)->p,  (*sk)->q,  (*sk)->g, (*sk)->x, NULL);
    mpz_inits((*sk)->c1, (*sk)->c2, NULL);
}


void montgomery_CRT_init(mpz_t rr, mpz_t rrr, mpz_t N) {
    mpz_set(rrr, rr);

    for (int i = 1; i <= N->_mp_size << 6; i++) {
        int overflow = mpn_lshift(rrr->_mp_d, rrr->_mp_d, rrr->_mp_size, 1);
        if (overflow > 0) {
            rrr->_mp_d[rrr->_mp_size] = overflow;
            rrr->_mp_size++;
        }

        if (mpz_cmp(rrr, N) >= 0) {
            mpz_sub(rrr, rrr, N);
        }
    }
}


int getBit(mpz_t x, int n) {
    if (x->_mp_size > n >> 6) {
        return (x->_mp_d[n >> 6] >> (n & 63)) & 1;
    }

    return 0;
}

void montgomery_init(mpz_t rr, uint64_t *o, mpz_t N) {
    *o = 1;
    for (int i = 1; i < 64; i++) {
        *o *= *o;
        *o *= N->_mp_d[0];
    }
    *o *= -1;

    mpz_set_ui(rr, 1);
    _mpz_realloc(rr, N->_mp_size); //FIX THIS

    for (int i = 1; i <= N->_mp_size << 7; i++) {
        int result = mpn_lshift(rr->_mp_d, rr->_mp_d, rr->_mp_size, 1);
        if (result) {
            rr->_mp_d[rr->_mp_size] = result;
            rr->_mp_size++;
        }

        if (mpz_cmp(rr, N) >= 0) {
            mpz_sub(rr, rr, N);
        }
    }
}


void montgomery_mul(mpz_t r, mpz_t x, mpz_t y, mpz_t N, uint64_t o) {
    mpz_t tmp, t, yi;
    mpz_init_set_ui(tmp, 0);
    mpz_inits(yi, t, NULL);

    for (int i = 0; i < N->_mp_size; i++) {
        ////////////////////////////////////////////////// FIX
        uint64_t u = (x->_mp_d[0] * ((y->_mp_size > i) ? y->_mp_d[i] : 0) + tmp->_mp_d[0]) * o;

        mpz_mul_ui(t, N, u);
        mpz_mul_ui(yi, x, (y->_mp_size > i) ? y->_mp_d[i] : 0);

        mpz_add(tmp, tmp, yi);
        mpz_add(tmp, tmp, t);

        mpn_rshift(tmp->_mp_d, tmp->_mp_d, tmp->_mp_size, 32);
        mpn_rshift(tmp->_mp_d, tmp->_mp_d, tmp->_mp_size, 32);

        tmp->_mp_size -= 1;
    }

    if (mpz_cmp(tmp, N) >= 0) {
        mpz_sub(tmp, tmp, N);
    }

    mpz_set(r, tmp);
    mpz_clears(tmp, yi, t, NULL);
}

void montgomery_reduction(mpz_t x, mpz_t N, uint64_t o) {
    mpz_t tmp;
    mpz_init(tmp);

    for (int i = 0; i < N->_mp_size; i++) {
        uint64_t u = x->_mp_d[0] * o;
        mpz_mul_ui(tmp, N, u);
        mpz_add(x, x, tmp);

        mpn_rshift(x->_mp_d, x->_mp_d, x->_mp_size, 32);
        mpn_rshift(x->_mp_d, x->_mp_d, x->_mp_size, 32);

        x->_mp_size -= 1;
    }
    if (mpz_cmp(x, N) >= 0) {
        mpz_sub(x, x, N);
    }

    mpz_clear(tmp);
}

////////////////////////////////////////////////////////////////////////
void montgomery_exp_mod(mpz_t r, mpz_t x, mpz_t y, mpz_t N, mpz_t rr, uint64_t o) {
    mpz_t tmp;
    mpz_t *T;

    T = malloc(sizeof(mpz_t) << (K_BITS - 1));
    mpz_init_set(T[0], x);
    mpz_init(tmp);
    montgomery_mul(tmp, x, x, N, o);
    for (int i = 1; i < 1 << (K_BITS - 1); i++) {
        mpz_init(T[i]);
        montgomery_mul(T[i], T[i - 1], tmp, N, o);
    }

    mpz_set_ui(tmp, 1);
    montgomery_mul(tmp, tmp, rr, N, o);
    int i = y->_mp_size << 6;

    while (i >= 0) {
        int u = 0;
        int l = i - K_BITS + 1;

        if (getBit(y, i)) {
            if (l < 0) l = 0;
            while (!getBit(y, l)) l++;
            int id = i >> 6;
            int ib = i & 63;
            int ld = l >> 6;
            int lb = l & 63;
            if (id == ld) u = (y->_mp_d[id] << (63 - ib)) >> (63 - ib + lb);
            else u = ( y->_mp_d[ld] >> lb) | (((y->_mp_d[id] << (63 - ib)) >> (63 - ib)) << (64 - lb));
        } else {
            l = i;
        }
        for (int j = 0; j < i - l + 1; j++) {
            montgomery_mul(tmp, tmp, tmp, N, o);
        }
        if (u != 0) {
            montgomery_mul(tmp, tmp, T[(u - 1) >> 1], N, o);
        }
        i = l - 1;
    }
    mpz_set(r, tmp);
    for (i = 0; i < 1 << (K_BITS - 1); i++) {
        mpz_clear(T[i]);
    }
    mpz_clear(tmp);
}

/////////////////////////////////
void exp_mod_crt(mpz_t r, mpz_t x, mpz_t y, mpz_t N) {
    mpz_t rr, rrr, x_tmp;
    uint64_t o;

    mpz_init_set(x_tmp, x);

    mpz_init(rr);
    mpz_init(rrr);

    montgomery_init(rr, &o, N);
    montgomery_CRT_init(rr, rrr, N);
    montgomery_mul(x_tmp, x_tmp, rrr, N, o);
    montgomery_reduction(x_tmp, N, o);

    montgomery_exp_mod(r, x_tmp, y, N, rr, o);
    montgomery_reduction(r, N, o);

    mpz_clears(rr, rrr, x_tmp, NULL);
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

void intToStr(char* str, mpz_t num) {
    int loc =0;
    //char* out = malloc(sizeof(char)*255);
    mpz_t pow, quot;
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

    int off = 0;
    for (int i=0; i< 255; i++) {
        if (*str == '0') {
            off++;
            memmove(str, &(str[1]), strlen(&(str[1])));
        }
    }
    str[loc - off] = '\0';
}

void strToInt(mpz_t num, char* str) {
    int pow = 0;
    int n = strcspn(str, "\n\r");
    str[n] = 0;
    mpz_t tmp;
    mpz_init(num);
    mpz_init(tmp);

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
void RSAEncrypt(RSA_public_key *pk, mpz_t cipher) {
    uint64_t o = 0;
    mpz_t rr;
    mpz_inits(cipher, rr, NULL);
    montgomery_init(rr, &o, pk->N);
    montgomery_mul(pk->m, pk->m, rr, pk->N, o);
    montgomery_exp_mod(cipher, pk->m, pk->e, pk->N, rr, o);
    montgomery_reduction(cipher, pk->N, o);
    //mpz_powm(cipher, pk->m, pk->e, pk->N);
}

//Uses CRT
void RSADecrypt(RSA_private_key *sk, mpz_t message) {
    mpz_t m1, m2;
    mpz_init(m1);
    mpz_init(m2);
    mpz_init(message);

    /////////////////////////////
    // m_p = c^d (mod p)
    exp_mod_crt(m1, sk->c, sk->d_p, sk->p);
    ///////////////////
    // m_q = c^d (mod q)
    exp_mod_crt(m2, sk->c, sk->d_q, sk->q);

    //////////////////
    // m = (m_p * q * q^-1 (mod p)) + (m_q * p * p^-1 (mod q)) (mod N)
    mpz_mul(m1, m1, sk->q);
    mpz_mul(m1, m1, sk->i_q);
    mpz_mul(m2, m2, sk->p);
    mpz_mul(m2, m2, sk->i_p);
    mpz_add(message, m1, m2);
    mpz_mod(message, message, sk->N);
}

void ElGamalEncrypt(ElGamal_public_key *pk, mpz_t c1, mpz_t c2) {
    mpz_t rr;
    uint64_t o;
    mpz_inits(pk->k, rr, c1, c2, NULL);

    //if (genRandomKey(pk->k, P_MOD_SIZE) < 0) {
    //    return;
    //}
    //mpz_mod(pk->k, pk->k, pk->p);

    //Enable for testing
    mpz_set_ui(pk->k, 1);

    /////////////////////////////////////
    // c_1 = g^w (mod p)
    // c_2 = m * h^w (mod p)
    // //////////////////////////////////////
    montgomery_init(rr, &o, pk->p);
    montgomery_mul(pk->g, pk->g, rr, pk->p, o);
    montgomery_exp_mod(c1, pk->g, pk->k, pk->p, rr, o);
    montgomery_reduction(c1, pk->p, o);

    montgomery_mul(pk->h, pk->h, rr, pk->p, o);
    montgomery_exp_mod(c2, pk->h, pk->k, pk->p, rr, o);
    montgomery_mul(pk->m, pk->m, rr, pk->p, o);
    montgomery_mul(c2, c2, pk->m, pk->p, o);
    montgomery_reduction(c2, pk->p, o);
}

void ElGamalDecryption(ElGamal_private_key *sk, mpz_t message) {
    mpz_t rr, tmp;
    uint64_t o;
    mpz_inits(rr, tmp, message, NULL);

    //////////////////////////////////
    // m = c_1^(q-x) * c_2 (mod p)
    montgomery_init(rr, &o, sk->p);
    mpz_sub(tmp, sk->q, sk->x);


    montgomery_mul(sk->c1, sk->c1, rr, sk->p, o);
    montgomery_exp_mod(message, sk->c1, tmp, sk->p, rr, o);
    montgomery_mul(sk->c2, sk->c2, rr, sk->p, o);
    montgomery_mul(message, message, sk->c2, sk->p, o);
    montgomery_reduction(message, sk->p, o);
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
        RSA_public_key *pk;
        init_RSA_pk(&pk);

        mpz_set(pk->N, fields[0]);
        mpz_set(pk->e, fields[1]);
        mpz_set(pk->m, fields[2]);

        RSAEncrypt(pk, cipher);

        char* out = malloc(sizeof(char)*255);
        intToStr(out, cipher);
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
        RSA_private_key *sk;
        init_RSA_sk(&sk);

        mpz_set(sk->N, fields[0]);
        mpz_set(sk->d, fields[1]);
        mpz_set(sk->p, fields[2]);
        mpz_set(sk->q, fields[3]);
        mpz_set(sk->d_p, fields[4]);
        mpz_set(sk->d_q, fields[5]);
        mpz_set(sk->i_p, fields[6]);
        mpz_set(sk->i_q, fields[7]);
        mpz_set(sk->c, fields[8]);

        RSADecrypt(sk, message);

        char* out = malloc(sizeof(char)*255);
        intToStr(out, message);
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
        mpz_t c1, c2;
        ElGamal_public_key *pk;
        init_ElGamal_pk(&pk);

        mpz_set(pk->p, fields[0]);
        mpz_set(pk->q, fields[1]);
        mpz_set(pk->g, fields[2]);
        mpz_set(pk->h, fields[3]);
        mpz_set(pk->m, fields[4]);

        ElGamalEncrypt(pk, c1, c2);

        char* out = malloc(sizeof(char)*255);
        intToStr(out, c1);
        fprintf(stdout, "%s\n", out);
        intToStr(out, c2);
        fprintf(stdout, "%s\n", out);
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
        ElGamal_private_key *sk;
        init_ElGamal_sk(&sk);

        mpz_set(sk->p, fields[0]);
        mpz_set(sk->q, fields[1]);
        mpz_set(sk->g, fields[2]);
        mpz_set(sk->x, fields[3]);
        mpz_set(sk->c1, fields[4]);
        mpz_set(sk->c2, fields[5]);

        ElGamalDecryption(sk, message);


        char* out = malloc(sizeof(char)*255);
        intToStr(out, message);
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

    if ( !strcmp( argv[ 1 ], "stage1" ) ) {
        stage1();
    } else if( !strcmp( argv[ 1 ], "stage2" ) ) {
        stage2();
    } else if( !strcmp( argv[ 1 ], "stage3" ) ) {
        stage3();
    } else if( !strcmp( argv[ 1 ], "stage4" ) ) {
        stage4();
    } else {
        abort();
    }

    return 0;
}
