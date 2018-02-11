// Joshua Van Leeuwen
// University of Bristol

#include "modmul.h"

#define WINDOW_SIZE 6
#define HEX_BASE 16
#define WORD_LENGTH 256
#define LIMB_SIZE 64

void init_RSA_pk(RSA_public_key **pk, mpz_t fields[3]) {
    *pk = (RSA_public_key*) malloc(sizeof(RSA_public_key));
    mpz_inits((*pk)->N, (*pk)->e,  NULL);
    mpz_init((*pk)->m);

    mpz_set((*pk)->N, fields[0]);
    mpz_set((*pk)->e, fields[1]);
    mpz_set((*pk)->m, fields[2]);
}

void init_RSA_sk(RSA_private_key **sk, mpz_t fields[9]) {
    *sk = (RSA_private_key*) malloc(sizeof(RSA_private_key));
    mpz_inits((*sk)->N,   (*sk)->d,   (*sk)->p,   (*sk)->q,   NULL);
    mpz_inits((*sk)->d_p, (*sk)->d_q, (*sk)->i_p, (*sk)->i_q, NULL);
    mpz_init( (*sk)->c);

    mpz_set((*sk)->N, fields[0]);
    mpz_set((*sk)->d, fields[1]);
    mpz_set((*sk)->p, fields[2]);
    mpz_set((*sk)->q, fields[3]);
    mpz_set((*sk)->d_p, fields[4]);
    mpz_set((*sk)->d_q, fields[5]);
    mpz_set((*sk)->i_p, fields[6]);
    mpz_set((*sk)->i_q, fields[7]);
    mpz_set((*sk)->c, fields[8]);

}

void init_ElGamal_pk(ElGamal_public_key **pk, mpz_t fields[5]) {
    *pk = (ElGamal_public_key*) malloc(sizeof(ElGamal_public_key));
    mpz_inits((*pk)->p, (*pk)->q, (*pk)->g, (*pk)->h, NULL);
    mpz_init( (*pk)->m);

    mpz_set((*pk)->p, fields[0]);
    mpz_set((*pk)->q, fields[1]);
    mpz_set((*pk)->g, fields[2]);
    mpz_set((*pk)->h, fields[3]);
    mpz_set((*pk)->m, fields[4]);
}

void init_ElGamal_sk(ElGamal_private_key **sk, mpz_t fields[6]) {
    *sk = (ElGamal_private_key*) malloc(sizeof(ElGamal_private_key));
    mpz_inits((*sk)->p,  (*sk)->q,  (*sk)->g, (*sk)->x, NULL);
    mpz_inits((*sk)->c1, (*sk)->c2, NULL);

    mpz_set((*sk)->p, fields[0]);
    mpz_set((*sk)->q, fields[1]);
    mpz_set((*sk)->g, fields[2]);
    mpz_set((*sk)->x, fields[3]);
    mpz_set((*sk)->c1, fields[4]);
    mpz_set((*sk)->c2, fields[5]);
}

int max(int x, int y) {
    if (x > y) {
        return x;
    }

    return y;
}

int bit_at_position(mpz_t x, int n) {
    int limb = n / LIMB_SIZE;
    int pos = n % LIMB_SIZE;
    return ((x->_mp_d[limb] >> pos) & 1);
}

int get_word(mp_limb_t x, int start, int end) {
    return (x << (LIMB_SIZE - 1 - end)) >> (LIMB_SIZE - 1 - end + start);
}

void mont_omega(mp_limb_t *o, mpz_t N) {
    *o = 1;
    for (int i = 1; i < LIMB_SIZE; i++) {
        *o *= *o * N->_mp_d[0];
    }
    *o = -*o;
}

void mont_ro(mpz_t ro, mpz_t N, int power) {
    mpz_init(ro);
    mpz_set_ui(ro, 1);

    for (int i=1; i <= N->_mp_size * LIMB_SIZE * power; i++){
        mpz_add(ro, ro, ro);
        if (mpz_cmp(ro, N) >= 0) {
            mpz_sub(ro, ro, N);
        }
    }
}

void mont_mul(mpz_t r, mpz_t x, mpz_t y, mpz_t N, mp_limb_t o) {
    mpz_t yix, uiN, ans;
    mpz_inits(ans, yix, uiN, NULL);
    mpz_set_ui(ans, 0);

    for (int i=0; i < N->_mp_size; i++) {
        mp_limb_t ui = (ans->_mp_d[0] + (y->_mp_d[i] * x->_mp_d[0])) * o;

        mpz_mul_ui(yix, x, y->_mp_d[i]);
        mpz_mul_ui(uiN, N, ui);

        mpz_add(ans, ans, yix);
        mpz_add(ans, ans, uiN);
        mpn_rshift(ans->_mp_d, ans->_mp_d, ans->_mp_size, LIMB_SIZE);

        ans->_mp_size --;
    }

    if (mpz_cmp(ans, N) >= 0) {
        mpz_sub(ans, ans, N);
    }

    mpz_set(r, ans);
}

void mont_red(mpz_t r, mpz_t t, mpz_t N, mp_limb_t o) {
    mpz_t uiN;
    mpz_init(uiN);
    mpz_set(r, t);

    for (int i=0; i < N->_mp_size; i++) {
        mp_limb_t ui = r->_mp_d[0] * o;

        mpz_mul_ui (uiN, N, ui);
        mpz_add(r, r, uiN);

        mpn_rshift(r->_mp_d, r->_mp_d, r->_mp_size, LIMB_SIZE);
        r->_mp_size -= 1;
    }

    if (mpz_cmp(r, N) >= 0) {
        mpz_sub(r, r, N);
    }
}

void mont_pow_mod(mpz_t r, mpz_t x, mpz_t y, mpz_t N, mpz_t ro2, mp_limb_t o) {
    mpz_t T[32], t, tmp;
    mpz_inits(t, tmp, NULL);

    // initiate T
    mont_mul(tmp, x, x, N, o);
    mpz_init_set(T[0], x);
    for (int i = 1; i < 32; i++) {
        mpz_init(T[i]);
        mont_mul(T[i], T[i - 1], tmp, N, o);
    }

    // initiate t
    mpz_set_ui(t, 1);
    mont_mul(t, t, ro2, N, o);

    // Set i to the number of bits in y - 1
    int i = (y->_mp_size * LIMB_SIZE) - 1;
    int u, l;

    while( i >= 0) {
        if (bit_at_position(y, i) == 0) {
            l = i;
            u = 0;

        } else {

            l = max(i - WINDOW_SIZE + 1, 0);

            while (bit_at_position(y, l) == 0) {
                l++;
            }

            int i_limb = i / LIMB_SIZE;
            int i_bits = i % LIMB_SIZE;
            int l_limb = l / LIMB_SIZE;
            int l_bits = l % LIMB_SIZE;

            // Check whether i and l are in the same limb
            if (i_limb == l_limb) {
                u = get_word(y->_mp_d[i_limb], l_bits, i_bits);
            } else {
                u = (get_word(y->_mp_d[i_limb], 0, i_bits) << (LIMB_SIZE - l_bits)) | get_word(y->_mp_d[l_limb], l_bits, LIMB_SIZE - 1);
            }
        }

        for (int j = 0; j < i - l + 1; j++) {
            mont_mul(t, t, t, N, o);
        }

        if (u != 0) {
            mont_mul(t, t, T[(u-1)/2], N, o);
        }

        i = l - 1;
    }

    mpz_set(r, t);
}

void exp_mod_crt(mpz_t r, mpz_t x, mpz_t y, mpz_t N) {
    mpz_t ro2, ro3, tmp;
    mp_limb_t o;
    mpz_inits(ro2, ro3, tmp, NULL);

    mpz_set(tmp, x);

    mont_omega(&o, N);
    mont_ro(ro2, N, 2);
    mont_ro(ro3, N, 3);

    mont_mul(tmp, tmp, ro3, N, o);
    mont_red(tmp, tmp, N, o);

    mont_pow_mod(r, tmp, y, N, ro2, o);
    mont_red(r, r, N, o);
}


int hexToZ(char ch) {
    if (ch >= '0' && ch <= '9')
        return ch - '0';
    if (ch >= 'A' && ch <= 'F')
        return ch - 'A' + 10;
    if (ch >= 'a' && ch <= 'f')
        return ch - 'a' + 10;
    return -1;
}

char zToHex(int n) {
    if (n > 9) {
        return 'A' + (n - 10);
    }
    return n + '0';
}

void zToStr(char* str, mpz_t num) {
    mpz_t power, quot;
    mpz_inits(power, quot, NULL);

    int loc = 0;
    for(int i=(WORD_LENGTH - 1); i>= 0; i--) {
        mpz_ui_pow_ui(power, HEX_BASE, i);

        if (mpz_cmp(num, power) < 0) {
            str[loc] = '0';

        } else {
            mpz_tdiv_q(quot, num, power);
            mpz_mod(num, num, power);

            int lquote = mpz_get_ui(quot);

            char c = zToHex(lquote);
            str[loc] = c;
        }
        loc++;
    }

    int off = 0;
    for (int i = 0; i < 255; i++) {
        if (*str == '0') {
            off++;
            memmove(str, &(str[1]), strlen(&(str[1])));
        }
    }
    str[loc - off] = '\0';
}

void strToZ(mpz_t num, char* str) {
    int pow = 0;
    int n = strcspn(str, "\n\r");
    str[n] = 0;
    mpz_t tmp;
    mpz_init(num);
    mpz_init(tmp);

    for (int i=(n - 1); i >= 0; i--) {
        mpz_ui_pow_ui(tmp, HEX_BASE, pow);
        mpz_mul_si(tmp, tmp, hexToZ(str[i]));
        mpz_add(num, num, tmp);
        mpz_init(tmp);
        pow++;
    }
}

void print_number(mpz_t x) {
    char* out = malloc(sizeof(char)*255);
    zToStr(out, x);
    fprintf( stdout, "%s\n", out);
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

        strToZ(field[i], line);
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
    mp_limb_t o;
    mpz_t ro;
    mpz_inits(cipher, ro, NULL);

    mont_omega(&o, pk->N);
    mont_ro(ro, pk->N, 2);

    mont_mul(pk->m, pk->m, ro, pk->N, o);
    mont_pow_mod(cipher, pk->m, pk->e, pk->N, ro, o);
    mont_red(cipher, cipher, pk->N, o);
}

//Uses CRT
void RSADecrypt(RSA_private_key *sk, mpz_t message) {
    mpz_t m1, m2;
    mpz_t ro2_p, ro2_q;
    mpz_t ro3_p, ro3_q;
    mp_limb_t o_p, o_q;
    mpz_inits(m1, m2, message, ro2_p, ro2_q, ro3_p, ro3_q, NULL);

    mont_omega(&o_p, sk->p);
    mont_omega(&o_q, sk->q);
    mont_ro(ro2_p, sk->p, 2);
    mont_ro(ro3_p, sk->p, 3);
    mont_ro(ro2_q, sk->q, 2);
    mont_ro(ro3_q, sk->q, 3);

    mont_mul(m1, sk->c, ro3_p, sk->p, o_p);
    mont_red(m1, m1, sk->p, o_p);
    mont_pow_mod(m1, m1, sk->d_p, sk->p, ro2_p, o_p);
    mont_red(m1, m1, sk->p, o_p);

    mont_mul(m2, sk->c, ro3_q, sk->q, o_q);
    mont_red(m2, m2, sk->q, o_q);
    mont_pow_mod(m2, m2, sk->d_q, sk->q, ro2_q, o_q);
    mont_red(m2, m2, sk->q, o_q);

    mpz_mul(m1, m1, sk->q);
    mpz_mul(m1, m1, sk->i_q);
    mpz_mul(m2, m2, sk->p);
    mpz_mul(m2, m2, sk->i_p);

    mpz_add(message, m1, m2);
    mpz_mod(message, message, sk->N);
}

void ElGamalEncrypt(ElGamal_public_key *pk, mpz_t c1, mpz_t c2) {
    mpz_t ro2;
    mp_limb_t o;
    mpz_inits(pk->k, ro2, c1, c2, NULL);

    //if (genRandomKey(pk->k, P_MOD_SIZE) < 0) {
    //    return;
    //}
    //mpz_mod(pk->k, pk->k, pk->p);

    //Enable for testing
    mpz_set_ui(pk->k, 1);

    mont_omega(&o, pk->p);
    mont_ro(ro2, pk->p, 2);

    mont_mul(pk->g, pk->g, ro2, pk->p, o);
    mont_pow_mod(c1, pk->g, pk->k, pk->p, ro2, o);

    mont_mul(pk->h, pk->h, ro2, pk->p, o);
    mont_pow_mod(c2, pk->h, pk->k, pk->p, ro2, o);
    mont_mul(pk->m, pk->m, ro2, pk->p, o);
    mont_mul(c2, c2, pk->m, pk->p, o);

    mont_red(c1, c1, pk->p, o);
    mont_red(c2, c2, pk->p, o);
}

void ElGamalDecryption(ElGamal_private_key *sk, mpz_t message) {
    mp_limb_t o;
    mpz_t ro2, qx;
    mpz_inits(ro2, qx, message, NULL);

    mont_omega(&o, sk->p);
    mont_ro(ro2, sk->p, 2);
    mpz_sub(qx, sk->q, sk->x);

    mont_mul(sk->c1, sk->c1, ro2, sk->p, o);
    mont_pow_mod(message, sk->c1, qx, sk->p, ro2, o);

    mont_mul(sk->c2, sk->c2, ro2, sk->p, o);
    mont_mul(message, message, sk->c2, sk->p, o);

    mont_red(message, message, sk->p, o);
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
        init_RSA_pk(&pk, fields);

        RSAEncrypt(pk, cipher);

        print_number(cipher);
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
        init_RSA_sk(&sk, fields);

        RSADecrypt(sk, message);

        print_number(message);
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
        init_ElGamal_pk(&pk, fields);

        ElGamalEncrypt(pk, c1, c2);

        print_number(c1);
        print_number(c2);
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
        init_ElGamal_sk(&sk, fields);

        ElGamalDecryption(sk, message);

        print_number(message);
    }
}

int main( int argc, char* argv[] ) {
    if( 2 != argc ) {
        fprintf(stderr, "Expected 1 argument; got=%d\n", argc-1);
        return 1;
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
        fprintf(stderr, "Unknown argument: %s\n", argv[1]);
        return 1;
    }

    return 0;
}
