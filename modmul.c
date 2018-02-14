
///////////////////////////////////////////////////////////
//                                                       //
//                 Joshua Van Leeuwen                    //
//                                                       //
//                University of Bristol                  //
//                                                       //
///////////////////////////////////////////////////////////

#include "modmul.h"

#define WINDOW_SIZE 6
#define HEX_BASE 16
#define NUMBER_OF_LIMBS 32
#define LIMB_SIZE 64
#define WORD_LENGTH 256
#define LAMDA_SIZE 1024

int max(int x, int y);
int generate_random_number(mpz_t seed, const size_t size);
int readGroup(int size, mpz_t field[size]);
int hexToZ(char ch);
char zToHex(int n);
void zToStr(char* str, mpz_t num);
void strToZ(mpz_t num, char* str);
void print_number(mpz_t x);

void init_RSA_pk(RSA_public_key **pk, mpz_t fields[3]);
void init_RSA_sk(RSA_private_key **sk, mpz_t fields[9]);
void init_ElGamal_pk(ElGamal_public_key **pk, mpz_t fields[5]);
void init_ElGamal_sk(ElGamal_private_key **sk, mpz_t fields[6]);

void RSAEncrypt(RSA_public_key *pk, mpz_t cipher);
void RSADecrypt(RSA_private_key *sk, mpz_t message);
void ElGamalEncrypt(ElGamal_public_key *pk, mpz_t c1, mpz_t c2, mpz_t key);
void ElGamalDecryption(ElGamal_private_key *sk, mpz_t message);

void mnt_ro2(mpz_t ro2, mpz_t N);
void mnt_omega(mp_limb_t *o, mpz_t N);
void mnt_mul(mpz_t r, mpz_t x, mpz_t y, mpz_t N, mp_limb_t o);
void mnt_pow_mod(mpz_t r, mpz_t x, mpz_t y, mpz_t N, mpz_t ro2, mp_limb_t o);
void mnt_red(mpz_t r, mpz_t t, mpz_t N, mp_limb_t o);

int BBS_init(BBS* bbs);
int BBS_next(BBS *bbs);
int BBS_check_prime(mpz_t p);

void stage1();
void stage2();
void stage3();
void stage4();

void init_RSA_pk(RSA_public_key **pk, mpz_t fields[3]) {
    *pk = (RSA_public_key*) malloc(sizeof(RSA_public_key));
    mpz_inits((*pk)->N, (*pk)->e, (*pk)->ro2,  NULL);
    mpz_init((*pk)->m);
    (*pk)->o = malloc(sizeof(mp_limb_t));

    mpz_set((*pk)->N, fields[0]);
    mpz_set((*pk)->e, fields[1]);
    mpz_set((*pk)->m, fields[2]);

    mnt_ro2((*pk)->ro2, (*pk)->N);
    mnt_omega((*pk)->o, (*pk)->N);
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
    mpz_inits((*pk)->p, (*pk)->q, (*pk)->g, (*pk)->h, (*pk)->ro2, NULL);
    mpz_init( (*pk)->m);
    (*pk)->o = malloc(sizeof(mp_limb_t));

    mpz_set((*pk)->p, fields[0]);
    mpz_set((*pk)->q, fields[1]);
    mpz_set((*pk)->g, fields[2]);
    mpz_set((*pk)->h, fields[3]);
    mpz_set((*pk)->m, fields[4]);

    mnt_ro2((*pk)->ro2, (*pk)->p);
    mnt_omega((*pk)->o, (*pk)->p);
}

void init_ElGamal_sk(ElGamal_private_key **sk, mpz_t fields[6]) {
    *sk = (ElGamal_private_key*) malloc(sizeof(ElGamal_private_key));
    mpz_inits((*sk)->p,  (*sk)->q,  (*sk)->g, (*sk)->x, (*sk)->ro2, NULL);
    mpz_inits((*sk)->c1, (*sk)->c2, NULL);
    (*sk)->o = malloc(sizeof(mp_limb_t));

    mpz_set((*sk)->p, fields[0]);
    mpz_set((*sk)->q, fields[1]);
    mpz_set((*sk)->g, fields[2]);
    mpz_set((*sk)->x, fields[3]);
    mpz_set((*sk)->c1, fields[4]);
    mpz_set((*sk)->c2, fields[5]);

    mnt_ro2((*sk)->ro2, (*sk)->p);
    mnt_omega((*sk)->o, (*sk)->p);
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

void mnt_omega(mp_limb_t *o, mpz_t N) {
    *o = 1;
    for (int i = 1; i < LIMB_SIZE; i++) {
        *o *= *o * N->_mp_d[0];
    }
    *o = -*o;
}

void mnt_ro2(mpz_t ro2, mpz_t N) {
    mpz_init(ro2);
    mpz_set_ui(ro2, 1);

    for (int i=1; i <= N->_mp_size * LIMB_SIZE * 2; i++){
        mpz_add(ro2, ro2, ro2);

        if (mpz_cmp(N, ro2) < 0) {
            mpz_sub(ro2, ro2, N);
        }
    }
}

void mnt_mul(mpz_t r, mpz_t x, mpz_t y, mpz_t N, mp_limb_t o) {
    mpz_t yix, uiN, res;
    mpz_inits(res, yix, uiN, NULL);
    mpz_set_ui(res, 0);

    for (int i = 0; i < N->_mp_size; i++) {
        mp_limb_t ui = (res->_mp_d[0] + (y->_mp_d[i] * x->_mp_d[0])) * o;

        mpz_mul_ui(yix, x, y->_mp_d[i]);
        mpz_mul_ui(uiN, N, ui);

        mpz_add(res, res, yix);
        mpz_add(res, res, uiN);

        for (int j = 0; j < res->_mp_size - 1; j++) {
            res->_mp_d[j] = res->_mp_d[j + 1];
        }
        res->_mp_size -= 1;
    }

    if (mpz_cmp(res, N) >= 0) {
        mpz_sub(res, res, N);
    }

    mpz_set(r, res);
}

void mnt_red(mpz_t r, mpz_t t, mpz_t N, mp_limb_t o) {
    mpz_t uiN, bi;
    mpz_inits(uiN, bi, NULL);
    mpz_set(r, t);

    // Divide by b at each iteration
    for (int i = 0; i < N->_mp_size; i++) {
        mp_limb_t ui = r->_mp_d[0] * o;

        mpz_mul_ui(uiN, N, ui);
        mpz_add(r, r, uiN);

        for (int j = 0; j < r->_mp_size-1; j++) {
            r->_mp_d[j] = r->_mp_d[j + 1];
        }
        r->_mp_size -= 1;
    }

    if (mpz_cmp(r, N) >= 0) {
        mpz_sub(r, r, N);
    }
}

void mnt_pow_mod(mpz_t r, mpz_t x, mpz_t y, mpz_t N, mpz_t ro2, mp_limb_t o) {
    mpz_t T[NUMBER_OF_LIMBS], t, tmp;
    mpz_inits(t, tmp, NULL);

    // initiate T
    mnt_mul(tmp, x, x, N, o);
    mpz_init_set(T[0], x);
    for (int i = 1; i < NUMBER_OF_LIMBS; i++) {
        mpz_init(T[i]);
        mnt_mul(T[i], T[i - 1], tmp, N, o);
    }

    // initiate t
    mpz_set_ui(t, 1);
    mnt_mul(t, t, ro2, N, o);

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
            // If not then get bits from each limb appropriately
            if (i_limb == l_limb) {
                u = get_word(y->_mp_d[i_limb], l_bits, i_bits);
            } else {
                u = (get_word(y->_mp_d[i_limb], 0, i_bits) << (LIMB_SIZE - l_bits)) | get_word(y->_mp_d[l_limb], l_bits, LIMB_SIZE - 1);
            }
        }

        for (int j = 0; j <= i - l; j++) {
            mnt_mul(t, t, t, N, o);
        }

        if (u != 0) {
            mnt_mul(t, t, T[(u-1)/2], N, o);
        }

        i = l - 1;
    }

    mpz_set(r, t);
}

int hexToZ(char ch) {
    if (ch >= '0' && ch <= '9') {
        return ch - '0';
    }

    if (ch >= 'A' && ch <= 'F') {
        return ch - 'A' + 10;
    }

    if (ch >= 'a' && ch <= 'f') {
        return ch - 'a' + 10;
    }

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

int generate_random_number(mpz_t seed, const size_t size) {
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

    mpz_import(seed, size, 1, sizeof(char), 0, 0, dat);

    return 0;
}

void RSAEncrypt(RSA_public_key *pk, mpz_t cipher) {
    mpz_inits(cipher, NULL);

    //c = m^e mod (N)
    mnt_mul(pk->m, pk->m, pk->ro2, pk->N, *pk->o);
    mnt_pow_mod(cipher, pk->m, pk->e, pk->N, pk->ro2, *pk->o);
    mnt_red(cipher, cipher, pk->N, *pk->o);
}

// Uses Chinese Remainder Theorem
void RSADecrypt(RSA_private_key *sk, mpz_t message) {
    mpz_t m1, m2;
    mpz_inits(m1, m2, message, NULL);

    mpz_powm(m1, sk->c, sk->d_p, sk->p);
    mpz_powm(m2, sk->c, sk->d_q, sk->q);

    // qInv = (1/q) mod p
    // h = qInv.(m1 - m2) mod p
    // m = m2 + h.q
    mpz_sub(m1, m1, m2);
    mpz_mul(m1, sk->i_q, m1);
    mpz_mod(m1, m1, sk->p);
    mpz_mul(m1, m1, sk->q);
    mpz_add(message, m1, m2);
}

void ElGamalEncrypt(ElGamal_public_key *pk, mpz_t c1, mpz_t c2, mpz_t key) {
    mpz_inits(pk->k, c1, c2, NULL);

    if(mpz_cmp(key, pk->p) <= 0) {
        fprintf(stderr, "Warning: ElGamal Encrypt key may not be large enough.\n");
    } else {
        mpz_mod(key, key, pk->p);
    }
    mpz_set(pk->k, key);

    //Enable for testing
#if defined(DEBUG)
    mpz_set_ui(pk->k, 1);
#endif

    mnt_mul(pk->g, pk->g, pk->ro2, pk->p, *pk->o);
    mnt_pow_mod(c1, pk->g, pk->k, pk->p, pk->ro2, *pk->o);

    mnt_mul(pk->h, pk->h, pk->ro2, pk->p, *pk->o);
    mnt_pow_mod(c2, pk->h, pk->k, pk->p, pk->ro2, *pk->o);
    mnt_mul(pk->m, pk->m, pk->ro2, pk->p, *pk->o);
    mnt_mul(c2, c2, pk->m, pk->p, *pk->o);

    mnt_red(c1, c1, pk->p, *pk->o);
    mnt_red(c2, c2, pk->p, *pk->o);
}

void ElGamalDecryption(ElGamal_private_key *sk, mpz_t message) {
    mpz_t qx;
    mpz_inits(qx, message, NULL);

    mpz_sub(qx, sk->q, sk->x);

    mnt_mul(sk->c1, sk->c1, sk->ro2, sk->p, *sk->o);
    mnt_pow_mod(message, sk->c1, qx, sk->p, sk->ro2, *sk->o);

    mnt_mul(sk->c2, sk->c2, sk->ro2, sk->p, *sk->o);
    mnt_mul(message, message, sk->c2, sk->p, *sk->o);

    mnt_red(message, message, sk->p, *sk->o);
}

int BBS_check_prime(mpz_t p) {
    if (mpz_probab_prime_p(p, 100) == 0) {
        return 0;
    }

    mpz_t test;
    mpz_init(test);
    mpz_mod_ui(test, p, 4);
    if (mpz_cmp_ui(test, 3) == 0) {
        return 1;
    }

    return 0;
}

int BBS_next(BBS *bbs) {
    mnt_mul(bbs->s, bbs->s, bbs->ro2, bbs->N, *bbs->o);
    mnt_pow_mod(bbs->s, bbs->s, bbs->two, bbs->N, bbs->ro2, *bbs->o);
    mnt_red(bbs->s, bbs->s, bbs->N, *bbs->o);

    return (int) (bbs->s->_mp_d[0] & 1);
}

int BBS_init(BBS* bbs) {
    mpz_t p, q, N, s;
    mpz_t gcd, ro2;
    mpz_inits(p, q, N, s, gcd, ro2, NULL);
    mpz_inits(bbs->p, bbs->q, bbs->N, bbs->s, bbs->ro2, bbs->two, NULL);
    mpz_set_ui(p, 0);
    mpz_set_ui(q, 0);
    mpz_set_ui(gcd, 0);
    mpz_set_ui(bbs->two, 2);

    while(BBS_check_prime(p) == 0) {
        if(generate_random_number(p, sizeof(mp_limb_t)*NUMBER_OF_LIMBS) < 0) {
            return -1;
        }
        mpz_nextprime(p, p);
    }

    while(BBS_check_prime(q) == 0) {
        if (generate_random_number(q, sizeof(mp_limb_t)*NUMBER_OF_LIMBS) < 0) {
            return -1;
        }
        mpz_nextprime(q, q);
    }

    mpz_mul(N, p, q);

    while(!(mpz_cmp_ui(gcd, 1) == 0)) {
        if (generate_random_number(s, sizeof(mp_limb_t)*NUMBER_OF_LIMBS) < 0) {
            return -1;
        }
        mpz_mod(s, s, N);
        if (mpz_cmp(s, N) >= 0) {
            mpz_sub(s, s, N);
        }

        mpz_gcd(gcd, s, N);
    }

    bbs->o = malloc(sizeof(mp_limb_t));
    mnt_omega(bbs->o, N);
    mnt_ro2(ro2, N);

    mpz_set(bbs->p, p);
    mpz_set(bbs->q, q);
    mpz_set(bbs->N, N);
    mpz_set(bbs->s, s);
    mpz_set(bbs->ro2, ro2);

    BBS_next(bbs);

    return 0;
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

    BBS bbs;
    if (BBS_init(&bbs) < 0) {
        return;
    }

    while(readGroup(exp_size, fields) != -1) {
        BBS_next(&bbs);

        mpz_t c1, c2;
        ElGamal_public_key *pk;
        init_ElGamal_pk(&pk, fields);

        ElGamalEncrypt(pk, c1, c2, bbs.s);

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
