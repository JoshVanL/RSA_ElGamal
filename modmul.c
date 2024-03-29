
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

int max(int x, int y);
int readGroup(int size, mpz_t field[size]);
int hexToZ(char ch);
char zToHex(int n);
void zToStr(char* str, mpz_t num);
void strToZ(mpz_t num, char* str);
void print_number(mpz_t x);
int concanonate(int x, int y);

void init_RSA_pk(RSA_public_key **pk, mpz_t fields[3]);
void init_RSA_sk(RSA_private_key **sk, mpz_t fields[9]);
void init_ElGamal_pk(ElGamal_public_key **pk, mpz_t fields[5]);
void init_ElGamal_sk(ElGamal_private_key **sk, mpz_t fields[6]);

void RSAEncrypt(RSA_public_key *pk, mpz_t cipher);
void RSADecrypt(RSA_private_key *sk, mpz_t message);
void ElGamalEncrypt(ElGamal_public_key *pk, mpz_t c1, mpz_t c2, mpz_t key);
void ElGamalDecryption(ElGamal_private_key *sk, mpz_t message);

void mnt_ro2(mpz_t ro2, mpz_t N);
void mnt_omega(mp_limb_t* o, mpz_t N);
void mnt_mul(mpz_t r, mpz_t x, mpz_t y, mpz_t N, mp_limb_t o);
void mnt_pow_mod(mpz_t r, mpz_t x, mpz_t y, mpz_t N, mpz_t ro2, mp_limb_t o);
void mnt_red(mpz_t r, mpz_t t, mpz_t N, mp_limb_t o);

int BBS_init(BBS* bbs);
int BBS_next(BBS* bbs);
int BBS_check_prime(mpz_t p);
int generate_random_number(mpz_t seed, const size_t size);

void stage1();
void stage2();
void stage3();
void stage4();

// Initiate RSA public key
void init_RSA_pk(RSA_public_key **pk, mpz_t fields[3]) {
    *pk = (RSA_public_key*) malloc(sizeof(RSA_public_key));
    mpz_init((*pk)->N);
    mpz_init((*pk)->e);
    mpz_init((*pk)->ro2);
    mpz_init((*pk)->m);
    (*pk)->o = malloc(sizeof(mp_limb_t));

    mpz_set((*pk)->N, fields[0]);
    mpz_set((*pk)->e, fields[1]);
    mpz_set((*pk)->m, fields[2]);

    mnt_ro2((*pk)->ro2, (*pk)->N);
    mnt_omega((*pk)->o, (*pk)->N);
}

// Initiate RSA secret key
void init_RSA_sk(RSA_private_key **sk, mpz_t fields[9]) {
    *sk = (RSA_private_key*) malloc(sizeof(RSA_private_key));
    mpz_init((*sk)->N);
    mpz_init((*sk)->d);
    mpz_init((*sk)->p);
    mpz_init((*sk)->q);
    mpz_init((*sk)->d_p);
    mpz_init((*sk)->d_q);
    mpz_init((*sk)->i_p);
    mpz_init((*sk)->i_q);
    mpz_init((*sk)->ro2_p);
    mpz_init((*sk)->ro2_q);
    mpz_init( (*sk)->c);
    (*sk)->o_p = malloc(sizeof(mp_limb_t));
    (*sk)->o_q = malloc(sizeof(mp_limb_t));

    mpz_set((*sk)->N,   fields[0]);
    mpz_set((*sk)->d,   fields[1]);
    mpz_set((*sk)->p,   fields[2]);
    mpz_set((*sk)->q,   fields[3]);
    mpz_set((*sk)->d_p, fields[4]);
    mpz_set((*sk)->d_q, fields[5]);
    mpz_set((*sk)->i_p, fields[6]);
    mpz_set((*sk)->i_q, fields[7]);
    mpz_set((*sk)->c,   fields[8]);

    mnt_ro2((*sk)->ro2_p, (*sk)->p);
    mnt_ro2((*sk)->ro2_q, (*sk)->q);
    mnt_omega((*sk)->o_p, (*sk)->p);
    mnt_omega((*sk)->o_q, (*sk)->q);
}

// Initiate ElGamal public key
void init_ElGamal_pk(ElGamal_public_key **pk, mpz_t fields[5]) {
    *pk = (ElGamal_public_key*) malloc(sizeof(ElGamal_public_key));
    mpz_init((*pk)->p);
    mpz_init((*pk)->q);
    mpz_init((*pk)->g);
    mpz_init((*pk)->h);
    mpz_init((*pk)->ro2);
    mpz_init((*pk)->m);
    (*pk)->o = malloc(sizeof(mp_limb_t));

    mpz_set((*pk)->p, fields[0]);
    mpz_set((*pk)->q, fields[1]);
    mpz_set((*pk)->g, fields[2]);
    mpz_set((*pk)->h, fields[3]);
    mpz_set((*pk)->m, fields[4]);

    mnt_ro2((*pk)->ro2, (*pk)->p);
    mnt_omega((*pk)->o, (*pk)->p);
}

// Initiate ElGamal secret key
void init_ElGamal_sk(ElGamal_private_key **sk, mpz_t fields[6]) {
    *sk = (ElGamal_private_key*) malloc(sizeof(ElGamal_private_key));
    mpz_init((*sk)->p);
    mpz_init((*sk)->q);
    mpz_init((*sk)->g);
    mpz_init((*sk)->x);
    mpz_init((*sk)->ro2);
    mpz_init((*sk)->c1);
    mpz_init((*sk)->c2);
    (*sk)->o = malloc(sizeof(mp_limb_t));

    mpz_set((*sk)->p,  fields[0]);
    mpz_set((*sk)->q,  fields[1]);
    mpz_set((*sk)->g,  fields[2]);
    mpz_set((*sk)->x,  fields[3]);
    mpz_set((*sk)->c1, fields[4]);
    mpz_set((*sk)->c2, fields[5]);

    mnt_ro2((*sk)->ro2, (*sk)->p);
    mnt_omega((*sk)->o, (*sk)->p);
}

// Calculate max of two ints
int max(int x, int y) {
    if (x > y) {
        return x;
    }

    return y;
}

// Get the bit of a given position (0|1)
int bit_at_position(mpz_t x, int n) {
    int limb = n / LIMB_SIZE;
    int pos = n % LIMB_SIZE;
    return ((x->_mp_d[limb] >> pos) & 1);
}

// Get the bits of a given limb by some range
int get_word(mp_limb_t x, int start, int end) {
    return (x << (LIMB_SIZE - 1 - end)) >> (LIMB_SIZE - 1 - end + start);
}

// bitwise or of two ints
int concanonate(int x, int y) {
    return x | y;
}

// Calculate the omega (Montgomery)
void mnt_omega(mp_limb_t *o, mpz_t N) {
    *o = 1;
    for (int i = 1; i < LIMB_SIZE; i++) {
        *o *= *o * N->_mp_d[0];
    }
    *o = -*o;
}

// Calculate rho for Montgomery (squared)
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

// Perform Montgomery Multiplication
void mnt_mul(mpz_t r, mpz_t x, mpz_t y, mpz_t N, mp_limb_t o) {
    mpz_t yix, uiN, res;
    mpz_init(res);
    mpz_init(yix);
    mpz_init(uiN);
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

// Perform Montgomery Reduction
void mnt_red(mpz_t r, mpz_t t, mpz_t N, mp_limb_t o) {
    mpz_t uiN, bi;
    mpz_init(uiN);
    mpz_init(bi);
    mpz_set(r, t);

    // Divide by b at each iteration instead of after loop
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

// Perform Montgomery exponential using "windowed" exponentiation
void mnt_pow_mod(mpz_t r, mpz_t x, mpz_t y, mpz_t N, mpz_t ro2, mp_limb_t o) {
    mpz_t T[NUMBER_OF_LIMBS], t_hat, x0, x_hat;
    mpz_init(t_hat);
    mpz_init(x_hat);
    mpz_init(x0);

    // initiate t_hat
    mpz_set_ui(t_hat, 1);
    mnt_mul(t_hat, t_hat, ro2, N, o);

    // initiate x_hat
    mnt_mul(x_hat, x, ro2, N, o);

    // initiate T
    mnt_mul(x0, x_hat, x_hat, N, o);
    mpz_init_set(T[0], x_hat);
    for (int i = 1; i < NUMBER_OF_LIMBS; i++) {
        mpz_init(T[i]);
        mnt_mul(T[i], T[i - 1], x0, N, o);
    }

    // Set i to the number of bits in y - 1
    int i = (y->_mp_size * LIMB_SIZE) - 1;
    int u, l;

    while( i >= 0 ) {
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
            if ( i_limb == l_limb ) {
                u = get_word(y->_mp_d[i_limb], l_bits, i_bits);
            } else {
                int ii = get_word(y->_mp_d[i_limb], 0, i_bits) << (LIMB_SIZE - l_bits);
                int ll = get_word(y->_mp_d[l_limb], l_bits, LIMB_SIZE - 1);
                u = concanonate(ii, ll);
            }
        }

        for (int j = 0; j < i - l + 1; j++) {
            mnt_mul(t_hat, t_hat, t_hat, N, o);
        }

        if (u != 0) {
            mnt_mul(t_hat, t_hat, T[(u-1)/2], N, o);
        }

        i = l - 1;
    }

    mpz_set(r, t_hat);
}

// Convert a hexadecimal char to it's integer representation
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

// Convert a number to it's hexadecimal representation
char zToHex(int n) {
    if (n > 9) {
        return 'A' + (n - 10);
    }
    return n + '0';
}

// Convert a number to it's hexadecimal string representation
void zToStr(char* str, mpz_t num) {
    mpz_t power, quot;
    mpz_init(power);
    mpz_init(quot);

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

// Convert a hexadecimal string to it's number representation
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

// Print a number to Stdout in it's hexadecimal representation
void print_number(mpz_t x) {
    char* out = malloc(sizeof(char)*255);
    zToStr(out, x);
    fprintf( stdout, "%s\n", out);
}


// Read a size number of hexadecimal lines
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

// Generate a size number of random bits, reading from /dev/random
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

// Perform RSA encryption
void RSAEncrypt(RSA_public_key *pk, mpz_t cipher) {
    mpz_init(cipher);

    //c = m^e mod (N)
    mnt_pow_mod(cipher, pk->m, pk->e, pk->N, pk->ro2, *pk->o);
    mnt_red(cipher, cipher, pk->N, *pk->o);
}

// Perform RSA decryption using the Chinese Remainder Theorem (CRT)
void RSADecrypt(RSA_private_key *sk, mpz_t message) {
    mpz_t m1, m2;
    mpz_init(m1);
    mpz_init(m2);
    mpz_init(message);

    mpz_mod(m1, sk->c, sk->p);
    mnt_pow_mod(m1, m1, sk->d_p, sk->p, sk->ro2_p, *sk->o_p);
    mnt_red(m1, m1, sk->p, *sk->o_p);

    mpz_mod(m2, sk->c, sk->q);
    mnt_pow_mod(m2, m2, sk->d_q, sk->q, sk->ro2_q, *sk->o_q);
    mnt_red(m2, m2, sk->q, *sk->o_q);

    mpz_sub(m1, m1, m2);
    mpz_mul(m1, sk->i_q, m1);
    mpz_mod(m1, m1, sk->p);
    mpz_mul(m1, m1, sk->q);
    mpz_add(message, m1, m2);
}

// Perform ElGamal encryption using a randomly generated k
void ElGamalEncrypt(ElGamal_public_key *pk, mpz_t c1, mpz_t c2, mpz_t key) {
    mpz_init(pk->k);
    mpz_init(c1);
    mpz_init(c2);

    if(mpz_cmp(key, pk->p) <= 0) {
        fprintf(stderr, "Warning: ElGamal Encrypt key may not be large enough.\n");
    } else {
        mpz_mod(key, key, pk->p);
    }
    mpz_set(pk->k, key);

    // enable for testing
#if defined(DEBUG)
    mpz_set_ui(pk->k, 1);
#endif

    mnt_pow_mod(c1, pk->g, pk->k, pk->p, pk->ro2, *pk->o);
    mnt_pow_mod(c2, pk->h, pk->k, pk->p, pk->ro2, *pk->o);

    mnt_mul(pk->m, pk->m, pk->ro2, pk->p, *pk->o);
    mnt_mul(c2, c2, pk->m, pk->p, *pk->o);

    mnt_red(c1, c1, pk->p, *pk->o);
    mnt_red(c2, c2, pk->p, *pk->o);
}

// Perform ElGamal decryption
void ElGamalDecryption(ElGamal_private_key *sk, mpz_t message) {
    mpz_t qx;
    mpz_init(qx);
    mpz_init(message);

    mpz_sub(qx, sk->q, sk->x);

    mnt_pow_mod(message, sk->c1, qx, sk->p, sk->ro2, *sk->o);

    mnt_mul(sk->c2, sk->c2, sk->ro2, sk->p, *sk->o);
    mnt_mul(message, message, sk->c2, sk->p, *sk->o);

    mnt_red(message, message, sk->p, *sk->o);
}

// Ensure BBS p or q is prime and satisfies = 3 (mod 4)
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

// Generate the next BBS state
int BBS_next(BBS *bbs) {
    mnt_pow_mod(bbs->s, bbs->s, bbs->two, bbs->N, bbs->ro2, *bbs->o);
    mnt_red(bbs->s, bbs->s, bbs->N, *bbs->o);

    return (int) (bbs->s->_mp_d[0] & 1);
}

// Initiate BBS with initial state
int BBS_init(BBS* bbs) {
    mpz_t p, q, N, s;
    mpz_t gcd, ro2;
    mpz_init(p);
    mpz_init(q);
    mpz_init(N);
    mpz_init(s);
    mpz_init(gcd);
    mpz_init(ro2);
    mpz_init(bbs->p);
    mpz_init(bbs->q);
    mpz_init(bbs->N);
    mpz_init(bbs->s);
    mpz_init(bbs->ro2);
    mpz_init(bbs->two);
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

// Perform first stage (RSA encryption)
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

// Perform second stage (RSA decryption)
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

// Perform third stage (ElGamal encryption)
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

// Perform fourth stage (ElGamal decryption)
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

// Parse CLI arguments
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
