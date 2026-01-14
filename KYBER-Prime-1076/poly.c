#include <stdint.h>
#include "params.h"
#include "poly.h"
#include "ntt.h"
#include "reduce.h"
#include "cbd.h"
#include "symmetric.h"

/*************************************************
* Name:        poly_compress
*
* Description: Compression and subsequent serialization of a polynomial
**************************************************/
void poly_compress_5(uint8_t r[KYBER_POLYCOMPRESSEDBYTES_CV], poly *a) {
    unsigned int i, j;
    uint8_t t[8];
    poly_caddq(a);

    for (i = 0; i < KYBER_N / 8; i++) {
        for (j = 0; j < 8; j++) {
            t[j] = ((((uint32_t)a->coeffs[8*i+j] << 5) + q/2) / q) & 31;
        }

        r[0] = (t[0] >> 0) | (t[1] << 5);
        r[1] = (t[1] >> 3) | (t[2] << 2) | (t[3] << 7);
        r[2] = (t[3] >> 1) | (t[4] << 4);
        r[3] = (t[4] >> 4) | (t[5] << 1) | (t[6] << 6);
        r[4] = (t[6] >> 2) | (t[7] << 3);
        r += 5;
    }
    t[0] = ((((uint32_t)a->coeffs[264] << 5) + q/2) / q) & 31;
    t[1] = ((((uint32_t)a->coeffs[265] << 5) + q/2) / q) & 31;
    t[2] = ((((uint32_t)a->coeffs[266] << 5) + q/2) / q) & 31;
    t[3] = ((((uint32_t)a->coeffs[267] << 5) + q/2) / q) & 31;
    t[4] = ((((uint32_t)a->coeffs[268] << 5) + q/2) / q) & 31;
    r[0] = (t[0] >> 0) | (t[1] << 5);
    r[1] = (t[1] >> 3) | (t[2] << 2) | (t[3] << 7);
    r[2] = (t[3] >> 1) | (t[4] << 4);
    r[3] = (t[4] >> 4);
}

void poly_compress_12(uint8_t r[KYBER_POLYCOMPRESSEDBYTES_CU], poly *a) {
    unsigned int i, j;
    uint16_t t[2];
    poly_caddq(a);

    for (i = 0; i < KYBER_N / 2; i++) {
        for (j = 0; j < 2; j++) {
            t[j] = ((((uint32_t)a->coeffs[2*i+j] << 12) + q/2) / q) & 0xfff;
        }

        r[0] = (t[0] >> 0);
        r[1] = (t[0] >> 8) | (t[1] << 4);
        r[2] = (t[1] >> 4);
        r += 3;
    }
    t[0] = ((((uint32_t)a->coeffs[268] << 12) + q/2) / q) & 0xfff;
    r[0] = (t[0] >> 0);
    r[1] = (t[0] >> 8);
}

/*************************************************
* Name:        poly_decompress
*
* Description: De-serialization and subsequent decompression of a polynomial;
*              approximate inverse of poly_compress
**************************************************/
void poly_decompress_5(poly *r, const uint8_t a[KYBER_POLYCOMPRESSEDBYTES_CV]) {
    unsigned int i,j;
    uint8_t t[8];
    for (i = 0; i < KYBER_N / 8; i++) {
        t[0] = (a[0] >> 0);
        t[1] = (a[0] >> 5) | (a[1] << 3);
        t[2] = (a[1] >> 2);
        t[3] = (a[1] >> 7) | (a[2] << 1);
        t[4] = (a[2] >> 4) | (a[3] << 4);
        t[5] = (a[3] >> 1);
        t[6] = (a[3] >> 6) | (a[4] << 2);
        t[7] = (a[4] >> 3);
        a += 5;

        for (j = 0; j < 8; j++)
            r->coeffs[8 * i + j] = ((uint32_t) (t[j] & 31) * q + 16) >> 5;
    }
    t[0] = (a[0] >> 0);
    t[1] = (a[0] >> 5) | (a[1] << 3);
    t[2] = (a[1] >> 2);
    t[3] = (a[1] >> 7) | (a[2] << 1);
    t[4] = (a[2] >> 4) | (a[3] << 4);
    r->coeffs[264] = ((uint32_t) (t[0] & 31) * q + 16) >> 5;
    r->coeffs[265] = ((uint32_t) (t[1] & 31) * q + 16) >> 5;
    r->coeffs[266] = ((uint32_t) (t[2] & 31) * q + 16) >> 5;
    r->coeffs[267] = ((uint32_t) (t[3] & 31) * q + 16) >> 5;
    r->coeffs[268] = ((uint32_t) (t[4] & 31) * q + 16) >> 5;
}

void poly_decompress_12(poly *r, const uint8_t a[KYBER_POLYCOMPRESSEDBYTES_CU]) {
    unsigned int i, j;
    uint16_t t[2];
    for (i = 0; i < KYBER_N / 2; i++) {
        t[0] = (a[0] >> 0) | ((uint16_t) a[1] << 8);
        t[1] = (a[1] >> 4) | ((uint16_t) a[2] << 4);
        a += 3;

        for (j = 0; j < 2; j++)
            r->coeffs[2 * i + j] = ((uint32_t) (t[j] & 0xFFF) * q + 2048) >> 12;
    }
    t[0] = (a[0] >> 0) | ((uint16_t) a[1] << 8);
    r->coeffs[268] = ((uint32_t) (t[0] & 0xFFF) * q + 2048) >> 12;
}

/*************************************************
* Name:        poly_tobytes
*
* Description: Serialization of a polynomial
**************************************************/
void poly_tobytes(uint8_t r[KYBER_POLYBYTES], poly *a){
    unsigned int i, j;
    uint16_t t[8];
    poly_caddq(a);

    for (i = 0; i < KYBER_N / 8; i++) {
        for (j = 0; j < 8; j++) {
            t[j] = a->coeffs[8 * i + j];
        }

        r[0] = (t[0] >> 0);
        r[1] = (t[0] >> 8) | (t[1] << 5);
        r[2] = (t[1] >> 3);
        r[3] = (t[1] >> 11) | (t[2] << 2);
        r[4] = (t[2] >> 6) | (t[3] << 7);
        r[5] = (t[3] >> 1);
        r[6] = (t[3] >> 9) | (t[4] << 4);
        r[7] = (t[4] >> 4);
        r[8] = (t[4] >> 12) | (t[5] << 1);
        r[9] = (t[5] >> 7) | (t[6] << 6);
        r[10] = (t[6] >> 2);
        r[11] = (t[6] >> 10) | (t[7] << 3);
        r[12] = (t[7] >> 5);
        r += 13;
    }
    t[0] = a->coeffs[264];
    t[1] = a->coeffs[265];
    t[2] = a->coeffs[266];
    t[3] = a->coeffs[267];
    t[4] = a->coeffs[268];

    r[0] = (t[0] >> 0);
    r[1] = (t[0] >> 8) | (t[1] << 5);
    r[2] = (t[1] >> 3);
    r[3] = (t[1] >> 11) | (t[2] << 2);
    r[4] = (t[2] >> 6) | (t[3] << 7);
    r[5] = (t[3] >> 1);
    r[6] = (t[3] >> 9) | (t[4] << 4);
    r[7] = (t[4] >> 4);
    r[8] = (t[4] >> 12);
}

/*************************************************
* Name:        poly_frombytes
*
* Description: De-serialization of a polynomial;
*              inverse of poly_tobytes
**************************************************/
void poly_frombytes(poly *r, const uint8_t a[KYBER_POLYBYTES]) {
    unsigned int i, j;
    uint16_t t[8];
    for (i = 0; i < KYBER_N / 8; i++) {
        t[0] = (a[0] >> 0) | ((uint16_t) a[1] << 8);
        t[1] = (a[1] >> 5) | ((uint16_t) a[2] << 3) | ((uint16_t) a[3] << 11);
        t[2] = (a[3] >> 2) | ((uint16_t) a[4] << 6);
        t[3] = (a[4] >> 7) | ((uint16_t) a[5] << 1) | ((uint16_t) a[6] << 9);
        t[4] = (a[6] >> 4) | ((uint16_t) a[7] << 4) | ((uint16_t) a[8] << 12);
        t[5] = (a[8] >> 1) | ((uint16_t) a[9] << 7);
        t[6] = (a[9] >> 6) | ((uint16_t) a[10] << 2) | ((uint16_t) a[11] << 10);
        t[7] = (a[11] >> 3) | ((uint16_t) a[12] << 5);
        a += 13;

        for (j = 0; j < 8; j++)
            r->coeffs[8 * i + j] = t[j] & 0x1FFF;
    }
    t[0] = (a[0] >> 0) | ((uint16_t) a[1] << 8);
    t[1] = (a[1] >> 5) | ((uint16_t) a[2] << 3) | ((uint16_t) a[3] << 11);
    t[2] = (a[3] >> 2) | ((uint16_t) a[4] << 6);
    t[3] = (a[4] >> 7) | ((uint16_t) a[5] << 1) | ((uint16_t) a[6] << 9);
    t[4] = (a[6] >> 4) | ((uint16_t) a[7] << 4) | ((uint16_t) a[8] << 12);

    r->coeffs[264] = t[0] & 0x1FFF;
    r->coeffs[265] = t[1] & 0x1FFF;
    r->coeffs[266] = t[2] & 0x1FFF;
    r->coeffs[267] = t[3] & 0x1FFF;
    r->coeffs[268] = t[4] & 0x1FFF;
}

/*************************************************
* Name:        polyeta_tobytes
*
* Description: Serialization of a polynomial sampled in cbd;
**************************************************/
void polyeta_tobytes(uint8_t r[KYBER_SECRETBYTES], const poly *a) {
    unsigned int i, j;
    uint8_t t[8];

    for(i = 0; i < KYBER_N/8; i++) {
        for (j = 0; j < 8; j++) {
            t[j] = KYBER_ETA1 - a->coeffs[8*i+j];
        }

        r[0]  = (t[0] >> 0) | (t[1] << 3) | (t[2] << 6);
        r[1]  = (t[2] >> 2) | (t[3] << 1) | (t[4] << 4) | (t[5] << 7);
        r[2]  = (t[5] >> 1) | (t[6] << 2) | (t[7] << 5);
        r += 3;
    }
    t[0] = KYBER_ETA1 - a->coeffs[264];
    t[1] = KYBER_ETA1 - a->coeffs[265];
    t[2] = KYBER_ETA1 - a->coeffs[266];
    t[3] = KYBER_ETA1 - a->coeffs[267];
    t[4] = KYBER_ETA1 - a->coeffs[268];

    r[0]  = (t[0] >> 0) | (t[1] << 3) | (t[2] << 6);
    r[1]  = (t[2] >> 2) | (t[3] << 1) | (t[4] << 4);
  
}

/*************************************************
* Name:        polyeta_frombytes
*
* Description: De-serialization of a polynomial sampled in cbd;
*              inverse of polyeta_tobytes
*************************************************/
void polyeta_frombytes(poly *r, const uint8_t a[KYBER_SECRETBYTES]) {
    unsigned int i,j;
    uint8_t t[8];

    for (i = 0; i < KYBER_N/8; i++) {
        t[0] = (a[0] >> 0) & 7;
        t[1] = (a[0] >> 3) & 7;
        t[2] = ((a[0] >> 6) | (a[1] << 2)) & 7;
        t[3] = (a[1] >> 1) & 7;
        t[4] = (a[1] >> 4) & 7;
        t[5] = ((a[1] >> 7) | (a[2] << 1)) & 7;
        t[6] = (a[2] >> 2) & 7;
        t[7] = (a[2] >> 5) & 7;
        a += 3;

        for (j = 0; j < 8; j++)
            r->coeffs[8 * i + j] = KYBER_ETA1 - t[j];
    }
    t[0] = (a[0] >> 0) & 7;
    t[1] = (a[0] >> 3) & 7;
    t[2] = ((a[0] >> 6) | (a[1] << 2)) & 7;
    t[3] = (a[1] >> 1) & 7;
    t[4] = (a[1] >> 4) & 7;

    r->coeffs[264] = KYBER_ETA1 - t[0];
    r->coeffs[265] = KYBER_ETA1 - t[1];
    r->coeffs[266] = KYBER_ETA1 - t[2];
    r->coeffs[267] = KYBER_ETA1 - t[3];
    r->coeffs[268] = KYBER_ETA1 - t[4];
 
}

/*************************************************
* Name:        poly_frommsg
*
* Description: Convert 32-byte message to polynomial
**************************************************/
void poly_frommsg(poly *r, const uint8_t msg[KYBER_INDCPA_MSGBYTES]) {
    unsigned int i, j;
    int16_t mask;
  
    for (i = 0; i < 32; i++) {
        for (j = 0; j < 8; j++) {
            mask = -(int16_t) ((msg[i] >> j) & 1);
            r->coeffs[8 * i + j] = mask & ((q + 1) / 2);
        }
    }

    for(i = 256; i < 269; i++)
        r->coeffs[i] = 0;
}

/*************************************************
* Name:        poly_tomsg
*
* Description: Convert polynomial to 32-byte message
**************************************************/
void poly_tomsg(uint8_t msg[KYBER_INDCPA_MSGBYTES], const poly *a) {
    unsigned int i, j;
    uint16_t t;

    for (i = 0; i < 32; i++) {
        msg[i] = 0;
        for (j = 0; j < 8; j++) {
            t = a->coeffs[8 * i + j];
            t += ((int16_t) t >> 15) & q; 
            t = (((t << 1) + q / 2) / q) & 1;
            msg[i] |= t << j;
        }
    }
}

/*************************************************
* Name:        poly_getnoise_eta1
*
* Description: Sample a polynomial deterministically from a seed and a nonce,
*              with output polynomial close to centered binomial distribution
*              with parameter KYBER_ETA1
**************************************************/
void poly_getnoise_eta1(poly *r, const uint8_t seed[KYBER_SYMBYTES], uint8_t nonce) {
    uint8_t buf[202];
    prf(buf, sizeof(buf), seed, nonce);
    cbd3(r, buf);
}

/*************************************************
* Name:        poly_getnoise_eta2
*
* Description: Sample a polynomial deterministically from a seed and a nonce,
*              with output polynomial close to centered binomial distribution
*              with parameter KYBER_ETA2
**************************************************/
void poly_getnoise_eta2(poly *r, const uint8_t seed[KYBER_SYMBYTES], uint8_t nonce) {
    uint8_t buf[202];
    prf(buf, sizeof(buf), seed, nonce);
    cbd3(r, buf);
}


/*************************************************
* Name:        poly_ntt
*
* Description: Computes cyclic number-theoretic transform (NTT) of
*              a polynomial in place;
**************************************************/
void poly_ntt(poly_p *r, poly *a) {
    ntt(r->coeffs, a->coeffs);
}

/*************************************************
* Name:        poly_invntt
*
* Description: Computes inverse of cyclic number-theoretic transform (NTT)
*              of a polynomial in place;
**************************************************/
void poly_invntt(poly *r, poly_p *a) {
    invntt(r->coeffs, a->coeffs);
}
/*************************************************
* Name:        poly_pointwisemul
*
* Description: Multiplication of two polynomials in NTT domain
**************************************************/
void polyp_pointwisemul(poly_p *r, const poly_p *a, const poly_p *b)
{
  poly_pointwise(r->coeffs, a->coeffs, b->coeffs);
}
/*************************************************
* Name:        poly_reduce
*
* Description: Applies Barrett reduction to all coefficients of a polynomial
*              for details of the Barrett reduction see comments in reduce.c.
**************************************************/
void poly_reduce(poly *r) {
    unsigned int i;
    for (i = 0; i < KYBER_N; i++)
        r->coeffs[i] = barrett_reduceq(r->coeffs[i]);
}

/*************************************************
* Name:        poly_caddq
*
* Description: For all coefficients of in/out polynomial add q if
*              coefficient is negative.
**************************************************/
void poly_caddq(poly *r) {
  unsigned int i;
  for(i = 0; i < KYBER_N; i++)
    r->coeffs[i] = caddq(r->coeffs[i]);
}

/*************************************************
* Name:        poly_add
*
* Description: Add two polynomials; no modular reduction is performed.
**************************************************/
void poly_add(poly *r, const poly *a, const poly *b) {
    unsigned int i;
    for (i = 0; i < KYBER_N; i++)
        r->coeffs[i] = a->coeffs[i] + b->coeffs[i];
}

void polyp_add(poly_p *r, const poly_p *a, const poly_p *b) {
    unsigned int i;
    for (i = 0; i < N_PRIME; i++)
        r->coeffs[i] = a->coeffs[i] + b->coeffs[i];
}

/*************************************************
* Name:        poly_sub
*
* Description: Subtract two polynomials; no modular reduction is performed.
**************************************************/
void poly_sub(poly *r, const poly *a, const poly *b) {
    unsigned int i;
    for (i = 0; i < KYBER_N; i++)
        r->coeffs[i] = a->coeffs[i] - b->coeffs[i];
}
/*************************************************
* Name:        poly_con
*
* Description: Compression and subsequent serialization of the second ciphertext.
**************************************************/
void poly_con(uint8_t r[KYBER_POLYCOMPRESSEDBYTES_CV], poly *v, const uint8_t m[KYBER_INDCPA_MSGBYTES]) {
    int16_t mcoeffs;

    unsigned int i, j;
    int16_t mask;
    int16_t qround = 3489;

    for (i = 0; i < 32; i++) {
        for (j = 0; j < 8; j++) {
            mask = -(int16_t) ((m[i] >> j) & 1);
            mcoeffs = mask & qround;
            v->coeffs[8 * i + j] += mcoeffs;
            v->coeffs[8 * i + j] = barrett_reduceq(v->coeffs[8 * i + j]);
        }
    }

    uint8_t t[8];
    poly_caddq(v);

    for (i = 0; i < 256 / 8; i++) {
        for (j = 0; j < 8; j++) {
            t[j] = ((((uint32_t)v->coeffs[8*i+j] << 5) + q/2) / q) & 31;
        }

        r[0] = (t[0] >> 0) | (t[1] << 5);
        r[1] = (t[1] >> 3) | (t[2] << 2) | (t[3] << 7);
        r[2] = (t[3] >> 1) | (t[4] << 4);
        r[3] = (t[4] >> 4) | (t[5] << 1) | (t[6] << 6);
        r[4] = (t[6] >> 2) | (t[7] << 3);
        r += 5;
    }

}
/*************************************************
* Name:        poly_rec
*
* Description: The optimized decryption.
**************************************************/
void poly_rec(uint8_t m[KYBER_INDCPA_MSGBYTES], const uint8_t c2[KYBER_POLYCOMPRESSEDBYTES_CV], const poly *su) {
    int16_t c2_coeffs[KYBER_N];
    int64_t tmp;

    unsigned int i,j;
    uint8_t t[8];
    for (i = 0; i < 256 / 8; i++) {
        t[0] = (c2[0] >> 0);
        t[1] = (c2[0] >> 5) | (c2[1] << 3);
        t[2] = (c2[1] >> 2);
        t[3] = (c2[1] >> 7) | (c2[2] << 1);
        t[4] = (c2[2] >> 4) | (c2[3] << 4);
        t[5] = (c2[3] >> 1);
        t[6] = (c2[3] >> 6) | (c2[4] << 2);
        t[7] = (c2[4] >> 3);
        c2 += 5;

        for (j = 0; j < 8; j++)
            c2_coeffs[8 * i + j] = t[j] & 31;
    }

    for (i = 0; i < 32; i++) {
        m[i] = 0;
        for (j = 0; j < 8; j++) {
              tmp = ((int32_t)c2_coeffs[8 * i + j] + 8) * q - ((int32_t)su->coeffs[8 * i + j] << 5);
              tmp = (tmp * 9617) >> 30;
              //tmp = tmp / 111632 - ((uint64_t)tmp >> 63);
              tmp = tmp & 1;
              m[i] |= tmp << j;
        }
    }
}







