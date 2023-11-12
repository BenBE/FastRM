#include "reedmuller.h"

// Based on https://stackoverflow.com/a/45969046
static size_t binomial(size_t n, size_t k)
{
    if (k == 0) {
        return 1;
    }

    if (k > n) {
        return 0;
    }

    if (k > n / 2) {
        return binomial(n, n - k);
    }

    size_t result = 1;

    for (size_t i = 0; i < k; i++) {
        result = (result * (n - i)) / (i + 1);
    }

    return result;
}

// Adapted from https://stackoverflow.com/a/45969046
static size_t binomial_sum(size_t n, size_t k)
{
    if (k == 0) {
        return 1;
    }

    if (k > n) {
        return n * n;
    }

    size_t sum = 1;
    size_t bc = 1;

    for (size_t i = 0; i < k; i++) {
        bc = (bc * (n - i)) / (i + 1);
        sum += bc;
    }

    return sum;
}

bool reedmuller_isvalid(uint8_t r, uint8_t m)
{
    if (m >= __CHAR_BIT__ * sizeof(size_t)) {
        return false;
    }

    if (r == 0 || r >= m) {
        return false;
    }

    return true;
}

size_t reedmuller_rawsize(uint8_t r, uint8_t m)
{
    if (!reedmuller_isvalid(r, m)) {
        return 0;
    }

    return binomial_sum(m, r);
}

size_t reedmuller_blocksize(uint8_t r, uint8_t m)
{
    if (!reedmuller_isvalid(r, m)) {
        return 0;
    }

    return 1ull << m;
}

size_t reedmuller_maxerror(uint8_t r, uint8_t m)
{
    if (!reedmuller_isvalid(r, m)) {
        return 0;
    }

    return (1ull << (m - r - 1)) - 1ull;
}

static uint8_t reedmuller_getvalue_bit(size_t value, size_t idx)
{
    if (idx >= __CHAR_BIT__ * sizeof(value)) {
        return 0;
    }

    return !!(value & (1 << idx));
}

static uint8_t reedmuller_getbuffer_bit(const uint8_t* buf, size_t size, size_t idx)
{
    if (idx >= 8 * size) {
        return 0;
    }

    return reedmuller_getvalue_bit(buf[idx / 8], idx % 8);
}

static void reedmuller_setbuffer_bit(uint8_t* buf, size_t size, size_t idx, uint8_t bit)
{
    if (idx >= 8 * size) {
        return;
    }

    uint8_t mask = (uint8_t)(1 << (idx % 8));

    buf[idx / 8] = (buf[idx / 8] & ~mask) | (!!bit * mask);
}

// A function to generate a combination of bits given an index and a maximum number of bits to set
static size_t reedmuller_bitmask(uint8_t r, uint8_t m, size_t idx) {
    if (!idx) {
        return 0;
    }

    size_t bit = 0;

    // For each number of bits
    for (size_t count = 0; count <= r; count++) {
        // Too many bits
        if (count > r) {
            return 0;
        }

        // Add another bit
        size_t perms = binomial(m, count);
        if (idx >= perms) {
            bit <<= 1;
            bit |= 1;
            idx -= perms;
            continue;
        }

        for (; idx > 0; idx--) {
            size_t t = bit | (bit - 1); // t gets v's least significant 0 bits set to 1
            bit = (t + 1) | (((~t & -~t) - 1) >> (__builtin_ctz(bit) + 1));
        }

        return bit;
    }

    return 0;
}

static uint8_t reedmuller_parity_eval(uint8_t r, uint8_t m, size_t i, const uint8_t* buf, size_t size, size_t bitoffs)
{
    if (!reedmuller_isvalid(r, m)) {
        return 0;
    }

    const size_t rsize = reedmuller_rawsize(r, m);

    uint8_t outbit = 0;

    for (size_t idx = 0; idx < rsize; idx++) {
        size_t mask = reedmuller_bitmask(r, m, idx);
        uint8_t bit = reedmuller_getbuffer_bit(buf, size, bitoffs + idx);

        outbit ^= bit * ((i & mask) == mask);
    }

    return outbit;
}

static bool reedmuller_encode_internal(uint8_t r, uint8_t m, const uint8_t* ibuf, size_t isize, size_t ibit, uint8_t* obuf, size_t osize, size_t obit)
{
    if (!reedmuller_isvalid(r, m)) {
        return false;
    }

    const size_t bsize = reedmuller_blocksize(r, m);

    if (obit > osize || osize - obit < bsize) {
        return false;
    }

    for (size_t idx = 0; idx < bsize; idx++) {
        reedmuller_setbuffer_bit(
            obuf, osize, obit + idx,
            reedmuller_parity_eval(r, m, idx, ibuf, isize, ibit)
        );
    }

    return true;
}

bool reedmuller_encode(uint8_t r, uint8_t m, const uint8_t* ibuf, size_t isize, uint8_t* obuf, size_t osize)
{
    if (!reedmuller_isvalid(r, m)) {
        return false;
    }

    isize *= __CHAR_BIT__;
    osize *= __CHAR_BIT__;

    const size_t rsize = reedmuller_rawsize(r, m);
    const size_t bsize = reedmuller_blocksize(r, m);

    bool result = true;

    size_t ibit = 0;
    size_t obit = 0;

    while (result && ibit < isize && obit + bsize <= osize) {
        result &= reedmuller_encode_internal(r, m, ibuf, isize, ibit, obuf, osize, obit);

        ibit += rsize;
        obit += bsize;
    }

    return result;
}

static size_t reedmuller_maskapply(size_t mask, size_t index)
{
    size_t result = 0;

    while (mask) {
        size_t lsb = mask ^ (mask & (mask - 1));
        mask ^= lsb;

        lsb *= index & 1;
        index >>= 1;

        result |= lsb;
    }

    return result;
}

static bool reedmuller_decode_internal(uint8_t r, uint8_t m, uint8_t* ibuf, size_t isize, size_t ibit, uint8_t* obuf, size_t osize, size_t obit)
{
    if (!reedmuller_isvalid(r, m)) {
        return false;
    }

    const size_t bsize = reedmuller_blocksize(r, m);
    const size_t rsize = reedmuller_rawsize(r, m);

    const size_t max_err = reedmuller_maxerror(r, m);

    if (obit > osize || osize - obit < rsize) {
        return false;
    }

    for (size_t terms = 0, tidx = rsize - 1; terms < rsize; terms++, tidx--) {
        size_t tmask = reedmuller_bitmask(r, m, tidx);
        size_t tbit = __builtin_popcount(tmask);
        size_t tcount = 1 << tbit;

        size_t epmask = ((1 << m) - 1) ^ tmask;
        size_t epbit = __builtin_popcount(epmask);
        size_t epcount = 1 << epbit;

        size_t eval[2] = {0, 0};

        // For all evaluation points
        for (size_t epidx = 0; epidx < epcount; epidx++) {
            size_t eppattern = reedmuller_maskapply(epmask, epidx);
            size_t epresult = 0;
            for (size_t tbidx = 0; tbidx < tcount; tbidx++) {
                size_t tbpattern = reedmuller_maskapply(tmask, tbidx) | eppattern;
                epresult ^= reedmuller_getbuffer_bit(ibuf, isize, ibit + tbpattern);
            }
            epresult &= 1;
            eval[epresult]++;
        }

        if (eval[0] > max_err && eval[1] > max_err) {
            return false;
        }

        reedmuller_setbuffer_bit(obuf, osize, obit + tidx, eval[1] > eval[0]);

        if (eval[1] > eval[0]) {
            for (size_t iidx = 0; iidx < bsize; iidx++) {
                if ((iidx & tmask) != tmask) {
                    continue;
                }

                reedmuller_setbuffer_bit(
                    ibuf, isize, ibit + iidx,
                    !reedmuller_getbuffer_bit(ibuf, isize, ibit + iidx)
                );
            }
        }
    }

    return true;
}

bool reedmuller_decode(uint8_t r, uint8_t m, uint8_t* ibuf, size_t isize, uint8_t* obuf, size_t osize)
{
    if (!reedmuller_isvalid(r, m)) {
        return false;
    }

    isize *= __CHAR_BIT__;
    osize *= __CHAR_BIT__;

    const size_t rsize = reedmuller_rawsize(r, m);
    const size_t bsize = reedmuller_blocksize(r, m);

    bool result = true;

    size_t ibit = 0;
    size_t obit = 0;

    while (result && ibit + bsize <= isize && obit + rsize <= osize) {
        result &= reedmuller_decode_internal(r, m, ibuf, isize, ibit, obuf, osize, obit);

        if (isize < bsize / 8) {
            break;
        }

        ibit += bsize;
        obit += rsize;
    }

    return result;
}
