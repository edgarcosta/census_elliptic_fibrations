/* Match models under all 6144 coordinate changes preserving a2 = 0.
 * A key packs degree 2, 6, 8, 12 forms as a1<<29 | a3<<22 | a4<<13 | a6.
 * Bit j is the coefficient of s^j t^(d-j); products are carry-less.
 * See docs/native-full-f2.md for the action and file arguments.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static uint32_t clmul(uint32_t a, uint32_t b) {
    uint32_t r = 0;
    while (b) {
        if (b & 1) r ^= a;
        a <<= 1;
        b >>= 1;
    }
    return r;
}

/* GL(2,F_2): M = [[a,b],[c,d]], f(s,t) -> f(a s + b t, c s + d t) */
static const int MATS[6][4] = {
    {1, 0, 0, 1}, {0, 1, 1, 0}, {1, 1, 0, 1},
    {1, 0, 1, 1}, {0, 1, 1, 1}, {1, 1, 1, 0}};
static uint32_t *TAB[6][13];

static void build_tables(void) {
    int degs[4] = {2, 6, 8, 12};
    for (int m = 0; m < 6; m++) {
        uint32_t ls = (uint32_t)(MATS[m][0] << 1 | MATS[m][1]); /* a s + b t */
        uint32_t lt = (uint32_t)(MATS[m][2] << 1 | MATS[m][3]); /* c s + d t */
        for (int k = 0; k < 4; k++) {
            int d = degs[k];
            uint32_t img[13];
            for (int j = 0; j <= d; j++) { /* s^j t^(d-j) */
                uint32_t v = 1;
                for (int e = 0; e < j; e++) v = clmul(v, ls);
                for (int e = 0; e < d - j; e++) v = clmul(v, lt);
                img[j] = v;
            }
            uint32_t n = 1u << (d + 1);
            TAB[m][d] = malloc(n * sizeof(uint32_t));
            for (uint32_t f = 0; f < n; f++) {
                uint32_t v = 0;
                for (int j = 0; j <= d; j++)
                    if (f >> j & 1) v ^= img[j];
                TAB[m][d][f] = v;
            }
        }
    }
}

#define HBITS 23
#define HSIZE (1u << HBITS)
static uint32_t *hkey; /* key + 1, 0 = empty; keys < 2^32 - 1 */
static uint8_t *hcell;
static uint32_t *hidx;

static uint32_t hpos(uint32_t key) {
    uint32_t h = (uint32_t)(((uint64_t)key * 0x9E3779B97F4A7C15ull) >> (64 - HBITS));
    while (hkey[h] && hkey[h] != key + 1) h = (h + 1) & (HSIZE - 1);
    return h;
}

int main(int argc, char **argv) {
    if (argc != 5) {
        fprintf(stderr, "usage\n");
        return 2;
    }
    build_tables();
    hkey = calloc(HSIZE, sizeof(uint32_t));
    hcell = calloc(HSIZE, 1);
    hidx = calloc(HSIZE, sizeof(uint32_t));
    static uint8_t pairset[1 << 10];
    FILE *f = fopen(argv[1], "r");
    char line[256];
    if (!fgets(line, sizeof line, f)) return 1;
    uint32_t nsm = 0, nall = 0;
    unsigned long long key;
    int cell, prefix;
    while (fscanf(f, "%llu,%d,%d", &key, &cell, &prefix) == 3) {
        uint32_t h = hpos((uint32_t)key);
        if (hkey[h]) {
            fprintf(stderr, "duplicate representative key %llu\n", key);
            return 1;
        }
        hkey[h] = (uint32_t)key + 1;
        hcell[h] = (uint8_t)cell;
        hidx[h] = cell == 0 ? nsm++ : 0xFFFFFFFFu;
        pairset[key >> 22] = 1;
        nall++;
    }
    fclose(f);
    fprintf(stderr, "representatives %u, smooth minimal %u\n", nall, nsm);
    uint32_t *hits = calloc(nsm, sizeof(uint32_t));

    FILE *in = fopen(argv[2], "r");
    FILE *out = fopen(argv[3], "w");
    fprintf(out, "row,published_key,n_smooth_min,native_key,M,p,w,other_cells\n");
    unsigned long long uk;
    long row = 0;
    long nmatched1 = 0, nzero = 0, nmulti = 0;
    while (fscanf(in, "%llu", &uk) == 1) {
        row++;
        uint32_t a1 = (uint32_t)(uk >> 29) & 7, a3 = (uint32_t)(uk >> 22) & 127;
        uint32_t a4 = (uint32_t)(uk >> 13) & 511, a6 = (uint32_t)uk & 8191;
        uint32_t found[64];
        int nf = 0;
        uint32_t fm = 0, fp = 0, fw = 0;
        unsigned othermask = 0;
        for (uint32_t p = 0; p < 8; p++) {
            uint32_t r = clmul(p, p) ^ clmul(a1, p);
            uint32_t a3n = a3 ^ clmul(r, a1);
            uint32_t r2 = clmul(r, r), r3 = clmul(r2, r);
            uint32_t a4b = a4 ^ clmul(p, a3) ^ r2 ^ clmul(clmul(r, p), a1);
            uint32_t a6b = a6 ^ clmul(r, a4) ^ r3;
            for (int m = 0; m < 6; m++) {
                uint32_t a1m = TAB[m][2][a1], a3m = TAB[m][6][a3n];
                if (!pairset[a1m << 7 | a3m]) continue;
                for (uint32_t w = 0; w < 128; w++) {
                    uint32_t a4n = a4b ^ clmul(w, a1);
                    uint32_t a6n = a6b ^ clmul(w, a3) ^ clmul(w, w) ^ clmul(clmul(r, w), a1);
                    uint32_t k = a1m << 29 | a3m << 22 | TAB[m][8][a4n] << 13 | TAB[m][12][a6n];
                    uint32_t h = hpos(k);
                    if (!hkey[h]) continue;
                    if (hcell[h] != 0) {
                        othermask |= 1u << hcell[h];
                        continue;
                    }
                    int seen = 0;
                    for (int i = 0; i < nf; i++)
                        if (found[i] == k) seen = 1;
                    if (!seen && nf < 64) {
                        found[nf++] = k;
                        if (nf == 1) fm = (uint32_t)m, fp = p, fw = w;
                    }
                }
            }
        }
        if (nf == 1) {
            nmatched1++;
            hits[hidx[hpos(found[0])]]++;
        } else if (nf == 0)
            nzero++;
        else
            nmulti++;
        fprintf(out, "%ld,%llu,%d,%u,%u,%u,%u,", row, uk, nf, nf ? found[0] : 0, fm, fp, fw);
        int first = 1;
        for (int c = 1; c < 6; c++)
            if (othermask >> c & 1) {
                fprintf(out, "%s%d", first ? "" : ";", c);
                first = 0;
            }
        fprintf(out, "\n");
    }
    fclose(out);
    FILE *hb = fopen(argv[4], "wb");
    fwrite(hits, sizeof(uint32_t), nsm, hb);
    fclose(hb);
    fprintf(stderr, "rows %ld: one smooth minimal key %ld, none %ld, several %ld\n", row,
            nmatched1, nzero, nmulti);
    return 0;
}
