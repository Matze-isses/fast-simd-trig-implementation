#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>   // add this include

#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <inttypes.h>

#include "test_interface.h"
#include "trig_simd.h"

static inline uint64_t dbl_to_u64(double x) {
    uint64_t u;
    memcpy(&u, &x, sizeof u);
    return u;
}

static inline double u64_to_dbl(uint64_t u) {
    double x;
    memcpy(&x, &u, sizeof x);
    return x;
}

/* Change the raw IEEE-754 bit pattern by +/- delta (LSB steps). */
static inline double tweak_lsb(double x, int delta) {
    uint64_t u = dbl_to_u64(x);

    // Optional safety: skip NaN/Inf (all exponent bits = 1)
    if ((u & 0x7ff0000000000000ULL) == 0x7ff0000000000000ULL) {
        return x;
    }

    u = (uint64_t)((int64_t)u + (int64_t)delta);
    return u64_to_dbl(u);
}

static void print_coeff_hex(const double *c, int n) {
    for (int i = 0; i < n; i++) {
        // C99 hex-float literal, pasteable back into C
        printf("0x coeff[%2d] = %a\n", i, c[i]);
    }
}


int main(void)
{
    int n = 256;

    /* Those values are 3 ULP off */
    double test_data[256] = {
        0x1.0ad052644c126p+0,
        0x1.12aebe1b349c4p+0,
        0x1.0a8116fd9ce68p+0,
        0x1.1182d643f0e6cp+0,
        0x1.1149449f820a5p+0,
        0x1.0ea30e60bd412p+0,
        0x1.06b2187de96f1p+0,
        0x1.1324948140f3ap+0,
        0x1.07f48880ba7d3p+0,
        0x1.1192c0b56b4fbp+0,
        0x1.0ef3a0c4b8556p+0,
        0x1.0ce126980c900p+0,
        0x1.1364341f6ab17p+0,
        0x1.17a9984777d73p+0,
        0x1.08ab1305257d7p+0,
        0x1.0e1f663d5df00p+0,
        0x1.10cfc4ea643d9p+0,
        0x1.116d626051330p+0,
        0x1.0c2d4e529313bp+0,
        0x1.18a7fd7bf8e0bp+0,
        0x1.0cd7d2ed207ecp+0,
        0x1.0c92d5e2b9798p+0,
        0x1.1440bcd3ed930p+0,
        0x1.12d41beb429a5p+0,
        0x1.14359a33adf5cp+0,
        0x1.131710aaef904p+0,
        0x1.0c4bc2351a4bep+0,
        0x1.13dcb2919ea9ap+0,
        0x1.1339494c5963bp+0,
        0x1.0e8f1a0cf44f2p+0,
        0x1.0c5c3b0a396cfp+0,
        0x1.0eeafe90840d2p+0,
        0x1.100206759e220p+0,
        0x1.0ddf0d1f81ee6p+0,
        0x1.053d14b645d4bp+0,
        0x1.13fc9e0d06c76p+0,
        0x1.04ec0da455f76p+0,
        0x1.13c6222931488p+0,
        0x1.11a45eeb8caecp+0,
        0x1.10cfadead4d75p+0,
        0x1.127ce4a01d33ep+0,
        0x1.0791921b95d33p+0,
        0x1.0c5bbb6e239a5p+0,
        0x1.13e57e65f9dd7p+0,
        0x1.1293aa38fe888p+0,
        0x1.12a64154dc27fp+0,
        0x1.0d942355ff5a9p+0,
        0x1.12d5be71f7a4ap+0,
        0x1.12dec11d10118p+0,
        0x1.180fbe720a584p+0,
        0x1.0fee566b5e80fp+0,
        0x1.1460c75c300acp+0,
        0x1.1449bac2ef7cdp+0,
        0x1.13839781f6968p+0,
        0x1.08988cc10e76ap+0,
        0x1.13bd0a669e430p+0,
        0x1.122c30db94e06p+0,
        0x1.0b6a54bf933e4p+0,
        0x1.0513547b64830p+0,
        0x1.13396bd508849p+0,
        0x1.11b0b49a8e29bp+0,
        0x1.0f08977131116p+0,
        0x1.05868976c9722p+0,
        0x1.13accfbcd54bfp+0,
        0x1.0d59ce4c8e90ap+0,
        0x1.11261fa630966p+0,
        0x1.06857ed8e9512p+0,
        0x1.1390a680bcb86p+0,
        0x1.0eafc241f4b52p+0,
        0x1.145f4ce2270d5p+0,
        0x1.0af5537b1e460p+0,
        0x1.0d7139ed320ecp+0,
        0x1.1111ca3175eb1p+0,
        0x1.145ad53445bf3p+0,
        0x1.0c6e0b2e5f994p+0,
        0x1.0c9820e9c661ep+0,
        0x1.0af0d545e5b63p+0,
        0x1.0be5de7a2c475p+0,
        0x1.0e5fbb238cbcap+0,
        0x1.139c5c55978e2p+0,
        0x1.0d8fbd11bcc52p+0,
        0x1.113cde0ad7845p+0,
        0x1.0b46e0b6e090ep+0,
        0x1.13fbb35e04a33p+0,
        0x1.0b0562255c7f5p+0,
        0x1.12b11d1d6d301p+0,
        0x1.020b2345c3f52p+0,
        0x1.10708fad3313fp+0,
        0x1.f6c8a7b67a0abp-1,
        0x1.1482ac3751ea2p+0,
        0x1.1470f864d7f56p+0,
        0x1.106f9f78f2daap+0,
        0x1.09ca49ef18c6dp+0,
        0x1.11a77a8bdb478p+0,
        0x1.08eaa27605221p+0,
        0x1.142deb6a8442ap+0,
        0x1.10d3b51cc7553p+0,
        0x1.1063e71b02956p+0,
        0x1.106542efd8174p+0,
        0x1.12d824183dd58p+0,
        0x1.1097c71e5870bp+0,
        0x1.138a9e157d1acp+0,
        0x1.0fe85ac5044afp+0,
        0x1.0fbec2bb46484p+0,
        0x1.0accdf24ae3a6p+0,
        0x1.0e8f4d735a1d9p+0,
        0x1.00103861b8063p+0,
        0x1.197e835882459p+0,
        0x1.133b9d219e19ap+0,
        0x1.194f05dee4f58p+0,
        0x1.135d1e4ddf36ep+0,
        0x1.12697b02488bep+0,
        0x1.11f70ab7aa174p+0,
        0x1.12c012c56eb92p+0,
        0x1.0ac18b47c05c5p+0,
        0x1.08fbbe65565c1p+0,
        0x1.072e85f8f424bp+0,
        0x1.0b6669a7a9b73p+0,
        0x1.12cec9c3c6f79p+0,
        0x1.111ae20e31a02p+0,
        0x1.0dc604563a7e8p+0,
        0x1.0a59347c4e571p+0,
        0x1.0b94ea092a117p+0,
        0x1.13059f2797859p+0,
        0x1.0c4003fc30cf6p+0,
        0x1.1362c26924ef0p+0,
        0x1.1496679696c54p+0,
        0x1.12f76f72b7924p+0,
        0x1.0a1f4c539e271p+0,
        0x1.10929169eaffep+0,
        0x1.0700ec412e22dp+0,
        0x1.0cfcb8b748569p+0,
        0x1.0f9ab8f956a66p+0,
        0x1.135523b69f218p+0,
        0x1.137333ba9fdd3p+0,
        0x1.1489e8044cd39p+0,
        0x1.09b6f1c5b9b91p+0,
        0x1.1119003a3d670p+0,
        0x1.10fa88794544bp+0,
        0x1.14ae8ac5ca15cp+0,
        0x1.0f42c28ccf2f9p+0,
        0x1.1078f7a9d676dp+0,
        0x1.0bff31fa6aa0dp+0,
        0x1.123b7fe9aef74p+0,
        0x1.12dc9e8bb71d2p+0,
        0x1.140ce4638ddd4p+0,
        0x1.11aaa557613cep+0,
        0x1.0faaa77e70267p+0,
        0x1.149cc506ba928p+0,
        0x1.1476d44cc6e03p+0,
        0x1.10f32b0700126p+0,
        0x1.12784bb9fca22p+0,
        0x1.008958a20dd11p+0,
        0x1.11cf8f0f73e59p+0,
        0x1.115dd2d21abf9p+0,
        0x1.13f6314e355f0p+0,
        0x1.0e29e7c9e58f0p+0,
        0x1.135f5523854d9p+0,
        0x1.13b40b77566cbp+0,
        0x1.1279f4dac341ap+0,
        0x1.0428176ebde88p+0,
        0x1.0f4587bc5cf2ep+0,
        0x1.12b857a1d8f84p+0,
        0x1.1481df0fd8b6bp+0,
        0x1.12b871a45665cp+0,
        0x1.05ad40111c965p+0,
        0x1.0e28c94b3a0e3p+0,
        0x1.0c778e7c011a9p+0,
        0x1.048236d3b4e80p+0,
        0x1.0a54fb24b8a15p+0,
        0x1.1314dd2a24d86p+0,
        0x1.1342b5c085e86p+0,
        0x1.13eeb117be9f0p+0,
        0x1.19980fa56330bp+0,
        0x1.0e5868133d41ap+0,
        0x1.12eed6e483244p+0,
        0x1.1213502bf4ff9p+0,
        0x1.10155cc7891d2p+0,
        0x1.119b266f82d5fp+0,
        0x1.104e1aa78cfb3p+0,
        0x1.13f394d8a354bp+0,
        0x1.103d60bf5b436p+0,
        0x1.06ce4694cc1bfp+0,
        0x1.0f9ff4945a761p+0,
        0x1.083c0bd06bd43p+0,
        0x1.198d83548ab5ep+0,
        0x1.0ae6844acc871p+0,
        0x1.0f96d46decc6bp+0,
        0x1.08cce5b572af2p+0,
        0x1.1802321b4b92ep+0,
        0x1.130a3928e89f2p+0,
        0x1.10d7b3aa49f1ep+0,
        0x1.0cc2c42df684ap+0,
        0x1.1255fb7b26f19p+0,
        0x1.0b45e9c5ce552p+0,
        0x1.0f3ed31f8b07cp+0,
        0x1.0ecc3778d43edp+0,
        0x1.0f09173f540c3p+0,
        0x1.0b3141c1d936bp+0,
        0x1.1271bca3b14fep+0,
        0x1.0f036fc6210e2p+0,
        0x1.13708bb45f4ccp+0,
        0x1.127187babde11p+0,
        0x1.0e7e267200cf7p+0,
        0x1.0f985f4ed7fd2p+0,
        0x1.16e0a69a1c6fbp+0,
        0x1.10145069fad6fp+0,
        0x1.134ddd4c2d26dp+0,
        0x1.1345e3e4551dap+0,
        0x1.16d012b7c3c47p+0,
        0x1.0fda6700ebf35p+0,
        0x1.0dc357c022beep+0,
        0x1.1306edc605748p+0,
        0x1.13e6fd04501fcp+0,
        0x1.14a88876bb641p+0,
        0x1.0741950c12e86p+0,
        0x1.1371e68b3e747p+0,
        0x1.144cb32e82664p+0,
        0x1.145cd732ed38bp+0,
        0x1.12f529dbc98e0p+0,
        0x1.0c20bb27629e1p+0,
        0x1.1238fee925424p+0,
        0x1.0ed355ab0516ep+0,
        0x1.114bfc5e43cfdp+0,
        0x1.13f088d192176p+0,
        0x1.1078f21a0fd3fp+0,
        0x1.12e7c35b879c1p+0,
        0x1.10fcaaa827e28p+0,
        0x1.12684728314b7p+0,
        0x1.12641f8f6d1b9p+0,
        0x1.1140c9b40feeep+0,
        0x1.090634cce8a49p+0,
        0x1.106a6b92e655cp+0,
        0x1.13c3b1e4b1805p+0,
        0x1.0ec4d9ceef40bp+0,
        0x1.14073a55ad3b8p+0,
        0x1.0eaa2a5202cd3p+0,
        0x1.0fb79a1ea6784p+0,
        0x1.0f721c347cb68p+0,
        0x1.11cb7fa9242d6p+0,
        0x1.129563290d756p+0,
        0x1.1082c63d36123p+0,
        0x1.0ebce6e30c299p+0,
        0x1.12c94a967de78p+0,
        0x1.0ad515a9ccb3fp+0,
        0x1.14c31fedee709p+0,
        0x1.12da7ff14438dp+0,
        0x1.0f660cb341beap+0,
        0x1.11acff04fc894p+0,
        0x1.10001df8c6b90p+0,
        0x1.fcf8b35940fedp-1,
        0x1.13cf314ec876ap+0,
        0x1.1001de5c80f29p+0,
        0x1.0e99dde94f364p+0,
        0x1.08599730c7105p+0,
        0x1.14a3d4a8b46eap+0
    };

    double taylor_poly[13] = {
        0.3333333333333333,
        0.13333333333333333,
        0.05396825396825397,
        0.021869488536155203,
        0.008863235529902197,
        0.003592128036572481,
        0.0014558343870513183,
        0.000590027440945586,
        0.00023912911424355248,
        9.691537956929451e-05,
        3.927832388331683e-05,
        1.5918905069328964e-05,
        6.451689215655431e-06
    };

    const int range = 0; // {-1,0,+1}
    const uint64_t base = (uint64_t)(2 * range + 1);

    uint64_t total = 1;
    for (int i = 0; i < 13; i++) total *= base;

    // global best (shared)
    double best_coeff[13];
    int best_shifts[13];
    int64_t best_total_ulp = INT64_MAX;
    uint64_t best_idx = 0;

    // progress printing cadence
    const uint64_t print_step = 10000ULL;

    #pragma omp parallel
    {
        // thread-local working buffers
        double coeff[13];
        int cur_shifts[13];

        double *res = (double*)malloc((size_t)n * sizeof(double));
        int64_t *ulp_error = (int64_t*)malloc((size_t)n * sizeof(int64_t));
        if (!res || !ulp_error) {
            fprintf(stderr, "malloc failed in thread %d\n", omp_get_thread_num());
            abort();
        }

        // thread-local best
        double best_coeff_t[13];
        int best_shifts_t[13];
        int64_t best_total_ulp_t = INT64_MAX;
        uint64_t best_idx_t = 0;

        #pragma omp for schedule(dynamic, 1024)
        for (uint64_t idx = 0; idx < total; idx++) {
            uint64_t x = idx;

            for (int i = 0; i < 13; i++) {
                int digit = (int)(x % base);
                int delta = digit - range;
                coeff[i] = tweak_lsb(taylor_poly[i], delta);
                x /= base;
                cur_shifts[i] = delta;
            }

            vfast_tan(test_data, res, coeff, (size_t)n);
            compare_results_tan_ulp_err_signed(test_data, res, ulp_error, (size_t)n);

            int64_t total_ulp = 0;
            for (int i = 0; i < n; i++) {
                int64_t e = ulp_error[i];
                total_ulp += (e < 0) ? -e : e;
            }

            // update thread-local best (no locking)
            if (total_ulp < best_total_ulp_t) {
                best_total_ulp_t = total_ulp;
                best_idx_t = idx;
                memcpy(best_coeff_t, coeff, sizeof(best_coeff_t));
                memcpy(best_shifts_t, cur_shifts, sizeof(best_shifts_t));
            }

            // progress: only thread 0 prints (avoids interleaving)
            if (omp_get_thread_num() == 0 && (idx % print_step == 0)) {
                double progress = 100.0 * (double)idx / (double)total;
                printf("progress=%6.2f%% idx=%" PRIu64 "\n", progress, idx);
            }
        }

        // merge thread best into global best once per thread
        #pragma omp critical
        {
            if (best_total_ulp_t < best_total_ulp) {
                best_total_ulp = best_total_ulp_t;
                best_idx = best_idx_t;
                memcpy(best_coeff, best_coeff_t, sizeof(best_coeff_t));
                memcpy(best_shifts, best_shifts_t, sizeof(best_shifts_t));

                printf("NEW BEST: idx=%" PRIu64 " total_ulp=%" PRId64 "\n",
                       best_idx, best_total_ulp);
            }
        }

        free(res);
        free(ulp_error);
    } // end parallel

    printf("\n=== BEST OVER ALL ===\n");
    printf("best_idx=%" PRIu64 " best_total_ulp=%" PRId64 "\n", best_idx, best_total_ulp);
    printf("best coeffs as C99 hex-floats:\n\n");
    for (int i = 0; i < 13; i++) {
        printf("%a,\n", best_coeff[i]);
    }

    printf("\nbest Shifts:\n\n");
    for (int i = 0; i < 13; i++) {
        printf("%d: %d,\n", i, best_shifts[i]);
    }

    return 0;
}
