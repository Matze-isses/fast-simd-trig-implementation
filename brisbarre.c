#include <float.h>
#include <limits.h>
#include <math.h>

static inline double quantize_pow2(double v, int k) {
  return nearbyint(v * ldexp(1.0, -k)) * ldexp(1.0, k);
}

static inline double mod_pi_over_2(double v) {
  // bring v mod (pi/2) into [0, pi/2)
  double r = fmod(v, M_PI_2);
  if (r < 0) r += M_PI_2;
  return r;
}

double t_hi(int i, double w) {
  // Th(i,w): nearest multiple of 2^-49 to ((2^(8i) * w) mod (pi/2))
  double tmp = mod_pi_over_2(ldexp(w, 8 * i));
  return quantize_pow2(tmp, -49);
}

double t_med(int i, double w) {
  // Tmed(i,w): nearest multiple of 2^-99 to (tmp - Th)
  double tmp = mod_pi_over_2(ldexp(w, 8 * i));
  double thi = t_hi(i, w);
  double rem = tmp - thi;
  return quantize_pow2(rem, -99);
}

double t_lo(int i, double w) {
  // Tlo(i,w): exact-ish residual
  double tmp = mod_pi_over_2(ldexp(w, 8 * i));
  return tmp - t_hi(i, w) - t_med(i, w);
}

/**
 * Range reduction Algorithm for 2^53 < x < 2^63 - 1
 * Based on Brisebarre et al., "A New Range-Reduction Algorithm" (2005).
 * Returns res in [0, 2*pi).
 */
double brisebarre_range_reduction(double x) {
  // For x >= 2^53, the double is already an integer exactly; cast safely.
  // We stay within signed 64-bit as requested: x < 2^63.
  uint64_t I = (uint64_t)x;

  double s_hi  = 0.0;
  double s_med = 0.0;
  double s_lo  = 0.0;

  int i = 7;
  int j = 56;

  while (i >= 0) {
    // Extract the current 7-bit chunk: w in [0,127]
    // (paper uses signed 7-bit via rounding; here we use nonnegative w with separate sign if needed)
    int w = (int)((I >> j) & 0x7FULL);  // 0..127

    if (w != 0) {
      // Sign is positive here (I >= 0 and w is nonnegative chunk)
      s_hi  += t_hi(i, (double)w);
      s_med += t_med(i, (double)w);
      s_lo  += t_lo(i, (double)w);
      // subtract the consumed part
      I -= ((uint64_t)w << j);
    }

    i -= 1;
    j -= 8;
  }

  // Combine the triple (simple, keeps your structure; Fast2sum not strictly necessary for returning y)
  double y = s_hi + s_med + s_lo;

  // Optional “second reduction step” toward tighter interval (paper §2.2):
  // If desired, you could subtract the nearest k*(pi/2) when |s_hi| > pi/4.
  // Keeping it simple and robust: normalize to [0, 2*pi).
  double res = fmod(y, 2.0 * M_PI);
  if (res < 0) res += 2.0 * M_PI;

  return res;
}
