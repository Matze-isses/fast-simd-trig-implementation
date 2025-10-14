#include "clock_utils.h"
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "bit_printing.h"

#define MAX_EXACT_INT (9007199254740992)
#define LONG_RANGE_END (57952155664616982739.0)

const double END_FRIENDLY_RANGE = M_PI * 2.0 * INT_MAX;
const double RANGE_MAX = M_PI * 2.0;
const double TWO_POW_NEG_49 = pow(2, -49);
const double TWO_POW_NEG_99 = pow(2, -99);


typedef struct {
  double value;
  const char *name;
  int quadrant;
} UglyCaseTriplet;


// ################################## brisbarre #########################################################################


void fast2sum(double a, double b, double *s, double *r) {
    *s = a + b;
    double z = *s - a;
    *r = b - z;
}

static inline double quantize_pow2(double v, int k) {
  return nearbyint(v * ldexp(1.0, -k)) * ldexp(1.0, k);
}

static inline double mod_pi(double v) {
  double r = fmod(v, RANGE_MAX);
  if (r < 0) r += RANGE_MAX;
  return r;
}

void calc_t(int i, int8_t w, double *t_hi, double *t_med, double *t_lo) {
  // printf("w:  %.d\n", w);
  // printf("uw: %.d\n", abs(w));

  double base = fmod(pow(2, 8 * i) * w, RANGE_MAX);
  *t_hi = round(base / TWO_POW_NEG_49) * TWO_POW_NEG_49;

  base -= *t_hi;
  *t_med = round(base / TWO_POW_NEG_99) * TWO_POW_NEG_99;

  *t_lo = base - *t_med;
}

int reduce_s_hi(double x) {
  int k = (int)round(x / RANGE_MAX);
  return k;
}

/**
 * Range reduction Algorithm for 2^53 < x < 2^63 - 1
 * Based on Brisebarre et al., "A New Range-Reduction Algorithm" (2005).
 * Returns res in [0, 2*pi).
 */
double brisebarre_range_reduction(double x) {
  // For x >= 2^53, the double is already an integer exactly; cast safely.
  // We stay within signed 64-bit as requested: x < 2^63.
  long I = (long) llround(x);
  double p = x - I;

  printf("Bits-I: ");
  print_bits_ulong(I); // for checking

  double s_hi  = p;
  double s_med = 0.0;
  double s_lo  = 0.0;

  double y_hi;
  double y_lo;

  int i = 7; // num iterations - 1
  int j = 56; // current first bit

  int8_t w;
  bool w_pos;
  double t_hi, t_med, t_lo;

  while (i >= 0) {
    // This is checked and will be the bit values of the integer part of x 
    // starting with the larger bits (first 8)
    w = (int8_t)(I >> j);
    w_pos = (w > 0);

    printf("Current Bits-w: ");
    print_bits_u8(w); // for checking

    if (w != 0) {
      // get the t values
      calc_t(i, abs(w), &t_hi, &t_med, &t_lo);

      printf("t_hi: %.17g; t_med: %.17g; t_lo: %.17g\n", t_hi, t_med, t_lo);

      if (w_pos) {
        s_hi  += t_hi;
        s_med += t_med;
      } else {
        s_hi  -= t_hi;
        s_med -= t_med;
      }
    }

    // subtract the consumed part
    I -= ((uint64_t)w << j);

    // prepare next iteration
    i -= 1;
    j -= 8;
  }

  for (int i = 0; i < 8; i++) {
    if (w_pos) {
      s_lo += t_lo;
    } else {
      s_lo -= t_lo;
    }
  }

  if (fabs(s_hi) >= RANGE_MAX) {
    int k = (int)(round(fabs(s_hi) / RANGE_MAX));

    double c_hi = (int)((k * RANGE_MAX) / pow(2, -49)) * pow(2, -49);
    double c_med = (int)((k * RANGE_MAX - c_hi) / pow(2, -99)) * pow(2, -99);
    double c_lo = (k * RANGE_MAX - c_hi - c_med);

    printf("k: %d; s_hi: %.10f; c_hi: %.10f; c_med: %.10f; c_lo: %.10f\n", k, s_hi, c_hi, c_med, c_lo);

    if (s_hi > 0) {
      s_hi += c_hi;
      s_med += c_med;
      s_lo += c_lo;
    } else {
      s_hi -= c_hi;
      s_med -= c_med;
      s_lo -= c_lo;
    }
  }

  if (fabs(s_hi) > pow(2, -14)){
    double tmp = s_med + s_lo;
    // fast2sum
    y_hi = s_hi + tmp;
    double z = y_hi - s_hi;
    y_lo = tmp - z;

  } else if (s_hi == 0.0) {
    // fast2sum
    y_hi = s_med + s_lo;
    double z = y_hi - s_med;
    y_lo = s_lo - z;

  } else {
    // fast2sum
    y_hi = s_hi + s_med;
    double z = y_hi - s_hi;
    y_lo = s_med - z + s_lo;
  }

  //printf("S_HI: %.17g\n", s_hi);
  // printf("(y_hi: %.10f; y_lo: %.10f", y_hi, y_lo);

  return y_hi + y_lo;
}

typedef int int4;
typedef union { int4 i[2]; double x; double d; } mynumber;

static const double s1 = -0x1.5555555555555p-3;   /* -0.16666666666666666     */
static const double s2 = 0x1.1111111110ECEp-7;    /*  0.0083333333333323288   */
static const double s3 = -0x1.A01A019DB08B8p-13;  /* -0.00019841269834414642  */
static const double s4 = 0x1.71DE27B9A7ED9p-19;   /*  2.755729806860771e-06   */
static const double s5 = -0x1.ADDFFC2FCDF59p-26;  /* -2.5022014848318398e-08  */
static const double aa = -0x1.5558000000000p-3;   /* -0.1666717529296875      */
static const double bb = 0x1.5555555556E24p-18;   /*  5.0862630208387126e-06  */
static const double big = 0x1.8000000000000p45;   /*  52776558133248          */
static const double hp0 = 0x1.921FB54442D18p0;    /*  1.5707963267948966      */
static const double hp1 = 0x1.1A62633145C07p-54;  /*  6.123233995736766e-17   */
static const double mp1 = 0x1.921FB58000000p0;    /*  1.5707963407039642      */
static const double mp2 = -0x1.DDE973C000000p-27; /* -1.3909067564377153e-08  */
static const double mp3 = -0x1.CB3B399D747F2p-55; /* -4.9789962505147994e-17  */
static const double pp3 = -0x1.CB3B398000000p-55; /* -4.9789962314799099e-17  */
static const double pp4 = -0x1.d747f23e32ed7p-83; /* -1.9034889620193266e-25  */
static const double hpinv = 0x1.45F306DC9C883p-1; /*  0.63661977236758138     */
static const double toint = 0x1.8000000000000p52; /*  6755399441055744        */

/* Reduce range of x to within PI/2 with abs (x) < 105414350.  The high part
   is written to *a, the low part to *da.  Range reduction is accurate to 136
   bits so that when x is large and *a very close to zero, all 53 bits of *a
   are correct.  */
static __always_inline int4
reduce_sincos (double x, double *a, double *da)
{
  mynumber v;

  double t = (x * hpinv + toint);
  double xn = t - toint;
  v.x = t;
  double y = (x - xn * mp1) - xn * mp2;
  int4 n = v.i[0] & 3;

  double b, db, t1, t2;
  t1 = xn * pp3;
  t2 = y - t1;
  db = (y - t2) - t1;

  t1 = xn * pp4;
  b = t2 - t1;
  db += (t2 - b) - t1;

  *a = b;
  *da = db;
  return n;
}


int get_quadrant(double x) {
  int quadrant = 0;
  double reduced_range;
  double w = fabs(x);

  if (w < 8.0) {
    reduced_range = fmod(x, RANGE_MAX);
    quadrant = floor(reduced_range / M_PI_2);
    quadrant = (quadrant + 4) % 4; // prevents problems with negative quadrants.
  }
  if (8.0 <= w && w < pow(2, 63) - 1) {
    printf("\n\n Starting Brisbarre for x=%.17g\n", w);

    reduced_range = brisebarre_range_reduction(w);

    printf("(Brisebarre %.10f) ", reduced_range);

    quadrant = floor(reduced_range / M_PI_2);
  }

  return quadrant % 4;
}


// ######################## Test with ugly values #######################################

// Example usage
int main(void) {

  UglyCaseTriplet ugly_tests[] = {
    {0.0, "0", 0}, 
    {M_PI_2, "pi / 2", 1},
    {M_PI, "pi", 2},
    {M_PI * 1.5, "1.5 * pi", 3},
    {M_PI * 2.0, "2.0 * pi", 0},

    {-M_PI_2, "- pi / 2", 3},
    {-M_PI, "- pi", 2},
    {-M_PI * 1.5, "- 1.5 * pi", 1},
    {-M_PI * 2.0, "- 2.0 * pi", 0},

    {nextafter(0.0, -INFINITY), "0 - eps", 3}, 
    {nextafter(M_PI_2, -INFINITY), "pi / 2 - eps", 0},
    {nextafter(M_PI, -INFINITY), "pi - eps", 1},
    {nextafter(M_PI * 1.5, -INFINITY), "1.5 * pi - eps", 2},
    {nextafter(M_PI * 2.0, -INFINITY), "2.0 * pi - eps", 3},

    {nextafter(0.0, INFINITY),  "0 + eps", 0}, 
    {nextafter(M_PI_2, INFINITY), "pi / 2 + eps", 1},
    {nextafter(M_PI, INFINITY), "pi + eps", 2},
    {nextafter(M_PI * 1.5, INFINITY), "1.5 * pi + eps", 3},
    {nextafter(M_PI * 2.0, INFINITY), "2.0 * pi + eps", 0},

    {nextafter(-M_PI_2, INFINITY), "-pi / 2 + eps", 3},
    {nextafter(-M_PI, INFINITY), "-pi + eps", 2},
    {nextafter(-M_PI * 1.5, INFINITY), "-1.5 * pi + eps", 1},
    {nextafter(-M_PI * 2.0, INFINITY), "-2.0 * pi + eps", 0},

    {nextafter(M_PI * 16.5, INFINITY), "16.5 * pi + eps", 1},
    {nextafter(M_PI * 17.0, INFINITY), "17.0 * pi + eps", 2},
    {nextafter(M_PI * 17.5, INFINITY), "17.5 * pi + eps", 3},
    {nextafter(M_PI * 18.0, INFINITY), "18.0 * pi + eps", 0},


    {nextafter(- M_PI * 16.5, INFINITY), "- 16.5 * pi + eps", 3},
    {nextafter(- M_PI * 17.0, INFINITY), "- 17.0 * pi + eps", 2},
    {nextafter(- M_PI * 17.5, INFINITY), "- 17.5 * pi + eps", 1},
    {nextafter(- M_PI * 18.0, INFINITY), "- 18.0 * pi + eps", 0},

    {44.5, "44.5", 0},
    {46.25, "46.25", 1},
    {47.5, "47.5", 2},
    {49.75, "49.75", 3},

    {255.0, "255.0", 2},

    {302.0, "302.0", 0},
    {305.5, "305.5", 2},

    // The following examples have at the start of all 8 columns as long values
    // 0, which should make the calculation of t exact. 
    {6316927.0, "6316927", 0},
    {6324095.0, "6324095", 3},
    {1077952576.0, "1077952576", 3},

    // The following have at the start of some 8 columns as long a 1, 
    // which should destroy the calculation of t
    {6356612.0, "6356612.0", 0},
    {6356612.0, "6356612.0", 0},
    {6356630.0, "6356630.0", 0},
    {6324223.0, "6324223.0", 1},
    {6356632.0, "6356632.0", 1},
    {6356991.0, "6356991.0", 2},
    {6356608.0, "6356608.0", 2},
    {6356628.0, "6356628.0", 3},
    {6356635.0, "6356635.0", 3},
    {23788.5, "23788.5", 0},

    //  {pow(2, 57) * M_PI, "pi * 2^58", 0},

    //  // Those are the edge cases where it stops working 
    //  {nextafter(M_PI * 2.0 * INT_MAX, -INFINITY), "2 * pi * intmax - eps", 3},

    //  {M_PI * 2.0 * INT_MAX, "2 * pi * intmax", 0},
    //  {nextafter(M_PI * 2.0 * INT_MAX, INFINITY), "2 * pi * intmax + eps", 0},
    //  {M_PI * 2.0 * INT_MAX + 1, "2 * pi * intmax + 1", 0},
    //  {M_PI * 2.0 * INT_MAX + M_PI, "pi * (2 * intmax + 1)", 2},

    //  // Converted large number to bits then took those and put them into a bit to decimal calculator (https://numeral-systems.com/ieee-754-converter/)
    //  // then put the result into an full precision calculator and checked the quadrant (https://www.mathsisfun.com/calculator-precision.html)
    {180899997887870864251481058976852284292172291974620700999680.0, "approx 1.8 * 10^59", 0},
    //  {9999999978010827105118843667605086770705495962419200.0, "approx 9 * 10^51", 1},
    //  {999999997801082606665947196063956106460623011840.0, "approx 9 * 10^47", 1},
    //  {679899997801082676275550042794554161600728662016.0, "approx 6.8 * 10^47", 1},
    //  {579899997887998780515705500125536825728609841840128.0, "approx 5.8 * 10^50", 1},
    //  {3698999978879987753059268883445507785324923154088429879296.0, "approx 3.7 * 10^57", 2},
    //  {1408999978878708665079569693857962203230771643241038251115788338092413288448.0, "approx 1.4 * 10^75", 2},
    //  {140899997887870848106493420455200807689536941778985583455291404091450082399540603781120.0, "approx 1.4 * 10^86", 2},
    //  {61000687969105996493672285664758005575629682910825211579427204492165120.0, "approx 6.1 * 10^70", 2}
    {nextafter(0.0, -INFINITY), "0 - eps", 3}, 
  };

  int n = sizeof(ugly_tests) / sizeof(ugly_tests[0]); 


  for (int i = 0; i < n; i++) {
    int q = get_quadrant(ugly_tests[i].value);
    printf("  %-25s is in quadrant %d", ugly_tests[i].name, q);
    printf(" Correct: %-5s (%d)\n", (q == ugly_tests[i].quadrant) ? "True" : "\033[31mFalse\033[0m", ugly_tests[i].quadrant); // ]]
                                                                                                                                  //
    double a;
    double da;
    int4 n = reduce_sincos(ugly_tests[i].value, &a, &da);
    printf("The result of the glibc-implementation is: Quadrant: %d; a: %.17g; da: %.17g  this is correct? %-5s (%d)\n\n", n, a, da, (n == ugly_tests[i].quadrant) ? "True" : "\033[31mFalse\033[0m", ugly_tests[i].quadrant);
  }

  // print_double_bits(ugly_tests[1].value);
  

  // random_test();
  return 0;
}
