#include "clock_utils.h"
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define MAX_EXACT_INT (9007199254740992)
#define LONG_RANGE_END (57952155664616982739.0)

const double END_FRIENDLY_RANGE = M_PI * 2.0 * INT_MAX;


typedef struct {
  double value;
  const char *name;
  int quadrant;
} UglyCaseTriplet;


// ################################## Printing Bits #########################################################################

void print_double_bits(double value) {
  uint64_t bits;
  union {
    double d;
    uint64_t u;
  } conv;

  conv.d = value;
  bits = conv.u;
  printf("%d ", (int)((bits >> 63) & 1));

  for (int i = 62; i >= 52; i--) {
    printf("%d", (int)((bits >> i) & 1));
  }

  printf(" ");

  for (int i = 51; i >= 0; i--) {
    printf("%d", (int)((bits >> i) & 1));
  }

  printf("\n");
}

void print_bits_ulong(unsigned long x) {
    int nbits = sizeof(unsigned long) * 8;
    for (int i = nbits - 1; i >= 0; i--) {
        putchar( (x & (1ul << i)) ? '1' : '0' );
        // optional: add a space or separator every 8 bits
        if (i % 8 == 0 && i != 0) putchar(' ');
    }
    printf("\n");
}

void print_bits_u8(uint8_t x) {
    for (int i = 7; i >= 0; i--) {
        putchar( (x & (1u << i)) ? '1' : '0' );
    }
    printf("\n");
}

// ################################## brisbarre #########################################################################


void fast2sum(double a, double b, double *s, double *r) {
    *s = a + b;
    double z = *s - a;
    *r = b - z;
}

static inline double quantize_pow2(double v, int k) {
  return nearbyint(v * ldexp(1.0, -k)) * ldexp(1.0, k);
}

static inline double mod_pi_over_2(double v) {
  // bring v mod (pi/2) into [0, pi/2)
  double r = fmod(v, 2.0 * M_PI);
  if (r < 0) r += 2.0 * M_PI;
  return r;
}

double t_hi(int i, int8_t w) {
  // Th(i,w): nearest multiple of 2^-49 to ((2^(8i) * w) mod (pi/2))
  double tmp = mod_pi_over_2(ldexp(w, 8 * i));
  return quantize_pow2(tmp, -49);
}

double t_med(int i, int8_t w) {
  // Tmed(i,w): nearest multiple of 2^-99 to (tmp - Th)
  double tmp = mod_pi_over_2(ldexp(w, 8 * i));
  double thi = t_hi(i, w);
  double rem = tmp - thi;
  return quantize_pow2(rem, -99);
}

double t_lo(int i, int8_t w) {
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
  long I = (long) llround(x);
  double p = x - I;

  // print_bits_ulong(I); // for checking

  double s_hi  = p;
  double s_med = 0.0;
  double s_lo  = 0.0;

  int i = 7; // num iterations - 1
  int j = 56; // current first bit
  int8_t w;

  while (i >= 0) {
    // This is checked and will be the bit values of the integer part of x 
    // starting with the larger bits (first 8)
    w = (int8_t)(I >> j);

    // print_bits_u8(w); // for checking

    s_hi  += t_hi(i, abs(w));
    s_med += t_med(i, abs(w));
    // subtract the consumed part
    I -= ((uint64_t)w << j);

    i -= 1;
    j -= 8;
  }

  for (int i = 0; i < 8; i++) {
    s_lo  = t_lo(i, abs(w));
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


int get_quadrant(double x) {
  int quadrant = 0;
  double reduced_range;
  double w = fabs(x);

  if (w < 8.0) {
    reduced_range = fmod(x, 2.0 * M_PI);
    quadrant = floor(reduced_range / M_PI_2);
    quadrant = (quadrant + 4) % 4; // prevents problems with negative quadrants.
  }
  if (8.0 <= w && w < pow(2, 63) - 1) {
    reduced_range = brisebarre_range_reduction(w);
    printf("(Brisebarre %.10f) ", reduced_range);

    quadrant = floor(reduced_range / M_PI_2);
    
    if (x < 0) {
      quadrant = abs(quadrant - 3);
    }
  }

  return quadrant;
}


// ######################## Test with uniform numbers for speed #######################################

// Helper: generate a uniform random double in [0, 1]
static inline double uniform01(void) {
  unsigned long long high = (unsigned long long)rand();
  unsigned long long low = (unsigned long long)rand();
  unsigned long long combined = (high << 31) | low;
  return (double)combined / (double)((1ULL << 62) - 1);
}

// Main function: fill vector with uniform random doubles in [lower, upper]
void fill_uniform(double lower, double upper, size_t n, double *vec) {
  if (upper < lower) {
    fprintf(stderr, "Error: upper bound < lower bound.\n");
    return;
  }
  double range = upper - lower;
  for (size_t i = 0; i < n; i++) {
    vec[i] = lower + uniform01() * range;
  }
}

int random_test() {
  srand((unsigned)time(NULL));

  size_t n = 1000000000;
  double lower = 6.5;
  double upper = 1000000000.0;

  double *vec = malloc(n * sizeof(double));

  if (!vec) {
    perror("malloc");
    return 1;
  }

  fill_uniform(lower, upper, n, vec);

  START_CLOCK;

  for (int i = 0; i < n; i++) {
    get_quadrant(vec[i]);
  }

  END_CLOCK("Quadrant Calculation WARMUP");

  START_CLOCK;

  for (int i = 0; i < n; i++) {
    get_quadrant(vec[i]);
  }

  END_CLOCK("Quadrant Calculation       ");

  free(vec);
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

    {23788.5, "23788.5", 0},

    {pow(2, 57) * M_PI, "pi * 2^58", 0},

    // Those are the edge cases where it stops working 
    {nextafter(M_PI * 2.0 * INT_MAX, -INFINITY), "2 * pi * intmax - eps", 3},

    {M_PI * 2.0 * INT_MAX, "2 * pi * intmax", 0},
    {nextafter(M_PI * 2.0 * INT_MAX, INFINITY), "2 * pi * intmax + eps", 0},
    {M_PI * 2.0 * INT_MAX + 1, "2 * pi * intmax + 1", 0},
    {M_PI * 2.0 * INT_MAX + M_PI, "pi * (2 * intmax + 1)", 2},

    // Converted large number to bits then took those and put them into a bit to decimal calculator (https://numeral-systems.com/ieee-754-converter/)
    // then put the result into an full precision calculator and checked the quadrant (https://www.mathsisfun.com/calculator-precision.html)
    {180899997887870864251481058976852284292172291974620700999680.0, "approx 1.8 * 10^59", 0},
    {9999999978010827105118843667605086770705495962419200.0, "approx 9 * 10^51", 1},
    {999999997801082606665947196063956106460623011840.0, "approx 9 * 10^47", 1},
    {679899997801082676275550042794554161600728662016.0, "approx 6.8 * 10^47", 1},
    {579899997887998780515705500125536825728609841840128.0, "approx 5.8 * 10^50", 1},
    {3698999978879987753059268883445507785324923154088429879296.0, "approx 3.7 * 10^57", 2},
    {1408999978878708665079569693857962203230771643241038251115788338092413288448.0, "approx 1.4 * 10^75", 2},
    {140899997887870848106493420455200807689536941778985583455291404091450082399540603781120.0, "approx 1.4 * 10^86", 2},
    {61000687969105996493672285664758005575629682910825211579427204492165120.0, "approx 6.1 * 10^70", 2}
  };

  int n = sizeof(ugly_tests) / sizeof(ugly_tests[0]); 

  for (int i = 0; i < n; i++) {
    int q = get_quadrant(ugly_tests[i].value);
    printf("  %-25s is in quadrant %d", ugly_tests[i].name, q);
    printf(" Correct: %-5s (%d)\n", (q == ugly_tests[i].quadrant) ? "True" : "\033[31mFalse\033[0m", ugly_tests[i].quadrant); // ]]
  }

  // print_double_bits(ugly_tests[1].value);


  // random_test();
  return 0;
}
