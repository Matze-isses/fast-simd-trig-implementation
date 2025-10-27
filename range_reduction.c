#include "tests/clock_utils.h"
#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "range_reduction.h"
#include "bit_printing.h"

#define MAX_EXACT_INT (9007199254740992)
#define LONG_RANGE_END (57952155664616982739.0)

const double END_FRIENDLY_RANGE = M_PI * 2.0 * INT_MAX;

const double RANGE_MAX = M_PI * 2.0;
const double TWO_POW_NEG_49 = pow(2, -49);
const double TWO_POW_NEG_99 = pow(2, -99);
const double ONE_OVER_RANGE = 1 / RANGE_MAX;
const double ONE_OVER_PI_2 = 1 / M_PI_2;
const int MAX_SIMD_DOUBLES = (int)(SIMD_LENGTH / 64);
const int MAX_SIMD_FLOAT = (int)(SIMD_LENGTH / 32);

const int TAYLOR_DEGREE = 12;
const int TAYLOR_LAST_COEFF = TAYLOR_DEGREE - 1;
const int TAYLOR_LOOP_INTERATIONS = TAYLOR_DEGREE - 2;
const double TAYLOR_COEFF_SIN[] = {
  0.70710678118654746,
  0.70710678118654757,
  -0.35355339059327373,
  -0.11785113019775793,
  0.029462782549439476,
  0.0058925565098878968,
  -0.00098209275164798252,
  -0.00014029896452114038,
  1.7537370565142544e-05,
  1.9485967294602834e-06,
  -1.948596729460283e-07,
  -1.7714515722366212e-08,
  1.4762096435305173e-09,
  1.1355458796388596e-10,
  -8.1110419974204248e-12,
  -5.4073613316136172e-13,
  3.3796008322585101e-14,
  1.98800048956383e-15,
  -1.1044447164243498e-16,
  -5.8128669285492101e-18
};

const double TAYLOR_COEFF_COS[] = {
0.70710678118654757,
-0.70710678118654746,
-0.35355339059327379,
0.11785113019775791,
0.029462782549439483,
-0.0058925565098878951,
-0.00098209275164798274,
0.00014029896452114036,
1.7537370565142548e-05,
-1.948596729460283e-06,
-1.9485967294602832e-07,
1.7714515722366209e-08,
1.4762096435305176e-09,
-1.1355458796388595e-10,
-8.1110419974204264e-12,
5.4073613316136162e-13,
3.3796008322585108e-14,
-1.9880004895638296e-15,
-1.10444471642435e-16,
5.8128669285492093e-18
};


int get_reduced_range(double x, int *quadrant, double *reduced_range) {
  int n;

  n = floor(x * ONE_OVER_RANGE);
  *reduced_range = x - n * RANGE_MAX;
  *quadrant = floor(*reduced_range * ONE_OVER_PI_2);
  *quadrant = (*quadrant < 0) ? ((*quadrant + 4) % 4) : *quadrant;
}

void get_simd_quadrant_float(float *src, float *quad, float *range) {
  SFLOAT x = LOAD_FLOAT_VEC(src);
  SFLOAT two_pi = LOAD_FLOAT((float) RANGE_MAX);
  SFLOAT one_over_2_pi = LOAD_FLOAT((float) ONE_OVER_RANGE);
  SFLOAT one_over_pi_2 = LOAD_FLOAT((float) ONE_OVER_PI_2);

  SFLOAT n = MUL_FLOAT_S(x, one_over_2_pi);
  n = FLOOR_FLOAT_S(n);

  SFLOAT range_multiple = MUL_FLOAT_S(n, two_pi);
  SFLOAT in_range = SUB_FLOAT_S(x, range_multiple); // in [0, 2 * pi]

  n = MUL_FLOAT_S(in_range, one_over_pi_2);
  n = FLOOR_FLOAT_S(n);

  SIMD_TO_FLOAT_VEC(quad, n);
  SIMD_TO_FLOAT_VEC(range, in_range);
}

void sin_simd_float(float *x, float *res, size_t n, float prec) {
  float *quadrants = malloc(MAX_SIMD_FLOAT * sizeof(float));
  float *reduced_range = malloc(MAX_SIMD_FLOAT * sizeof(float));
  float *partial_values = malloc(MAX_SIMD_FLOAT * sizeof(float));
  
  for (int i = 0; i < (int)n; i+= MAX_SIMD_FLOAT) {
    for (int j = 0; j < MAX_SIMD_FLOAT; j++) {
      partial_values[j] = x[i+j];
    }

    get_simd_quadrant_float(partial_values, quadrants, reduced_range);
    
    for (int j = 0; j < MAX_SIMD_FLOAT; j++) {
      res[i+j] = quadrants[j];
    }
  }

  free(quadrants);
  free(reduced_range);
  free(partial_values);
}


double taylor_eval(double x, double a, const double coeffs[], int n) {
    double result = coeffs[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        result = result * (x - a) + coeffs[i];
    }
    return result;
}

void sin_simd(double *input, double *res, size_t n, float prec) {
  double *partial_values = malloc(MAX_SIMD_DOUBLES * sizeof(double));
  double *quadrants = malloc(MAX_SIMD_DOUBLES * sizeof(double));

  const SDOUBLE two_pi = LOAD_DOUBLE(RANGE_MAX);
  const SDOUBLE one_over_2_pi = LOAD_DOUBLE(ONE_OVER_RANGE);
  const SDOUBLE one_over_pi_2 = LOAD_DOUBLE(ONE_OVER_PI_2);
  const SDOUBLE simd_pi_2 = LOAD_DOUBLE(M_PI_2);
  const SDOUBLE center_point = LOAD_DOUBLE(0.5 * M_PI_2);
  
  #pragma omp parallel for
  for (int i = 0; i < (int)n; i+= MAX_SIMD_DOUBLES) {
    for (int j = 0; j < MAX_SIMD_DOUBLES; j++) {
      partial_values[j] = input[i+j];
    }

    SDOUBLE x   = LOAD_DOUBLE_VEC(partial_values);


    // works but is potentially negative
    SDOUBLE ranges_away = MUL_DOUBLE_S(x, one_over_2_pi);
    SDOUBLE num_ranges_away = FLOOR_DOUBLE_S(ranges_away);
    SDOUBLE range_multiple = MUL_DOUBLE_S(num_ranges_away, two_pi);
    SDOUBLE in_outer_range = SUB_DOUBLE_S(x, range_multiple); // in [0, 2 * pi]

    SDOUBLE small_ranges_away = MUL_DOUBLE_S(in_outer_range, one_over_pi_2);
    SDOUBLE simd_quadrants = FLOOR_DOUBLE_S(small_ranges_away);
    SDOUBLE small_subtraction_amount = MUL_DOUBLE_S(simd_quadrants, simd_pi_2);
    SDOUBLE in_range = SUB_DOUBLE_S(in_outer_range, small_subtraction_amount);

    SDOUBLE centered_values = SUB_DOUBLE_S(in_range, center_point);


    SDOUBLE result = LOAD_DOUBLE(TAYLOR_COEFF_SIN[TAYLOR_LAST_COEFF]);

    for (int j = TAYLOR_LOOP_INTERATIONS; j >= 0; --j) {
      SDOUBLE coeff = LOAD_DOUBLE(TAYLOR_COEFF_SIN[j]);
      result = MUL_DOUBLE_S(result, centered_values);
      result = ADD_DOUBLE_S(result, coeff);
    }
    
    SIMD_TO_DOUBLE_VEC(quadrants, simd_quadrants); 
    SIMD_TO_DOUBLE_VEC(&res[i], result); 

    for (int j = 0; j < 4; j++) {
      switch ((int)quadrants[j]) {
        case 1:
          res[i+j] = sqrt(1 - res[i+j] * res[i+j]);
          break;
        case 2:
          res[i+j] = - res[i+j];
          break;
        case 3:
          res[i+j] = - sqrt(1 - res[i+j] * res[i+j]);
          break;
        default:
          break;
      }
    }
  }

 
  int num_left_over = (n % 4);

  #pragma omp parallel for
  for (int i = n - num_left_over; i < (int)n; i++) {
    double reduced_range;
    int quadrant;
    get_reduced_range(input[i], &quadrant, &reduced_range);
    res[i] = taylor_eval(reduced_range, 0.5 * M_PI_2, TAYLOR_COEFF_SIN, TAYLOR_DEGREE);

    if (quadrant == 1) {
      res[i] = sqrt(1 - (res[i] * res[i]));
    } else if (quadrant == 2) {
      res[i] = - res[i];
    } else if (quadrant == 3) {
      res[i] = - sqrt(1 - (res[i] * res[i]));
    }
  }
  

  free(partial_values);
  free(quadrants);
}

// gcc range_reduction.c bit_printing.c -o range_reduction -lm -mavx -O2 && ./range_reduction
