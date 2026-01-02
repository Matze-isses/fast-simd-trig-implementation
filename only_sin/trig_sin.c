#include "trig_sin.h"
#include <math.h>

const double TAYLOR_COEFF_SIN[] = {
  1,
  -0.16666666666666666,
  0.0083333333333333332,
  -0.00019841269841269841,
  2.7557319223985893e-06,
  -2.505210838544172e-08,
  1.6059043836821613e-10,
  -7.6471637318198164e-13,
  2.8114572543455206e-15,
  -8.2206352466243295e-18,
  1.9572941063391263e-20,
  -3.8681701706306835e-23,
  6.4469502843844736e-26,
  -9.183689863795546e-29,
  1.1309962886447718e-31,
  -1.2161250415535181e-34,
  1.1516335620771951e-37,
  -9.6775929586318907e-41,
  7.2654601791530714e-44,
  -4.9024697565135435e-47,
  2.9893108271424046e-50,
};

const double RANG_REDUCTION_CORRECTION = 3.8981718325193755e-17;
const double RANGE_MAX_SIN = M_PI * 2.0;
const double MED_RANGE_SIN = M_PI;
const double SMALL_RANGE_SIN = M_PI_2;

const double ONE_OVER_RANGE_SIN = 1 / RANGE_MAX_SIN;
const double ONE_OVER_MED_RANGE_SIN = 1 / MED_RANGE_SIN;
const double ONE_OVER_SMALL_RANGE_SIN = 1 / SMALL_RANGE_SIN;

const double RANGE_CENTER_SIN = 0;

const int SIZE_TAYLOR_COEFF = 20;

__m512d SIN(__m512d X) {
    int taylor_degree = 19;

    const int taylor_last_coeff = taylor_degree - 1;
    const int taylor_loop_iteration = taylor_degree - 2;

    const SDOUBLE range_reduction_correction = LOAD_DOUBLE(RANG_REDUCTION_CORRECTION);

    const SDOUBLE spi = LOAD_DOUBLE(M_PI);

    const SDOUBLE two_pi = LOAD_DOUBLE(RANGE_MAX_SIN);
    const SDOUBLE one_over_2_pi = LOAD_DOUBLE(ONE_OVER_RANGE_SIN);

    const SDOUBLE med_range = LOAD_DOUBLE(MED_RANGE_SIN);
    const SDOUBLE one_over_med_range = LOAD_DOUBLE(ONE_OVER_MED_RANGE_SIN);

    const SDOUBLE one_over_small_range = LOAD_DOUBLE(ONE_OVER_SMALL_RANGE_SIN);

    const SDOUBLE small_range = LOAD_DOUBLE(SMALL_RANGE_SIN);
    const SDOUBLE center_point = LOAD_DOUBLE(RANGE_CENTER_SIN);

    const SDOUBLE quadrant_multiplier = LOAD_DOUBLE(-2.0);
    const SDOUBLE ones = LOAD_DOUBLE(1.0);

    // works but is potentially negative
    const SDOUBLE ranges_away = MUL_DOUBLE_S(X, one_over_2_pi);
    const SDOUBLE num_ranges_away = FLOOR_DOUBLE_S(ranges_away);
    const SDOUBLE range_multiple = MUL_DOUBLE_S(num_ranges_away, two_pi);

    SDOUBLE in_outer_range = SUB_DOUBLE_S(X, range_multiple);
    SDOUBLE correction_term = MUL_DOUBLE_S(X, range_reduction_correction);
    in_outer_range = SUB_DOUBLE_S(in_outer_range, correction_term);

    // Gives Sign of the Result
    const SDOUBLE medium_ranges_away = MUL_DOUBLE_S(in_outer_range, one_over_med_range);
    const SDOUBLE sign = FLOOR_DOUBLE_S(medium_ranges_away);

    // Gives Quadrant of the result
    const SDOUBLE small_ranges_away = MUL_DOUBLE_S(in_outer_range, one_over_small_range);
    const SDOUBLE q = FLOOR_DOUBLE_S(small_ranges_away);
    const SDOUBLE q1 = ABS_PD(SUB_DOUBLE_S(q, ones));
    const SDOUBLE q2 = MUL_DOUBLE_S(q1, q1);
    const SDOUBLE q3 = ABS_PD(SUB_DOUBLE_S(q2, ones));

    // q3 * pi gives the mirroring points, where 0 and 2 do not need to be mirrored
    const SDOUBLE mirroring = MUL_DOUBLE_S(spi, q3);

    // all values mirroringare eighter in the first or in the thierd quadrant
    in_outer_range = SUB_DOUBLE_S(mirroring, in_outer_range);
    in_outer_range = ABS_PD(in_outer_range);

    const SDOUBLE small_ranges_away1 = MUL_DOUBLE_S(in_outer_range, one_over_small_range);
    const SDOUBLE q11 = FLOOR_DOUBLE_S(small_ranges_away1);
    const SDOUBLE initial_move = MUL_DOUBLE_S(small_range, q11);
    const SDOUBLE small_subtraction_amount = MUL_DOUBLE_S(q11, small_range);
    const SDOUBLE in_range = SUB_DOUBLE_S(in_outer_range, small_subtraction_amount);

    SDOUBLE result = LOAD_DOUBLE(TAYLOR_COEFF_SIN[taylor_last_coeff]);
    const SDOUBLE x_square = MUL_DOUBLE_S(in_range, in_range);

    for (int j = taylor_loop_iteration; j >= 0; j-=1) {
      SDOUBLE coeff = LOAD_DOUBLE(TAYLOR_COEFF_SIN[j]);
      result = FMADD_PD(result, x_square, coeff);
    }

    // to uneven the degrees
    result = MUL_DOUBLE_S(result, in_range);

    const SDOUBLE multiplied_quadrants = MUL_DOUBLE_S(sign, quadrant_multiplier);
    const SDOUBLE quadrant_evaluation = ADD_DOUBLE_S(multiplied_quadrants, ones);
    const SDOUBLE quadrant_evaluated_result = MUL_DOUBLE_S(result, quadrant_evaluation);

    return quadrant_evaluated_result;
}
