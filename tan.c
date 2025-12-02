#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>

#include <stdio.h>
#include <stdlib.h>
#include "trig_simd.h"

#include <string.h>   // just used for the error calculation
#include "./util/bit_printing.h"


double CORRECTION = 0.00000000000000006123233995736765;
double M_PI_8 = M_PI / 8;

// Working until 0.4
double TAYLOR_COEFF_TAN[] = {
  1.0000000000000000,
  0.33333333333333331,
  0.13333333333333333,
  0.053968253968253971,
  0.021869488536155203,
  0.0088632355299021973,
  0.0035921280365724811,
  0.0014558343870513183,
  0.00059002744094558595,
  0.00023912911424355248,
  9.6915379569294509e-05,
  3.9278323883316833e-05,
  1.5918905069328964e-05,
  6.4516892156554306e-06
};


#define N_FIRST_MID 17

static const double X_FIRST_MID[17] = {
    0.39269908169872414,
    0.3974311529539063,
    0.4054030240311012,
    0.4286863854186324,
    0.4576011079513348,
    0.487426901464288,
    0.5407991673806207,
    0.5723175108611952,
    0.6115367603429416,
    0.6435891799153872,
    0.6768825409954158,
    0.7002130267439945,
    0.7267833123145859,
    0.7490809891502078,
    0.7665916163094255,
    0.7796607724399203,
    0.7853981633974483
};

static const double Y_FIRST_MID[17] = {
    0.41421356237309503,
    0.4197684582474178,
    0.42917669979674156,
    0.4570320684430114,
    0.4924645633683938,
    0.5300875167665825,
    0.6005164649453173,
    0.6442430603066724,
    0.7012087722833138,
    0.7501376202189168,
    0.8035182551172894,
    0.8426526043282521,
    0.8891412771435833,
    0.9298813332611141,
    0.9630769477832821,
    0.9885905533499666,
    1.0000000000000000
};


double addon_lookup[256] = {
3.0531133177191805e-16,
2.220446049250313e-16,
3.885780586188048e-16,
1.2878587085651816e-13,
1.0187406473960436e-12,
3.785305402459471e-12,
6.422529175154068e-12,
1.87158066822235e-11,
3.1586400162098016e-11,
4.7865489349874224e-11,
5.759093202328813e-11,
8.885203683917098e-11,
1.1255296694656636e-10,
1.377497005350392e-10,
1.6403556291066934e-10,
1.9105284021492253e-10,
2.1850710130166817e-10,
2.4615698368535277e-10,
2.7380486766759304e-10,
3.0129632122566363e-10,
3.2850777653692376e-10,
3.5533898046224976e-10,
3.8172687233384295e-10,
3.950025861954032e-10,
4.329453462403876e-10,
4.5770931489386157e-10,
4.818874188572408e-10,
5.054627827405511e-10,
5.284351844991875e-10,
5.508051792446622e-10,
5.725776519582837e-10,
5.937618174911563e-10,
6.042903955005841e-10,
6.246012596022865e-10,
6.538847241444046e-10,
6.728271273459541e-10,
6.912447281237633e-10,
7.003763125013052e-10,
7.265643642284658e-10,
7.434972637554438e-10,
7.599663121027334e-10,
7.759869413703768e-10,
7.839218163496753e-10,
8.067353451934878e-10,
8.214945390605521e-10,
8.358609360215041e-10,
8.498507453325033e-10,
8.634722936662342e-10,
8.767406800558319e-10,
8.833129783170079e-10,
8.960836517246662e-10,
9.145477708472072e-10,
9.26521415145487e-10,
9.381975196731673e-10,
9.439786730069955e-10,
9.606979656240355e-10,
9.71541402883247e-10,
9.821266022669306e-10,
9.924600030686292e-10,
1.0025510421840522e-09,
1.0124086013973965e-09,
1.0220388979576e-09,
1.0314524789833968e-09,
1.0406473460733423e-09,
1.0496421509742504e-09,
1.058434229150862e-09,
1.067024579803899e-09,
1.071297051069564e-09,
1.083683809355307e-09,
1.091739587621987e-09,
1.099627278122739e-09,
1.1073483241474946e-09,
1.1111750408687726e-09,
1.118659831433888e-09,
1.1259874144187165e-09,
1.1366829699710479e-09,
1.1402061517173934e-09,
1.150480266609577e-09,
1.153864115366332e-09,
1.160494478291696e-09,
1.166988061740426e-09,
1.1733698457305763e-09,
1.1826883916654651e-09,
1.1887667517029854e-09,
1.194730203657457e-09,
1.2005891836253113e-09,
1.2063368082237957e-09,
1.2119809600363851e-09,
1.2175249697321533e-09,
1.2202723276288907e-09,
1.2256674564170567e-09,
1.2335623633674686e-09,
1.2387304515470987e-09,
1.2437976204537904e-09,
1.2487895162394125e-09,
1.253689596580898e-09,
1.2585094077977033e-09,
1.2632461743322665e-09,
1.2679061134335257e-09,
1.2724880038561537e-09,
1.2747616295882835e-09,
1.2792311654408195e-09,
1.2857948039624034e-09,
1.290087592309419e-09,
1.29431299011884e-09,
1.298475771349672e-09,
1.3025666101285083e-09,
1.3065992732208542e-09,
1.3105693197346113e-09,
1.3125394104918087e-09,
1.31832567085155e-09,
1.3221171935029474e-09,
1.3258526498916012e-09,
1.3277123844801508e-09,
1.333154364679956e-09,
1.3367298379307613e-09,
1.3402503551418476e-09,
1.3437194690268939e-09,
1.3471171955714567e-09,
1.3505111473577358e-09,
1.3538351550934635e-09,
1.3571133106182742e-09,
1.360343948597631e-09,
1.3635311768567249e-09,
1.3651144659121428e-09,
1.3682360799904814e-09,
1.3728336245577566e-09,
1.375850766649478e-09,
1.3788326036490162e-09,
1.3817642585678414e-09,
1.383314351954823e-09,
1.3875234294857819e-09,
1.390347281748916e-09,
1.3931330533623054e-09,
1.394517723518618e-09,
1.3985982372233252e-09,
1.4012797588947024e-09,
1.4039267526300137e-09,
1.4052423669141945e-09,
1.4091218192291421e-09,
1.4104049039787014e-09,
1.4141877668905067e-09,
1.4166749995325745e-09,
1.417911565937402e-09,
1.4215584265286907e-09,
1.4239557311057638e-09,
1.4263246139734065e-09,
1.4286645200201065e-09,
1.4298289219283333e-09,
1.432125085187863e-09,
1.4355215904870988e-09,
1.437755470234947e-09,
1.439962371563297e-09,
1.4410591608893242e-09,
1.4443005680320198e-09,
1.4464319741946952e-09,
1.4474921261609097e-09,
1.4506221779342354e-09,
1.4526841951578717e-09,
1.454724340987923e-09,
1.4567390627107102e-09,
1.4587325791737271e-09,
1.4607033360647392e-09,
1.4626529987182835e-09,
1.463624110797923e-09,
1.4664912617590176e-09,
1.4683788629454853e-09,
1.4693152250444541e-09,
1.4720945573643007e-09,
1.473922761618951e-09,
1.4757310928814604e-09,
1.4775222156870882e-09,
1.478411060240603e-09,
1.480175315649035e-09,
1.4819176996638816e-09,
1.4845000784191598e-09,
1.4862052699626815e-09,
1.4878821508190754e-09,
1.4895492617128525e-09,
1.4911979429044209e-09,
1.4928317471074593e-09,
1.4944471216082889e-09,
1.4960482852544033e-09,
1.4968463135645038e-09,
1.4992026509119682e-09,
1.499986579389656e-09,
1.502296731459296e-09,
1.5038199574490818e-09,
1.5053299717848745e-09,
1.5068252201544396e-09,
1.5083071458477093e-09,
1.5097761929538933e-09,
1.5105061645925844e-09,
1.5126661034869926e-09,
1.51409518256429e-09,
1.5155052768278665e-09,
1.5169109301993444e-09,
1.5182928247980954e-09,
1.5196692793040256e-09,
1.5210287473976791e-09,
1.5223775573502962e-09,
1.5230499084140092e-09,
1.5250404272748597e-09,
1.5256997887291845e-09,
1.5276554465870618e-09,
1.5289448596078614e-09,
1.5295820166016938e-09,
1.531496374163055e-09,
1.5327479285787149e-09,
1.5333726510746715e-09,
1.53522783374882e-09,
1.535842564237555e-09,
1.5376591111504467e-09,
1.5388679219796586e-09,
1.5400611896865257e-09,
1.540651273224114e-09,
1.542415084543336e-09,
1.543574490447952e-09,
1.5447330081741484e-09,
1.545875760733395e-09,
1.547007300040093e-09,
1.5481307347187112e-09,
1.54868973201161e-09,
1.550354400414733e-09,
1.5514484141831986e-09,
1.5525369878588435e-09,
1.5530794428286754e-09,
1.5541523623596731e-09,
1.555747641823757e-09,
1.5567996891618918e-09,
1.5578501821877921e-09,
1.5588808022215517e-09,
1.5599098679430767e-09,
1.5604211256459166e-09,
1.5614378678918683e-09,
1.5624427307514566e-09,
1.563941864901608e-09,
1.5649297413489194e-09,
1.5659105123688732e-09,
1.5668837338722597e-09,
1.567368457244811e-09,
1.5688083054854474e-09,
1.5697582123053166e-09,
1.5707031231215751e-09,
1.5716381529529144e-09,
1.5725684088252478e-09,
1.5734901159802916e-09,
1.5739520797808382e-09,
1.5748626847056357e-09,
1.576212493858975e-09,
1.5771118855312238e-09,
1.5775596384770552e-09,
1.5784434870269592e-09,
1.5793231167293698e-09,
1.5801955299821202e-09,
};

static double LAGRANGE_DEN_FIRST_MID[N_FIRST_MID][N_FIRST_MID];

void init_first_mid_lagrange_table(void)
{
    for (size_t i = 0; i < N_FIRST_MID; ++i) {
        for (size_t j = 0; j < N_FIRST_MID; ++j) {
            if (i == j) {
                LAGRANGE_DEN_FIRST_MID[i][j] = 0.0;
            } else {
                LAGRANGE_DEN_FIRST_MID[i][j] = 1.0 / (X_FIRST_MID[i] - X_FIRST_MID[j]);
            }
        }
    }
}

void start_of_range(double input, double *res) {
    double taylor = input;
    double x_square = taylor * taylor;
    double result = TAYLOR_COEFF_TAN[13];

    for (int j = 12; j >= 0; j-=1) {
      double coeff = TAYLOR_COEFF_TAN[j];
      result = result * x_square + coeff;
    }

    *res = result * taylor;

}

void end_of_range(double input, double *res) {
    double from_behind = M_PI_2 - input;
    double x_square = from_behind * from_behind;

    TAYLOR_COEFF_TAN[0] += CORRECTION * 1/from_behind;
    double result = TAYLOR_COEFF_TAN[13];

    for (int j = 12; j >= 0; j-=1) {
      double coeff = TAYLOR_COEFF_TAN[j];
      result = result * x_square + coeff;
    }

    result = result * from_behind;
    *res = 1 / result;
    TAYLOR_COEFF_TAN[0] = 1.0;
}

void first_mid_range(double input, double *res) {
    double result = 0.0;

    for (size_t i = 0; i < N_FIRST_MID; ++i) {
        double Li = 1.0;

        for (size_t j = 0; j < N_FIRST_MID; ++j) {
            if (j == i) continue;
            Li *= (input - X_FIRST_MID[j]) * LAGRANGE_DEN_FIRST_MID[i][j];
        }

        result += Y_FIRST_MID[i] * Li;
    }

    *res = result;
}

void sec_mid_range(double input, double *res) {
    double result = 0.0;
    double from_end = M_PI_2 - input;

    for (size_t i = 0; i < N_FIRST_MID; ++i) {
        double Li = 1.0;

        for (size_t j = 0; j < N_FIRST_MID; ++j) {
            if (j == i) continue;
            Li *= (from_end - X_FIRST_MID[j]) * LAGRANGE_DEN_FIRST_MID[i][j];
        }

        result += Y_FIRST_MID[i] * Li;
    }

    *res = 1/result;
}


#define PRINT_M128I(reg) do {                                   \
    int32_t vals[4];                                            \
    _mm_storeu_si128((__m128i*)vals, (reg));                    \
    printf(#reg " = [%d, %d, %d, %d]\n",                        \
           vals[0], vals[1], vals[2], vals[3]);                 \
} while (0)


void tan_simd(double *input, double *res, size_t n, int prec) {
  int simd_doubles = 4;
  // init_first_mid_lagrange_table();

  const SDOUBLE one_over_pi_8 = LOAD_DOUBLE(1/M_PI_8);
  const SDOUBLE m_pi_2 = LOAD_DOUBLE(M_PI_2);
  const SDOUBLE index_sum = LOAD_DOUBLE(M_PI_4+1.0/256.0);
  const SDOUBLE correction = LOAD_DOUBLE(CORRECTION);

  int last_taylor_coeff = 10;
  int taylor_loop_iteration = last_taylor_coeff - 1;
  double test_lookup[2] = {0.0, 0.1};

  const SDOUBLE half = LOAD_DOUBLE(0.5);
  const SDOUBLE one = LOAD_DOUBLE(1.0);
  const SDOUBLE two = LOAD_DOUBLE(2.0);
  
  for (int i = 0; i < (int) n; i += 4) {
    SDOUBLE x   = LOAD_DOUBLE_VEC(&input[i]);

    SDOUBLE from_behind = SUB_DOUBLE_S(m_pi_2, x);

    const SDOUBLE not_floored = MUL_DOUBLE_S(x, one_over_pi_8);
    const SDOUBLE quadrant = FLOOR_DOUBLE_S(not_floored);

    SDOUBLE in_q0 = SUB_DOUBLE_S(quadrant, two);
    in_q0 = ABS_PD(in_q0);
    in_q0 = MUL_DOUBLE_S(in_q0, half);
    in_q0 = FLOOR_DOUBLE_S(in_q0);

 
    SDOUBLE bot = SUB_DOUBLE_S(index_sum, x);
    SDOUBLE idouble = DIV_DOUBLE_S(one, bot);
    __m128i idx = _mm256_cvttpd_epi32(idouble);   

    SDOUBLE lin_addon = _mm256_i32gather_pd(addon_lookup, idx, 8);
    PRINT_FULL_M256D(lin_addon);
    x = ADD_DOUBLE_S(x, lin_addon);

    SDOUBLE in_q1 = SUB_DOUBLE_S(quadrant, one);
    in_q1 = ABS_PD(in_q1);
    in_q1 = MUL_DOUBLE_S(in_q1, half);
    in_q1 = FLOOR_DOUBLE_S(in_q1);

    
    /* ---- Calculation for first range ---- */
    const SDOUBLE x_square = MUL_DOUBLE_S(x, x);
    SDOUBLE result_q0 = LOAD_DOUBLE(TAYLOR_COEFF_TAN[last_taylor_coeff]);

    for (int j = taylor_loop_iteration; j >= 0; j-=1) {
      SDOUBLE coeff = LOAD_DOUBLE(TAYLOR_COEFF_TAN[j]);
      result_q0 = FMADD_PD(result_q0, x_square, coeff);
    }

    result_q0 = MUL_DOUBLE_S(result_q0, x);


    /* ---- Calculation for fourth range ---- */
    SDOUBLE from_behind_square = MUL_DOUBLE_S(from_behind, from_behind);

    // Calculation of the correction term
    SDOUBLE one_over_from_behind = DIV_DOUBLE_S(one, from_behind);
    SDOUBLE correction_term = MUL_DOUBLE_S(correction, one_over_from_behind);
    SDOUBLE first_coeff = ADD_DOUBLE_S(one, correction_term);

    SDOUBLE result_q3 = LOAD_DOUBLE(TAYLOR_COEFF_TAN[last_taylor_coeff]);

    // Important that first term is treated with correction
    for (int j = taylor_loop_iteration; j >= 1; j-=1) {
      SDOUBLE coeff = LOAD_DOUBLE(TAYLOR_COEFF_TAN[j]);
      result_q3 = FMADD_PD(result_q3, from_behind_square, coeff);
    }

    // add the corrected first coefficiant
    result_q3 = FMADD_PD(result_q3, from_behind_square, first_coeff);
    result_q3 = MUL_DOUBLE_S(result_q3, from_behind);
    result_q3 = DIV_DOUBLE_S(one, result_q3);


    /* ---- Final Addup ---- */


    /*
    PRINT_M256D(result_q0);
    PRINT_M256D(result_q1);
    PRINT_M256D(result_q2);
    PRINT_M256D(result_q3);
    */

    SDOUBLE result = LOAD_DOUBLE(0.0);
    result = FMADD_PD(result_q0, in_q0, result);
    // result = FMADD_PD(result_q0, in_q1, result);
    result = FMADD_PD(result_q3, in_q1, result);


    SIMD_TO_DOUBLE_VEC(&res[i], result);

  }

  int num_left_over = (n % 4);

  for (size_t i = n - num_left_over; i < (int)n; i++) {
    if (input[i] < M_PI_8) {
      start_of_range(input[i], &res[i]);

    } else if (input[i] < 2 * M_PI_8) {
      first_mid_range(input[i], &res[i]);

    } else if (input[i] < 3 * M_PI_8) {
      sec_mid_range(input[i], &res[i]);

    } else if (input[i] > 3 * M_PI_8){
      end_of_range(input[i], &res[i]);
    }
  }

  for (int i = 0; i < n; i++) {
      if (isnan(res[i])) {
        printf("NaN at index %d, input = %.17g\n", i, input[i]);
      }
  }

}


// gcc ./cmeasure/cbind_to_hw_thread.c ./cmeasure/cmeasure.c ./cmeasure/CrystalClockInC.c ./trig_simd.c ./tests/test_interface_tan.c ./tests/value_generation.c ./tests/trig_arb_comparison.c ./util/bit_printing.c ./tan.c -o test -lm -mavx -mavx2 -mfma -O2 -lflint -Wextra && ./test 10000000 0 1.570796326794896 10000000
//
// Compilation and run on laptop
//
//
// gcc ./cmeasure/cbind_to_hw_thread.c ./cmeasure/cmeasure.c ./cmeasure/CrystalClockInC.c ./trig_simd.c ./tests/test_interface_tan.c ./tests/value_generation.c ./tests/trig_arb_comparison.c ./util/bit_printing.c ./tan.c -o test -lm -mavx -mavx2 -mfma -O2 $(pkg-config --cflags --libs flint) -Wextra && ./test 10000000 0 1.570796326794896 10000000
//
