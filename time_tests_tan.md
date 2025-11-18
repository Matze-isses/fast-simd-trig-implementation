
# Time Needed for the Lagrange Polyromial (Quadrant 2 and 3)

n = 100000000

## Single Execution
TIME OC: 5954.621333

TIME CC: 6030.3210754342317


```c
void tan_simd(double *input, double *res, size_t n, int prec) {
  int simd_doubles = 4;
  init_first_mid_lagrange_table();

  const SDOUBLE one = LOAD_DOUBLE(1.0);
  
  for (int i = 0; i < (int) n; i += 4) {
    SDOUBLE x   = LOAD_DOUBLE_VEC(&input[i]);

    /* ---- Calculation for second range ---- */
    SDOUBLE result_q1 = LOAD_DOUBLE(0.0);

    for (size_t i = 0; i < N_FIRST_MID; ++i) {
        SDOUBLE Li = one;

        // to prevent if statements the loops are splitted
        for (size_t j = 0; j < i; ++j) {
            SDOUBLE x_first = LOAD_DOUBLE(X_FIRST_MID[j]);
            SDOUBLE lagrange_lookup = LOAD_DOUBLE(LAGRANGE_DEN_FIRST_MID[i][j]);

            SDOUBLE sub_x = SUB_DOUBLE_S(x, x_first);

            Li = MUL_DOUBLE_S(Li, sub_x);
            Li = MUL_DOUBLE_S(Li, lagrange_lookup);
        }

        for (size_t j = i+1; j < N_FIRST_MID; ++j) {
            SDOUBLE x_first = LOAD_DOUBLE(X_FIRST_MID[j]);
            SDOUBLE lagrange_lookup = LOAD_DOUBLE(LAGRANGE_DEN_FIRST_MID[i][j]);

            SDOUBLE sub_x = SUB_DOUBLE_S(x, x_first);

            Li = MUL_DOUBLE_S(Li, sub_x);
            Li = MUL_DOUBLE_S(Li, lagrange_lookup);
        }

        SDOUBLE y_first = LOAD_DOUBLE(Y_FIRST_MID[i]);

        SDOUBLE eval_inner = MUL_DOUBLE_S(y_first, Li);
        result_q1 = ADD_DOUBLE_S(result_q1, eval_inner);
    }

    /* ---- Calculation for thierd range ---- */
    SDOUBLE result_q2 = result_q1;

    for (size_t i = 0; i < N_FIRST_MID; ++i) {
        SDOUBLE Li = one;

        // to prevent if statements the loops are splitted
        for (size_t j = 0; j < i; ++j) {
            SDOUBLE x_first = LOAD_DOUBLE(X_FIRST_MID[j]);
            SDOUBLE lagrange_lookup = LOAD_DOUBLE(LAGRANGE_DEN_FIRST_MID[i][j]);

            SDOUBLE sub_x = SUB_DOUBLE_S(x, x_first);

            Li = MUL_DOUBLE_S(Li, sub_x);
            Li = MUL_DOUBLE_S(Li, lagrange_lookup);
        }

        for (size_t j = i+1; j < N_FIRST_MID; ++j) {
            SDOUBLE x_first = LOAD_DOUBLE(X_FIRST_MID[j]);
            SDOUBLE lagrange_lookup = LOAD_DOUBLE(LAGRANGE_DEN_FIRST_MID[i][j]);

            SDOUBLE sub_x = SUB_DOUBLE_S(x, x_first);

            Li = MUL_DOUBLE_S(Li, sub_x);
            Li = MUL_DOUBLE_S(Li, lagrange_lookup);
        }

        SDOUBLE y_first = LOAD_DOUBLE(Y_FIRST_MID[i]);

        SDOUBLE eval_inner = MUL_DOUBLE_S(y_first, Li);
        result_q2 = ADD_DOUBLE_S(result_q2, eval_inner);
    }

    result_q2 = DIV_DOUBLE_S(one, result_q2);
    SIMD_TO_DOUBLE_VEC(&res[i], result_q2);
  }

}
```


## Double Execution 

OWN TIME CC: 11908.752674338131

OWN TIME OC: 11935.889466000001

```c 

void tan_simd(double *input, double *res, size_t n, int prec) {
  int simd_doubles = 4;
  init_first_mid_lagrange_table();

  const SDOUBLE one = LOAD_DOUBLE(1.0);
  
  for (int i = 0; i < (int) n; i += 4) {
    SDOUBLE x   = LOAD_DOUBLE_VEC(&input[i]);

    /* ---- Calculation for second range ---- */
    SDOUBLE result_q1 = LOAD_DOUBLE(0.0);

    for (size_t i = 0; i < N_FIRST_MID; ++i) {
        SDOUBLE Li = one;

        // to prevent if statements the loops are splitted
        for (size_t j = 0; j < i; ++j) {
            SDOUBLE x_first = LOAD_DOUBLE(X_FIRST_MID[j]);
            SDOUBLE lagrange_lookup = LOAD_DOUBLE(LAGRANGE_DEN_FIRST_MID[i][j]);

            SDOUBLE sub_x = SUB_DOUBLE_S(x, x_first);

            Li = MUL_DOUBLE_S(Li, sub_x);
            Li = MUL_DOUBLE_S(Li, lagrange_lookup);
        }

        for (size_t j = i+1; j < N_FIRST_MID; ++j) {
            SDOUBLE x_first = LOAD_DOUBLE(X_FIRST_MID[j]);
            SDOUBLE lagrange_lookup = LOAD_DOUBLE(LAGRANGE_DEN_FIRST_MID[i][j]);

            SDOUBLE sub_x = SUB_DOUBLE_S(x, x_first);

            Li = MUL_DOUBLE_S(Li, sub_x);
            Li = MUL_DOUBLE_S(Li, lagrange_lookup);
        }

        SDOUBLE y_first = LOAD_DOUBLE(Y_FIRST_MID[i]);

        SDOUBLE eval_inner = MUL_DOUBLE_S(y_first, Li);
        result_q1 = ADD_DOUBLE_S(result_q1, eval_inner);
    }

    /* ---- Calculation for thierd range ---- */
    SDOUBLE result_q2 = result_q1;

    for (size_t i = 0; i < N_FIRST_MID; ++i) {
        SDOUBLE Li = one;

        // to prevent if statements the loops are splitted
        for (size_t j = 0; j < i; ++j) {
            SDOUBLE x_first = LOAD_DOUBLE(X_FIRST_MID[j]);
            SDOUBLE lagrange_lookup = LOAD_DOUBLE(LAGRANGE_DEN_FIRST_MID[i][j]);

            SDOUBLE sub_x = SUB_DOUBLE_S(x, x_first);

            Li = MUL_DOUBLE_S(Li, sub_x);
            Li = MUL_DOUBLE_S(Li, lagrange_lookup);
        }

        for (size_t j = i+1; j < N_FIRST_MID; ++j) {
            SDOUBLE x_first = LOAD_DOUBLE(X_FIRST_MID[j]);
            SDOUBLE lagrange_lookup = LOAD_DOUBLE(LAGRANGE_DEN_FIRST_MID[i][j]);

            SDOUBLE sub_x = SUB_DOUBLE_S(x, x_first);

            Li = MUL_DOUBLE_S(Li, sub_x);
            Li = MUL_DOUBLE_S(Li, lagrange_lookup);
        }

        SDOUBLE y_first = LOAD_DOUBLE(Y_FIRST_MID[i]);

        SDOUBLE eval_inner = MUL_DOUBLE_S(y_first, Li);
        result_q2 = ADD_DOUBLE_S(result_q2, eval_inner);
    }

    result_q2 = DIV_DOUBLE_S(one, result_q2);
    SIMD_TO_DOUBLE_VEC(&res[i], result_q2);
  }

  for (int i = 0; i < (int) n; i += 4) {
    SDOUBLE x   = LOAD_DOUBLE_VEC(&res[i]);

    /* ---- Calculation for second range ---- */
    SDOUBLE result_q1 = LOAD_DOUBLE(0.0);

    for (size_t i = 0; i < N_FIRST_MID; ++i) {
        SDOUBLE Li = one;

        // to prevent if statements the loops are splitted
        for (size_t j = 0; j < i; ++j) {
            SDOUBLE x_first = LOAD_DOUBLE(X_FIRST_MID[j]);
            SDOUBLE lagrange_lookup = LOAD_DOUBLE(LAGRANGE_DEN_FIRST_MID[i][j]);

            SDOUBLE sub_x = SUB_DOUBLE_S(x, x_first);

            Li = MUL_DOUBLE_S(Li, sub_x);
            Li = MUL_DOUBLE_S(Li, lagrange_lookup);
        }

        for (size_t j = i+1; j < N_FIRST_MID; ++j) {
            SDOUBLE x_first = LOAD_DOUBLE(X_FIRST_MID[j]);
            SDOUBLE lagrange_lookup = LOAD_DOUBLE(LAGRANGE_DEN_FIRST_MID[i][j]);

            SDOUBLE sub_x = SUB_DOUBLE_S(x, x_first);

            Li = MUL_DOUBLE_S(Li, sub_x);
            Li = MUL_DOUBLE_S(Li, lagrange_lookup);
        }

        SDOUBLE y_first = LOAD_DOUBLE(Y_FIRST_MID[i]);

        SDOUBLE eval_inner = MUL_DOUBLE_S(y_first, Li);
        result_q1 = ADD_DOUBLE_S(result_q1, eval_inner);
    }

    /* ---- Calculation for thierd range ---- */
    SDOUBLE result_q2 = result_q1;

    for (size_t i = 0; i < N_FIRST_MID; ++i) {
        SDOUBLE Li = one;

        // to prevent if statements the loops are splitted
        for (size_t j = 0; j < i; ++j) {
            SDOUBLE x_first = LOAD_DOUBLE(X_FIRST_MID[j]);
            SDOUBLE lagrange_lookup = LOAD_DOUBLE(LAGRANGE_DEN_FIRST_MID[i][j]);

            SDOUBLE sub_x = SUB_DOUBLE_S(x, x_first);

            Li = MUL_DOUBLE_S(Li, sub_x);
            Li = MUL_DOUBLE_S(Li, lagrange_lookup);
        }

        for (size_t j = i+1; j < N_FIRST_MID; ++j) {
            SDOUBLE x_first = LOAD_DOUBLE(X_FIRST_MID[j]);
            SDOUBLE lagrange_lookup = LOAD_DOUBLE(LAGRANGE_DEN_FIRST_MID[i][j]);

            SDOUBLE sub_x = SUB_DOUBLE_S(x, x_first);

            Li = MUL_DOUBLE_S(Li, sub_x);
            Li = MUL_DOUBLE_S(Li, lagrange_lookup);
        }

        SDOUBLE y_first = LOAD_DOUBLE(Y_FIRST_MID[i]);

        SDOUBLE eval_inner = MUL_DOUBLE_S(y_first, Li);
        result_q2 = ADD_DOUBLE_S(result_q2, eval_inner);
    }

    result_q2 = DIV_DOUBLE_S(one, result_q2);
    SIMD_TO_DOUBLE_VEC(&res[i], result_q2);
  }

}

```



# Calculation for non Corrected Taylor (Quadrant 1)

## Single Execution
OWN TIME OC: 98.426210999999995

OWN TIME CC: 98.002087624663702

```c
void tan_simd(double *input, double *res, size_t n, int prec) {
  int simd_doubles = 4;
  const SDOUBLE m_pi_2 = LOAD_DOUBLE(M_PI_2);

  int last_taylor_coeff = 13;
  int taylor_loop_iteration = 12;
  
  for (int i = 0; i < (int) n; i += 4) {
    SDOUBLE x   = LOAD_DOUBLE_VEC(&input[i]);
    SDOUBLE from_behind = SUB_DOUBLE_S(m_pi_2, x);
    
    /* ---- Calculation for first range ---- */
    const SDOUBLE x_square = MUL_DOUBLE_S(x, x);
    SDOUBLE result_q0 = LOAD_DOUBLE(TAYLOR_COEFF_TAN[last_taylor_coeff]);

    for (int j = taylor_loop_iteration; j >= 0; j-=1) {
      SDOUBLE coeff = LOAD_DOUBLE(TAYLOR_COEFF_TAN[j]);
      result_q0 = FMADD_PD(result_q0, x_square, coeff);
    }

    result_q0 = MUL_DOUBLE_S(result_q0, x);

    SIMD_TO_DOUBLE_VEC(&res[i], result_q0);

  }

}
```


## Double Execution

OWN TIME OC: 178.64845800000001

OWN TIME CC: 185.02178684241147

```c
void tan_simd(double *input, double *res, size_t n, int prec) {
  int simd_doubles = 4;
  const SDOUBLE m_pi_2 = LOAD_DOUBLE(M_PI_2);

  int last_taylor_coeff = 13;
  int taylor_loop_iteration = 12;
  
  for (int i = 0; i < (int) n; i += 4) {
    SDOUBLE x   = LOAD_DOUBLE_VEC(&input[i]);
    SDOUBLE from_behind = SUB_DOUBLE_S(m_pi_2, x);
    
    /* ---- Calculation for first range ---- */
    const SDOUBLE x_square = MUL_DOUBLE_S(x, x);
    SDOUBLE result_q0 = LOAD_DOUBLE(TAYLOR_COEFF_TAN[last_taylor_coeff]);

    for (int j = taylor_loop_iteration; j >= 0; j-=1) {
      SDOUBLE coeff = LOAD_DOUBLE(TAYLOR_COEFF_TAN[j]);
      result_q0 = FMADD_PD(result_q0, x_square, coeff);
    }

    result_q0 = MUL_DOUBLE_S(result_q0, x);

    SIMD_TO_DOUBLE_VEC(&res[i], result_q0);

  }

  for (int i = 0; i < (int) n; i += 4) {
    SDOUBLE x   = LOAD_DOUBLE_VEC(&res[i]);
    SDOUBLE from_behind = SUB_DOUBLE_S(m_pi_2, x);
    
    /* ---- Calculation for first range ---- */
    const SDOUBLE x_square = MUL_DOUBLE_S(x, x);
    SDOUBLE result_q0 = LOAD_DOUBLE(TAYLOR_COEFF_TAN[last_taylor_coeff]);

    for (int j = taylor_loop_iteration; j >= 0; j-=1) {
      SDOUBLE coeff = LOAD_DOUBLE(TAYLOR_COEFF_TAN[j]);
      result_q0 = FMADD_PD(result_q0, x_square, coeff);
    }

    result_q0 = MUL_DOUBLE_S(result_q0, x);

    SIMD_TO_DOUBLE_VEC(&res[i], result_q0);

  }

}
```


# Calculation for the Corrected Taylor (Quadrant 4)

## Single Execution 

OWN TIME OC: 105.78240599999999

OWN TIME CC: 104.86405872933535

```c 

void tan_simd(double *input, double *res, size_t n, int prec) {
  int simd_doubles = 4;

  const SDOUBLE m_pi_2 = LOAD_DOUBLE(M_PI_2);
  const SDOUBLE correction = LOAD_DOUBLE(CORRECTION);

  int last_taylor_coeff = 13;
  int taylor_loop_iteration = 12;

  const SDOUBLE one = LOAD_DOUBLE(1.0);
  
  for (int i = 0; i < (int) n; i += 4) {
    SDOUBLE x   = LOAD_DOUBLE_VEC(&input[i]);
    SDOUBLE from_behind = SUB_DOUBLE_S(m_pi_2, x);

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

    SIMD_TO_DOUBLE_VEC(&res[i], result_q3);

  }
}

```


## Double Execution

OWN TIME OC: 200.02065899999999

OWN TIME CC: 204.57104528577494


```c 
void tan_simd(double *input, double *res, size_t n, int prec) {
  int simd_doubles = 4;

  const SDOUBLE m_pi_2 = LOAD_DOUBLE(M_PI_2);
  const SDOUBLE correction = LOAD_DOUBLE(CORRECTION);

  int last_taylor_coeff = 13;
  int taylor_loop_iteration = 12;

  const SDOUBLE one = LOAD_DOUBLE(1.0);
  
  for (int i = 0; i < (int) n; i += 4) {
    SDOUBLE x   = LOAD_DOUBLE_VEC(&input[i]);
    SDOUBLE from_behind = SUB_DOUBLE_S(m_pi_2, x);

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

    SIMD_TO_DOUBLE_VEC(&res[i], result_q3);
  }

  for (int i = 0; i < (int) n; i += 4) {
    SDOUBLE x   = LOAD_DOUBLE_VEC(&res[i]);
    SDOUBLE from_behind = SUB_DOUBLE_S(m_pi_2, x);

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

    SIMD_TO_DOUBLE_VEC(&res[i], result_q3);

  }
}
```
