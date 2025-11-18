
# Time Needed for the Lagrange Polyromial

## Single Run
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


## Double Run


OWN TIME CC: 11908.752674338131
TIME OC:     11935.889466000001

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
