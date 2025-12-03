

# Testing in full Algorithm



This has times of 

2121
2094
2186

```c 
    /* obtaining bool vectors for the each quadrant */
    const SDOUBLE q_sub_2 = SUB_DOUBLE_S(quadrant, two);
    const SDOUBLE q_sub_1 = SUB_DOUBLE_S(quadrant, one);

    const SDOUBLE abs_q_sub_2 = ABS_PD(q_sub_2);
    const SDOUBLE abs_q_sub_1 = ABS_PD(q_sub_1);

    const SDOUBLE in_q0_0 = MUL_DOUBLE_S(abs_q_sub_2, half);
    const SDOUBLE in_q3_0 = MUL_DOUBLE_S(abs_q_sub_1, half);

    const SDOUBLE in_q1_0 = SUB_DOUBLE_S(abs_q_sub_1, two);
    const SDOUBLE in_q2_0 = SUB_DOUBLE_S(abs_q_sub_2, two);

    const SDOUBLE in_q0 = FLOOR_DOUBLE_S(in_q0_0);
    const SDOUBLE in_q3 = FLOOR_DOUBLE_S(in_q3_0);

    const SDOUBLE in_q1_1 = ABS_PD(in_q1_0);
    const SDOUBLE in_q2_1 = ABS_PD(in_q2_0);

    const SDOUBLE in_q1_2 = MUL_DOUBLE_S(in_q1_1, half);
    const SDOUBLE in_q2_2 = MUL_DOUBLE_S(in_q2_1, half);

    const SDOUBLE in_q1   = FLOOR_DOUBLE_S(in_q1_2);
    const SDOUBLE in_q2   = FLOOR_DOUBLE_S(in_q2_2);
```



2095
2076
2173

```c 
    /* obtaining bool vectors for the each quadrant */
    // 1 if quadrant == 0 else 0
    SDOUBLE in_q0 = SUB_DOUBLE_S(quadrant, two);
    in_q0 = ABS_PD(in_q0);
    in_q0 = MUL_DOUBLE_S(in_q0, half);
    in_q0 = FLOOR_DOUBLE_S(in_q0);

    // 1 if quadrant == 1 else 0
    SDOUBLE in_q1 = SUB_DOUBLE_S(quadrant, one);
    in_q1 = ABS_PD(in_q1);
    in_q1 = SUB_DOUBLE_S(in_q1, two);
    in_q1 = ABS_PD(in_q1);
    in_q1 = MUL_DOUBLE_S(in_q1, half);
    in_q1 = FLOOR_DOUBLE_S(in_q1);

    // 1 if quadrant == 2 else 0
    SDOUBLE in_q2 = SUB_DOUBLE_S(quadrant, two);
    in_q2 = ABS_PD(in_q2);
    in_q2 = SUB_DOUBLE_S(in_q2, two);
    in_q2 = ABS_PD(in_q2);
    in_q2 = MUL_DOUBLE_S(in_q2, half);
    in_q2 = FLOOR_DOUBLE_S(in_q2);


    // 1 if quadrant == 3 else 0
    SDOUBLE in_q3 = SUB_DOUBLE_S(quadrant, one);
    in_q3 = ABS_PD(in_q3);
    in_q3 = MUL_DOUBLE_S(in_q3, half);
    in_q3 = FLOOR_DOUBLE_S(in_q3);
```



## other test in algorithm
