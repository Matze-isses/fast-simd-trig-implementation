#include <stdio.h>
#include <flint/arb.h>
#include <flint/arf.h>

static double arb_to_double_correct(const arb_t x)
{
    slong prec = 128;
    arb_t t;
    double d;

    arb_init(t);

    for (;;)
    {
        arb_set_round(t, x, prec);

        if (arb_can_round_arf(t, 53, ARF_RND_NEAR))
        {
            d = arf_get_d(arb_midref(t), ARF_RND_NEAR);
            break;
        }

        prec *= 2;
    }

    arb_clear(t);
    return d;
}

static void split_multiple_of_pi_arb(const char *factor_str, slong work_prec)
{
    arb_t factor;
    arb_t pi;
    arb_t x;
    arb_t sum;
    arb_t rest;
    arb_t tmp;

    double a;
    double b;
    double c;

    arb_init(factor);
    arb_init(pi);
    arb_init(x);
    arb_init(sum);
    arb_init(rest);
    arb_init(tmp);

    if (arb_set_str(factor, factor_str, work_prec) != 0)
    {
        fprintf(stderr, "Could not parse factor: %s\n", factor_str);
    arb_clear(factor);
    arb_clear(pi);
    arb_clear(x);
    arb_clear(sum);
    arb_clear(rest);
    arb_clear(tmp);
    }

    arb_const_pi(pi, work_prec);
    arb_mul(x, factor, pi, work_prec);

    a = arb_to_double_correct(x);

    arb_set_d(sum, a);
    arb_sub(rest, x, sum, work_prec);

    b = arb_to_double_correct(rest);

    arb_set_d(tmp, b);
    arb_add(sum, sum, tmp, work_prec);
    arb_sub(rest, x, sum, work_prec);

    c = arb_to_double_correct(rest);

    arb_set_d(tmp, c);
    arb_add(sum, sum, tmp, work_prec);
    arb_sub(rest, x, sum, work_prec);

    printf("factor                 = %s\n", factor_str);

    printf("x (arb)                = ");
    arb_printn(x, 60, 0);
    printf("\n");

    printf("a (double)             = %a   %.17e\n", a, a);
    printf("b (double)             = %a   %.17e\n", b, b);
    printf("c (double)             = %a   %.17e\n", c, c);

    printf("a+b+c (arb)            = ");
    arb_printn(sum, 60, 0);
    printf("\n");

    printf("remaining error        = ");
    arb_printn(rest, 60, 0);
    printf("\n");
}


int main(void)
{
    split_multiple_of_pi_arb("1.0", 512);
    return 0;
}

/* Compilation by: 

gcc get_long_double_split.c -o test -lflint -lgmp -lm && ./test 

*/
