#include<stdio.h>
#include<assert.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>

int main () {
    gsl_rng * myrng = gsl_rng_alloc(gsl_rng_taus);
    double a;

    for(int i = 0; i < 100; i++) {
        a = 1.0 + gsl_ran_gaussian_tail(myrng, -1.0, 5);
        printf("%1.2f\t", a);
        assert(a > 0.0);
    }
    printf("\n");

    return 0;
}
