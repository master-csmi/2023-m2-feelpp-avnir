#include <math.h>

double wavelet(double t)
{
    double fc = 10000;
    return sin(2*M_PI*fc*t)*exp(-5*pow(fc*t-2,2));
}