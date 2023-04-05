//
// Created by Dan Wilkins on 3/18/21.
//
#include <iostream>
using namespace std;

#include <gsl/gsl_integration.h>

double integrand(double x, void* params)
{
    double po = *(double *) params;
    return pow(x, po);
}

int main()
{
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);

    double result, error;
    double expected = 1./3;
    double po = 2.0;

    gsl_function F;
    F.function = &integrand;
    F.params = &po;

    gsl_integration_qags (&F, 0, 1, 0, 1e-7, 1000, w, &result, &error);

    cout << "result          = " << result << endl;
    cout << "exact result    = " << expected << endl;
    cout << "estimated error = " << error << endl;
    cout << "actual error    = " << result - expected << endl;
    cout << "intervals       = " << w->size << endl;

    gsl_integration_workspace_free(w);

    return 0;
}
