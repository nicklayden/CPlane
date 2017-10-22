#include "integrators.hpp"
#include <vector>

ab4::ab4()
{
    std::cout << "Created ab4 integrator method." << std::endl;
}

std::vector<std::vector<double> > ab4::integrate(double x0, double y0, int Nsteps, double (*f)(double, double), double (*g)(double, double))
{
    double xn,xn1,xn2,xn3;
    double tn,tn1,tn2,tn3;
    double yn,yn1,yn2,yn3;
    xn = x0;
    yn = y0;
    tn = 0.0;
    std::vector<std::vector<double> > vals;


    // Initializing the AB4 method with euler steps:
    // SHOULD be using rk4 to start this... maybe later.
    tn1 = tn + h;
    yn1 = yn + h*f(xn,yn);
    xn1 = xn + h*g(xn,yn);

    tn2 = tn1 + h;
    yn2 = yn1 + h*f(xn1,yn1);
    xn2 = xn1 + h*g(xn1,yn1);

    tn3 = tn2 + h;
    yn3 = yn2 + h*f(xn2,yn2);
    xn3 = xn2 + h*g(xn2,yn2);

    // Solving the equations.
    // Maybe change the way that values are stored in vectors.
    // Perhaps a single column vector for each dimension is a better idea.
    // right now its a mess of 2d vectors stacked together.
    for (size_t i = 0; i < Nsteps; i++) {
        std::vector<double> vector2d;

        yn3 += (h/24.0)*( 55.0*(g(xn3,yn3)) - 59.0*(g(xn2,yn2)) + 37.0*(g(xn1,yn1)) - 9.0*(g(xn,yn)));
        xn3 += (h/24.0)*( 55.0*(f(xn3,yn3)) - 59.0*(f(xn2,yn2)) + 37.0*(f(xn1,yn1)) - 9.0*(f(xn,yn)));
        yn = yn1;
        tn = tn1;
        xn = xn1;

        yn1 = yn2;
        tn1 = tn2;
        xn1 = xn2;

        yn2 = yn3;
        tn2 = tn3;
        xn2 = xn3;

        tn3 += h;
        vector2d.push_back(xn3);
        vector2d.push_back(yn3);
        vals.push_back(vector2d);
    }
    return vals;
}