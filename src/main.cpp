/*
    SEE README.
    System of autonomous DE plane solver.

*/

#include <iostream>
#include "integrators.hpp"
#include "plotter.hpp"


double yprime(double x, double y);
double xprime(double x, double y);
std::vector<std::vector<double> > adams_bashforth_4step(double x0, double y0, int Nsteps, double (*f)(double, double), double (*g)(double, double));


int main()
{
    // gonna have to rework the integrators. it sucks right now.
    // std::string smethod = "ab4";
    // integrator amethod = pickmethod(smethod);

    std::vector<std::vector<double> > xypts,xypts2,xy3,xy4;
    std::vector<double> x,y,x2,y2,x3,y3,x4,y4;
    int Nsteps = 100000;


    xypts = adams_bashforth_4step(1,1,Nsteps,xprime,yprime);
    xypts2 =  adams_bashforth_4step(-1,-1,Nsteps,xprime,yprime);
    xy3 = adams_bashforth_4step(0.5,1.4,Nsteps,xprime,yprime);
    xy4 = adams_bashforth_4step(-0.1,-1.2,Nsteps,xprime,yprime);

    for (size_t i = 0; i < Nsteps; i++) {
        x.push_back(xypts[i][0]);
        y.push_back(xypts[i][1]);
        x2.push_back(xypts2[i][0]);   
        y2.push_back(xypts2[i][1]);
        x3.push_back(xy3[i][0]);
        y3.push_back(xy3[i][1]);
        x4.push_back(xy4[i][0]);    
        y4.push_back(xy4[i][1]);
    }

    Plot plt("TEST");
    plt.plot(x,y);
    plt.plot(x2,y2);
    plt.plot(x3,y3, sf::Color(0,0,0,255));
    plt.plot(x4,y4, sf::Color::Red);
    plt.show();


}




std::vector<std::vector<double> > adams_bashforth_4step(double x0, double y0, int Nsteps,
                             double (*f)(double, double), double (*g)(double, double))
{
    double xn,xn1,xn2,xn3;
    double tn,tn1,tn2,tn3;
    double yn,yn1,yn2,yn3;
    double h = 0.001;
    xn = x0;
    yn = y0;
    tn = 0.0;
    std::vector<std::vector<double> > vals;


    // Initializing the AB4 method with euler steps:
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

double xprime(double x, double y)
{
  return y;
}

double yprime(double x, double y)
{
  return -2*0.1*y - x + 1.2*x*(1 - x*x);
}
