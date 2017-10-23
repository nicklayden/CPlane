/*
    SEE README.
    System of autonomous DE plane solver.

*/

#include <iostream>
#include "integrators.hpp"
#include "plotter.hpp"


double yprime(double x, double y);
double xprime(double x, double y);
std::vector<std::vector<double> > adams_bashforth_4step(double x0, double y0, int Nsteps, double (*f)(double, double), double (*g)(double, double), bool timereverse=false);


int main()
{
    // gonna have to rework the integrators. it sucks right now.
    // std::string smethod = "ab4";
    // integrator amethod = pickmethod(smethod);

    std::vector<std::vector<double> > xypts,xypts2,xy3,xy4;
    std::vector<double> x,y,x2,y2,x3,y3,x4,y4;
    int Nsteps = 100000;


    // xypts = adams_bashforth_4step(1,1,Nsteps,xprime,yprime, true);
    xypts2 =  adams_bashforth_4step(1,1,Nsteps,xprime,yprime);
    // xy3 = adams_bashforth_4step(0.5,1.4,Nsteps,xprime,yprime);
    // xy4 = adams_bashforth_4step(-0.1,-1.2,Nsteps,xprime,yprime);


    // x = xypts[0];
    // y = xypts[1];
    

    std::cout << x.size() << std::endl;
    std::cout << y.size() << std::endl;

    Plot plt("TEST",2,2);
    // plt.setxlim(-.1,.1);
    // plt.setylim(-.1,.1);
    plt.EventLoop();
    // plt.plot(x,y);
    plt.plot(xypts2[0],xypts2[1]);
    // plt.plot(x3,y3, sf::Color::Red);
    // plt.plot(x4,y4, sf::Color::Red);
    plt.show();


}




std::vector<std::vector<double> > adams_bashforth_4step(double x0, double y0, int Nsteps,
                             double (*f)(double, double), double (*g)(double, double),
                             bool timereverse)
{
    double xn,xn1,xn2,xn3;
    double tn,tn1,tn2,tn3;
    double yn,yn1,yn2,yn3;
    double h = 0.001;
    xn = x0;
    yn = y0;
    tn = 0.0;
    std::vector<std::vector<double> > vals;
    std::vector<double> xvalues(Nsteps), yvalues(Nsteps);


    // Initializing the AB4 method with euler steps:
    if (!timereverse) {
        tn1 = tn + h;
        yn1 = yn + h*f(xn,yn);
        xn1 = xn + h*g(xn,yn);

        tn2 = tn1 + h;
        yn2 = yn1 + h*f(xn1,yn1);
        xn2 = xn1 + h*g(xn1,yn1);

        tn3 = tn2 + h;
        yn3 = yn2 + h*f(xn2,yn2);
        xn3 = xn2 + h*g(xn2,yn2);
    } else {
        tn1 = tn - h;
        yn1 = yn - h*f(xn,yn);
        xn1 = xn - h*g(xn,yn);
    
        tn2 = tn1 - h;
        yn2 = yn1 - h*f(xn1,yn1);
        xn2 = xn1 - h*g(xn1,yn1);
    
        tn3 = tn2 - h;
        yn3 = yn2 - h*f(xn2,yn2);
        xn3 = xn2 - h*g(xn2,yn2);
    
    }
    // Solving the equations.
    if (!timereverse) {
        for (size_t i = 0; i < Nsteps; i++) {

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

            xvalues[i] = xn3;
            yvalues[i] = yn3;
        }
    } else {
        for (size_t i = 0; i < Nsteps; i++) {

            yn3 -= (h/24.0)*( 55.0*(g(xn3,yn3)) - 59.0*(g(xn2,yn2)) + 37.0*(g(xn1,yn1)) - 9.0*(g(xn,yn)));
            xn3 -= (h/24.0)*( 55.0*(f(xn3,yn3)) - 59.0*(f(xn2,yn2)) + 37.0*(f(xn1,yn1)) - 9.0*(f(xn,yn)));
            yn = yn1;
            tn = tn1;
            xn = xn1;

            yn1 = yn2;
            tn1 = tn2;
            xn1 = xn2;

            yn2 = yn3;
            tn2 = tn3;
            xn2 = xn3;

            tn3 -= h;
            xvalues[i] = xn3;
            yvalues[i] = yn3;
        }
    }
    vals.push_back(xvalues); vals.push_back(yvalues);
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
