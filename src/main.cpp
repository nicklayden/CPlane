/*
    SEE README.
    System of autonomous DE plane solver.


    // JUANS THESIS EQUATION:
        // Forward orbit from initial conditions.
        // xypts = adams_bashforth_4step(1 ,1 - 0.05*j,Nsteps,xprime,yprime,0,beta[i]);   
        // xypts2 = adams_bashforth_4step(-1 + 0.1*j, 0 + 0.1*j,Nsteps,xprime,yprime,0,beta[i]);   
        // Backward orbit from initial conditions.
        // xy3 = adams_bashforth_4step(-1 + 0.1*j, 0.1*j,Nsteps,xprime,yprime,0.3,beta[i], true);   
        

    TO DO:
        - Dump {x(t), y(t), t, beta} data to file for external plotting
            - 3D plot of bifurcation diagrams as function of beta.
        - Vector field? 


*/

#include <iostream>
#include "integrators.hpp"
#include "plotter.hpp"
#include <string>
#include <sstream>
#include <complex>


double yprime(double x, double y, double a, double b);
double xprime(double x, double y, double a, double b);
std::vector<std::vector<double> > adams_bashforth_4step(double x0, double y0, int Nsteps, double (*f)(double, double,double,double), double (*g)(double, double,double,double), double a, double b, bool timereverse=false);
std::string Numbertostring(std::string text, double num);
std::string Complextostring(std::string name, std::complex<double> comp);


int main()
{
    //______________________________________________________________________________//
    /*
        Screen Text */

    sf::Font font;
    sf::Text text, eq1,eq2,eq3;
    if (!font.loadFromFile("OpenSans-Regular.ttf")) {
        return 1;
    }

    text.setFont(font);
    text.setString("Test");
    text.setCharacterSize(16);
    text.setColor(sf::Color::White);
    text.setPosition(-320,-350);


    eq1.setFont(font);
    eq1.setString("Test");
    eq1.setCharacterSize(16);
    eq1.setColor(sf::Color::White);
    eq1.setPosition(-100,-350);

    //______________________________________________________________________________//



    // containers and Number of steps to run the solver for.
    std::vector<std::vector<double> > xypts,xypts2,xy3, xypts3,xypts4;
    std::vector<double> xeqpoints(3), yeqpoints(3);
    xeqpoints[0] = 0; yeqpoints[0] = 0; // point (0,0)
    int Nsteps = 50000; // num time steps to integrate.
    std::complex<double> xeq; // x component of the eq points.
    std::complex<double> root_xeq; // complex sqrt from <complex> header


    // Range of alpha (damping) and beta (magnetic 'strength') values.
    std::vector<double> beta(100),alpha(100);
    for (size_t i = 0; i < beta.size(); i++) {
        beta[i] = 0.7 + 0.01*i;
        alpha[i] = 0.1 + 0.005*i;
    }

    // Three Plots to plot the data.
    Plot plt("CPlane",2,2);
    Plot plt2("x(t)",50,5,600,600);
    Plot plt3("y(t)",50,5,600,600);
    plt2.plotView.setCenter(50./2., 0);
    plt3.plotView.setCenter(50./2., 0);
    Plot plt4("bifurcation diagram",2,2,300,300);
    plt4.plotView.setCenter(1,0);

    // plt.mainwindow.setFramerateLimit(10);
    while(plt.mainwindow.isOpen()) {
        // SEG FAULT CORE DUMPED ON RANDOM CLOSES OF THE WINDOWS ???????? NOT CONSISTENTLY HAPPENING!!
        for (int i = 0; i < beta.size(); i++) {

            std::complex<double> xeq = 1 - 1/beta[i]; // x component of the eq points.
            std::complex<double> root_xeq = std::sqrt(xeq); // complex sqrt from <complex> header
            xeqpoints[1] = root_xeq.real(); // get the real part of the eq point.
            yeqpoints[1] = 0;
            xeqpoints[2] = -root_xeq.real();
            yeqpoints[2] = 0;
            for (int j = 0; j < 5; j++) {
                // MAGNETO ELASTIC BEAM
                xypts2 = adams_bashforth_4step(-j,-j, Nsteps, xprime, yprime,0.1,beta[i]);
                xypts = adams_bashforth_4step(j,j, Nsteps, xprime, yprime,0.1,beta[i]);


                plt.plot(xypts[0],xypts[1],sf::Color::Green); // phase plane
                plt.plot(xypts2[0],xypts2[1],sf::Color::Green); // phase plane
                plt2.plot(xypts[2],xypts[0],sf::Color::Blue); // X-T plane
                plt3.plot(xypts[2],xypts[1],sf::Color::Red);  // Y-T plane
            }

            // plt.set_xlabel(std::string("nips"), font);
            plt4.scatter(beta[i],xeqpoints[1],sf::Color::Red);
            plt4.scatter(beta[i],xeqpoints[2],sf::Color::Red);
            plt4.mainwindow.display();
            plt2.mainwindow.display();
            plt2.mainwindow.clear(sf::Color::Black);
            plt3.mainwindow.display();
            plt3.mainwindow.clear(sf::Color::Black);
            plt.scatter(xeqpoints,yeqpoints, sf::Color::Red);
            plt.mainwindow.setView(plt.axesView);
            text.setString(Numbertostring("Beta =",beta[i]));
            eq1.setString(Complextostring("x_eq",root_xeq));
            plt.mainwindow.draw(text);
            plt.mainwindow.draw(eq1);
            plt.mainwindow.display();
            plt.mainwindow.clear(sf::Color::Black);
        }
        plt4.mainwindow.clear();
        plt.EventLoop();  
        plt2.EventLoop();  
        plt3.EventLoop();  
        plt4.EventLoop();  
    }
    // plt.show();


}



/*
    Adams bashforth 4 step method. Works ok enough. Will switch to RK4 when I have time.
*/
std::vector<std::vector<double> > adams_bashforth_4step(double x0, double y0, int Nsteps,
                                                        double (*f)(double, double,double,double), 
                                                        double (*g)(double, double,double,double),
                                                        double a, double b,bool timereverse)
{
    double xn,xn1,xn2,xn3;
    double tn,tn1,tn2,tn3;
    double yn,yn1,yn2,yn3;
    double h = 0.001;
    xn = x0;
    yn = y0;
    tn = 0.0;
    std::vector<std::vector<double> > vals;
    std::vector<double> xvalues(Nsteps), yvalues(Nsteps), timevals(Nsteps);


    // Initializing the AB4 method with euler steps:
    if (!timereverse) {
        tn1 = tn + h;
        yn1 = yn + h*f(xn,yn,a,b);
        xn1 = xn + h*g(xn,yn,a,b);

        tn2 = tn1 + h;
        yn2 = yn1 + h*f(xn1,yn1,a,b);
        xn2 = xn1 + h*g(xn1,yn1,a,b);

        tn3 = tn2 + h;
        yn3 = yn2 + h*f(xn2,yn2,a,b);
        xn3 = xn2 + h*g(xn2,yn2,a,b);
    } else {
        tn1 = tn - h;
        yn1 = yn - h*f(xn,yn,a,b);
        xn1 = xn - h*g(xn,yn,a,b);
    
        tn2 = tn1 - h;
        yn2 = yn1 - h*f(xn1,yn1,a,b);
        xn2 = xn1 - h*g(xn1,yn1,a,b);
    
        tn3 = tn2 - h;
        yn3 = yn2 - h*f(xn2,yn2,a,b);
        xn3 = xn2 - h*g(xn2,yn2,a,b);
    
    }
    // Solving the equations.
    if (!timereverse) {
        for (size_t i = 0; i < Nsteps; i++) {

            yn3 += (h/24.0)*( 55.0*(g(xn3,yn3,a,b)) - 59.0*(g(xn2,yn2,a,b)) + 37.0*(g(xn1,yn1,a,b)) - 9.0*(g(xn,yn,a,b)));
            xn3 += (h/24.0)*( 55.0*(f(xn3,yn3,a,b)) - 59.0*(f(xn2,yn2,a,b)) + 37.0*(f(xn1,yn1,a,b)) - 9.0*(f(xn,yn,a,b)));
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
            timevals[i] = tn3;
            xvalues[i] = xn3;
            yvalues[i] = yn3;
        }
    } else {
        for (size_t i = 0; i < Nsteps; i++) {

            yn3 -= (h/24.0)*( 55.0*(g(xn3,yn3,a,b)) - 59.0*(g(xn2,yn2,a,b)) + 37.0*(g(xn1,yn1,a,b)) - 9.0*(g(xn,yn,a,b)));
            xn3 -= (h/24.0)*( 55.0*(f(xn3,yn3,a,b)) - 59.0*(f(xn2,yn2,a,b)) + 37.0*(f(xn1,yn1,a,b)) - 9.0*(f(xn,yn,a,b)));
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
            timevals[i] = tn3;
            xvalues[i] = xn3;
            yvalues[i] = yn3;
        }
    }
    vals.push_back(xvalues); vals.push_back(yvalues), vals.push_back(timevals);
    return vals;
}

double xprime(double x, double y, double a, double b)
{
    // juan thesis eqn : x*(2*x*x - y*y - 2) - (sqrt(6)/2)*b*y*y;
  return y;
}

double yprime(double x, double y, double a, double b)
{
    // juan thesis eqn : y*(2*x*x - y*y + 1 + (sqrt(6)/2)*b*x );
  return -2*a*y - x + b*x*(1-x*x);
}


std::string Numbertostring(std::string text, double num)
{
    std::ostringstream ss;
    ss << text << " " << num << " ";
    return ss.str();
}

std::string Complextostring(std::string name,std::complex<double> comp)
{
    std::ostringstream ss;
    ss << name << " = " << comp.real() << " +i " << comp.imag();
    return ss.str();
}