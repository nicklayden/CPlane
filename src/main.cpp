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
#include <boost/numeric/odeint/integrate/integrate_const.hpp>
#include <boost/numeric/odeint.hpp>

namespace ni = boost::numeric::odeint;


std::string Numbertostring(std::string text, double num);
std::string Complextostring(std::string name, std::complex<double> comp);
std::vector<std::vector<double> > transpose_copy(std::vector<std::vector<double> > data);
void eq_pt(std::vector<double>& x, std::vector<double>& y, double lambda);
void bif_pt1(std::vector<double>& x, std::vector<double>& y);
void bif_pt2(std::vector<double>& x, std::vector<double>& y);


struct push_back_state_and_time 
{
    std::vector<std::vector<double> >& m_states;
    std::vector<double>& m_times;

    push_back_state_and_time(std::vector<std::vector<double> >& states, std::vector<double>& times)
    :m_states(states), m_times(times) { }

    void operator()(const std::vector<double>& x, double t)
    {
        m_states.push_back(x);
        m_times.push_back(t);
    }


};

class sean_problem
{
public:
    double m_par;
    sean_problem(double par): m_par(par) {}
    // exponential potential with constant epsilon background

    void operator() (const std::vector<double>& x, std::vector<double>& dxdt, const double /* t */) {
        dxdt[0] = x[0]*(-2 + 2*x[0]*x[0] - x[1]*x[1] - x[2]*x[2]) + (sqrt(6.)/2.)*m_par*x[2]*x[2];
        dxdt[1] = x[1]*(1 + 2*x[0]*x[0] - x[1]*x[1] - x[2]*x[2]);
        dxdt[2] = x[2]*((sqrt(6.)/2.)*m_par*x[0] + 3*x[0]*x[0] + 1 - x[0]*x[0] - x[1]*x[1] - x[2]*x[2] );
    }
};


class dynamical_system
{
public:
    double  m_parameter;

    dynamical_system(double parameter) :m_parameter(parameter) {}

    void operator() (const std::vector<double>& x, std::vector<double>& dxdt, const double /*t*/){
        dxdt[0] = -2*x[0] - (x[0]*x[0] + x[1]*x[1]);
        dxdt[1] = -(1 - x[0])*x[1];
    }

};

class quintessence
{
public:
    double  lambda;

    quintessence(double parameter) :lambda(parameter) {}

    void operator() (const std::vector<double>& x, std::vector<double>& dxdt, const double /*t*/){
        dxdt[0] = x[0]*(2*x[0]*x[0] - x[1]*x[1] - 2) + (sqrt(6.)/2.)*lambda*x[1]*x[1] ;
        dxdt[1] = x[1]*(-lambda*(sqrt(6.)/2.)*x[0] + 2*x[0]*x[0] - x[1]*x[1] +1 );
    }

};

int main()
{
    //______________________________________________________________________________//
    /*
        Screen Text */

    // sf::Font font;
    // sf::Text text, eq1,eq2,eq3;
    // if (!font.loadFromFile("OpenSans-Regular.ttf")) {
    //     return 1;
    // }

    // text.setFont(font);
    // text.setString("Test");
    // text.setCharacterSize(16);
    // text.setColor(sf::Color::White);
    // text.setPosition(-320,-350);


    // eq1.setFont(font);
    // eq1.setString("Test");
    // eq1.setCharacterSize(16);
    // eq1.setColor(sf::Color::White);
    // eq1.setPosition(-100,-350);

    //______________________________________________________________________________//
    //______________________________________________________________________________//
    // BOOST NUMERIC.ODEINT METHOD OF INTEGRATION.

    std::vector<double> x(2);
    std::vector<std::vector<double> > x_vec, test_transpose,y2,test_transpose2;
    std::vector<double> times,xstate,ystate, t2;

    std::vector<double> circle, cx, cy;
    // std::vector<double> range;
    for (int i = 0; i < 200; i++) {
        float angle = M_PI*i/200;
        cx.push_back(cos(angle));
        cy.push_back(sin(angle));
    }


    // Integration method
    ni::runge_kutta_dopri5<std::vector<double> > stepper;

    std::vector<double> x_eq, y_eq;
    std::vector<double> x_bif, y_bif;

    // bifurcation points
    bif_pt1(x_bif,y_bif);
    bif_pt2(x_bif,y_bif);


    Plot plt_animate("Animation test for odeint",4,2, 800,400);
    // plt_animate.mainwindow.setFramerateLimit(15);
    plt_animate.plotView.setCenter(sf::Vector2f(0,0.5));

    // Plot slices
    // Plot xy("XY plane",4,2,800,400);
    // Plot zx("ZX plane",4,2,800,400);
    // Plot yz("YZ plane",4,2,800,400);
    // xy.plotView.setCenter(sf::Vector2f(0,0.5));
    // zx.plotView.setCenter(sf::Vector2f(0,0.5));
    // yz.plotView.setCenter(sf::Vector2f(0,0.5));



    while(plt_animate.mainwindow.isOpen()) {
        for (size_t i = 0; i < 250; i++)
        {
            double lambda = 0.01*i;
            for (int j = 0; j < 20; j++)
            {   
                quintessence system2(lambda); 
                
                x[0] = -sqrt(3)/2. + sqrt(3)*j/20; x[1] = 0.5; //+ 0.005*j; //x[2] = 0.2 + 0.005*j;
                ni::integrate_const(stepper,system2, x, 0.,5.,0.001,push_back_state_and_time(y2,t2)); // forward solution

                x[0] = -sqrt(3)/2 + sqrt(3)*j/20; x[1] = 0.5; //+ 0.005*j; //x[2] = 0.2 + 0.005*j;// reset initial conditions for back solution.
                ni::integrate_adaptive(stepper,system2, x, 5.,0.,-0.001, push_back_state_and_time(x_vec,times)); // backward solution

                // transposing Nx2 vector into 2XN
                test_transpose2 = transpose_copy(y2);
                test_transpose = transpose_copy(x_vec);

                // plot solution curves in phase plane
                plt_animate.plot(test_transpose2[0], test_transpose2[1], sf::Color::Red);
                plt_animate.plot(test_transpose[0], test_transpose[1], sf::Color::Red);

                // xy.plot(test_transpose2[0], test_transpose2[1], sf::Color::Red);
                // xy.plot(test_transpose[0], test_transpose[1], sf::Color::Red);

                // zx.plot(test_transpose2[0], test_transpose2[2], sf::Color::Red);
                // zx.plot(test_transpose[0], test_transpose[2], sf::Color::Red);
                
                // yz.plot(test_transpose2[2], test_transpose2[1], sf::Color::Red);
                // yz.plot(test_transpose[2], test_transpose[1], sf::Color::Red);


                // Clear vectors for next solution curve.
                y2.clear(); t2.clear(); times.clear(); x_vec.clear(); 
            }

            // // plot eq points (non trivial)
            eq_pt(x_eq,y_eq,lambda);
            plt_animate.plot(cx,cy);
            plt_animate.scatter(x_eq, y_eq, sf::Color::Cyan);
            plt_animate.scatter(x_bif, y_bif, sf::Color::Green);
            plt_animate.show();
            // xy.show();
            // zx.show();
            // yz.show();
            x_eq.clear(); y_eq.clear();
            x_eq.push_back(0);x_eq.push_back(1);x_eq.push_back(-1); 
            y_eq.push_back(0);y_eq.push_back(0);y_eq.push_back(0);
        
        }
        plt_animate.EventLoop();
        // xy.EventLoop();
        // zx.EventLoop();
        // yz.EventLoop();
    } // end while

} // end main



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

std::vector<std::vector<double> > transpose_copy(std::vector<std::vector<double> > data) {
    
    // Assuming all inner vectors are same length. we can allocate space in advance.
    std::vector<std::vector<double> > result(data[0].size(), std::vector<double>(data.size()));
    // std::cout << "i= " << data[0].size() << std::endl;
    // std::cout << "j= " << data.size() << std::endl;
    for (int i = 0; i < data[0].size(); i++) {
        for (int j = 0; j < data.size(); j++) {
            result[i][j] = data[j][i];
        }
    }
    return result;

}

void bif_pt2(std::vector<double>& x, std::vector<double>& y) {
    double x2,y2; // bif points for eq_point x2
    double lambda = sqrt(8./3.);
    // bif at sqrt(8/3)
    x2 = 2./(lambda*sqrt(6.));
    y2 = 2./(lambda*sqrt(3.));

    x.push_back(x2);
    y.push_back(y2);
}

void bif_pt1(std::vector<double>& x, std::vector<double>& y) {
    double x1,y1;
    double lambda=sqrt(2.);

    x1 = lambda/sqrt(6.);
    y1 = sqrt(1 - (lambda*lambda)/6.);
    x.push_back(x1);
    y.push_back(y1);
}

void eq_pt(std::vector<double>& x, std::vector<double>& y, double lambda) {
    // Returns the EQ PTS for the quintessence model with exponential curve (the basic one)
    //  first: x = lambda/sqrt(6),      y = sqrt(1 - lambda^2/6)
    // second: x = 2/(sqrt(6)lambda),   y = 2/(sqrt(3)lambda)

    double x1,x2,y1,y2;
    
    // if(lambda*lambda < 6){
        x1 = lambda/sqrt(6);
        y1 = sqrt(1 - lambda*lambda/6);
        x.push_back(x1);
        y.push_back(y1);
    // }
    if(lambda > sqrt(2)) {
    x2 = 2./(sqrt(6.)*lambda);
    y2 = 2./(sqrt(3.)*lambda);

    x.push_back(x2);
    y.push_back(y2);
    }
}