#pragma once
/*
    Classes for various dynamical systems. 
    These are systems of equations, using vector containers.
    Follow guidelines of Boost/numeric/odeint for defining systems in a class.


*/
#include <vector>

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