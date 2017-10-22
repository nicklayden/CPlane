/*
    Header file for the base abstract class integrator, and numerical methods.
*/

#include <iostream>
#include <string>


/*
    Pure virtual function step() must be overridden by all inheriting classes.
*/
class baseintegrator
{
    public:
        virtual void step() = 0;
};


/*
    Various templates for numerical integration methods.
*/
class rk4: public baseintegrator
{
    public:
        void step() { std::cout << "Stepping with rk4..." << std::endl;}
};

class euler: public baseintegrator
{
    public:
        void step() { std::cout << "Stepping with euler..." << std::endl;}
};

class ab4: public baseintegrator
{
    public:
        void step() { std::cout << "Stepping with adams bashforth 4 step..." << std::endl;}
};

/*
    Generic integrator class for instantiating any numerical method
*/
class integrator
{
    public:
        integrator(baseintegrator* ptr_method) :ptr_method(ptr_method) {}
        void step() { ptr_method->step();}
    private:
        baseintegrator* ptr_method;
};


/*
    Factory function for creating the instances
*/
integrator pickmethod(std::string s)
{
    if (s == "rk4")
    {
        std::cout << "Creating RK4 integrator.\n";
        return integrator(new rk4);
    }
    else if (s == "euler")
    {
        std::cout << "Creating Euler integrator.\n";
        return integrator(new euler);
    }
    else if (s == "ab4")
    {
        std::cout << "Creating AB4 integrator.\n";
        return integrator(new ab4);
    }
    else 
    {
        std::cerr << "ERROR: METHOD NOT FOUND.\n DEFAULTING TO EULER.\n";
        return integrator(new euler);
    }
}