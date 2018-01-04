# CPlane
Mini project for solving for and plotting the phase plane of an autonomous system of differential equations. Plot is 2D, maybe 3D if SFML is swapped out for another graphics solution. Using integration methods to solve a system of DEs, probably using 2 or 3 numerical methods. Coupled with boost/odeint, this program can also solve systems with more than 2 dimensions, guidelines for writing a class or a function for use in odeint integration functions are given in the boost/numeric/odeint documentation.
So far, the program needs to be recompiled for any new systems or parameters entered in. Might change later depending on time and motivation.

The idea to do this was PPLANE, but animated and shit because why not. And written in c++

# Animation example

![Changing magnetic strength](https://github.com/nicklayden/CPlane/blob/master/Peek%202017-10-23%2018-48.gif "Magneto Elastic Beam")

Magneto elastic beam problem. With damping term a=0.1. The magnetic force term is beta=[0.9,1.9]


## Things to make: 
* [x] 2D plot sublibrary (3d if im bored) 
* [x] numerical integration methods that take in functions (currently uses boost/odeint)
* [ ] plot vector field 
* [ ] dynamic plot (animate orbits as function of time) 
	- draw individual orbits dynamically.
	- maybe with a collection of vectors that change length as function of time.
* [ ] draw isoclines somehow. (idk)
	- probably gonna need to do implicit plotting of these

# Possible fun stuff to do:
* [ ] use GPU to calculate this stuff (marginal increase for double precision, maybe not worth) 

# External dependencies:
* [x] SFML 
* [ ] boost program options. for command line control.
* [x] boost odeint for integration.
