# CPlane
Mini project for solve for the phase plane of a system of differential equations. 2D, maybe 3D if it is cool enough. Using integration methods to solve a system of DEs, probably using 2 or 3 numerical methods.
Should also plot vector field, and isoclines,

Similar to PPLANE, but animated and shit because why not. And written in c++

# Animation example

![Small N](https://github.com/nicklayden/CPlane/blob/master/Peek 2017-10-23 18-48.giff "Nbody Simulation")


## Things to make: 
* [x] static 2d plot (3d if im bored) 
* [x] numerical integration methods that take in functions
* [ ] plot vector field 
* [ ] dynamic plot (animate orbits as function of time) 
	- draw individual orbits dynamically.
	- maybe with a collection of vectors that change length as function of time.
* [ ] draw isoclines somehow. (idk)
	- probably gonna need to do implicit plotting of these (dunno how to do that yet)

# Possible fun stuff to do:
* [ ] use GPU to calculate this stuff (marginal increase for double precision, maybe not worth) 

# External dependencies:
* [x] SFML (probably), maybe opengl or SDL 
* [ ] boost program options.
* [ ] maybe other boost stuff.
