# CPlane
Mini project for solve for the phase plane of a system of differential equations. 2D, maybe 3D if it is cool enough. Using integration methods to solve a system of DEs, probably using 2 or 3 numerical methods.
Should also plot vector field, and isoclines,

Similar to PPLANE, but animated and shit because why not. And written in c++



## Things to make: 
	* [x] static 2d plot (3d if im bored) # [checkbox:checked]
	* [x] numerical integration methods that take in functions # [checkbox:checked]
	* [ ] plot vector field # [checkbox:unchecked]
	* [ ] dynamic plot (animate orbits as function of time) # [checkbox:unchecked]
		- draw individual orbits dynamically.
		- maybe with a collection of vectors that change length as function of time.
	* [ ] draw isoclines somehow. (idk) # [checkbox:unchecked]
		- probably gonna need to do implicit plotting of these (dunno how to do that yet)

# Possible fun stuff to do:
	* [ ] use GPU to calculate this stuff (marginal increase for double precision, maybe not worth) # [checkbox:unchecked]

# External dependencies:
	* [x] SFML (probably), maybe opengl or SDL # [checkbox:checked]
	* [ ] boost program options.
	* [ ] maybe other boost stuff.
