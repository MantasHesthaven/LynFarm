# LynFarm _Ver. 1.0 (July 2022)_

**LynFarm is developed for calculating asphalt pavement responses under the loadings of farming vehicles on tracks and/or flotation tires.**

LynFarm is an add-on to the ALVA MATLAB-package developed by Asmus Skar, Sebastian Andersen and Julius Nielsen for pavement modeling.
[Link to ALVA](https://github.com/asmusskar/ALVA). [Read more about ALVA](https://doi.org/10.21105/joss.02548).

## The components of LynFarm are in short described below:

- 'user_input.m' - collects user inputs comprising of parameters about: the farming vehicle's geometry, the desired computational accuracy, and the linear elastic half-space resembling the pavement.

- 'main.m' - is an executional file taking in the user's inputs and initiating a sequence of sub-procedures in the right order. In the end, the file gathers all information and provides the output. A plot of a specific response can optionally be chosen too.

- 'filling_function.m' - runs a geometrical algorithm to fill a single rectangular patch with small identical micro-circles positioned in an intersecting yet non-overlapping lattice-like arrangement that optimally covers the patch's geometry.

- 'pavement_response.m' - calculates the responses in the defined half-space for a full tracked vehicle by the means of interpolation and superposition algorithms.

## Installation
- Download LynFarm and ALVA on your PC. [Link to ALVA](https://github.com/asmusskar/ALVA)
- Make the LynFarm a sub-folder in the ALVA folder.
- Open MATLAB
- Navigate to the LynFarm folder og open 'user_input.m' to define a vehicle.
- Open 'main.m' and launch the script to initiate LynFarm.
- ALVA is compatible with [OCTAVE](https://www.gnu.org/software/octave/index)

## Further work
This version of LynFarm is not final and further work is continuously being done to improve the performance.
Future versions are expected to be released.
