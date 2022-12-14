# LynFarm 
_Ver. 1.0 (July 2022)_

**LynFarm is developed for calculating asphalt pavement responses under the loadings of farming vehicles on tracks and/or flotation tires.**

**The code has been developed at the Technical University of Denmark in Kongens Lyngby. Hence the name. **

LynFarm is an add-on to the ALVA MATLAB-package developed by Asmus Skar, Sebastian Andersen and Julius Nielsen for pavement modeling.
[Link to ALVA](https://github.com/asmusskar/ALVA).

[Read more about ALVA](https://doi.org/10.21105/joss.02548).

## LynFarm roadmap

<div>
<img src="images/Roadmap_LynFarm_Ver.1.1.png" width="100%">
</div>
<p>
 <b>Figure 1:</b> LynFarm roadmap on a timeline. 
</p>

**The components are in short described below:**

- `user_input.m` - collects user inputs comprising of parameters about: the farming vehicle's geometry, the desired computational accuracy, and the linear elastic half-space resembling the pavement.

- `main.m` - is an executional file taking in the user's inputs and initiating a sequence of sub-procedures in the right order. In the end, the file gathers all information and provides the output. A plot of a specific response can optionally be chosen too.

- **ALVA** - This command initiates ALVA as a sub-routine.

- `filling_function.m` - runs a geometrical algorithm to fill a single rectangular patch with small identical micro-circles positioned in an intersecting yet non-overlapping lattice-like arrangement that optimally covers the patch's geometry.

- `pavement_response.m` & `grid assembly.m`- calculates the responses in a discrete half-space for a full farming vehicle by the means of interpolation and superposition algorithms. The discrete half-space is for presentation assembled as a finite 2D-grid containing many nodes. Each node contains the magnitudes of each response type. The `grid assembly.m` file collects each patch's contribution w.r.t. response magnitudes to each node element in the 2D-grid, and assembles it all onto a single common space.

- **PLOT** - This command is optional and can be turned on to construct and review a finite 2D-grid of a selected response for the farming vehicle once all calculations are done.

- **OUTPUT** - Last step in the code where results are presented in command window (and the plot is shown if selected to).

## Installation
- Download LynFarm and ALVA on your PC. [Link to ALVA](https://github.com/asmusskar/ALVA).
- Make the LynFarm a sub-folder in the ALVA folder.
- Open MATLAB.
- Navigate to the LynFarm folder and open `user_input.m` to define a vehicle via geometry inputs.
- Launch the `user_input.m` script to initiate LynFarm.
- _Optional: To run one of the example vehicles, open the **examples** folder and launch the script._
- LynFarm & ALVA are compatible with [OCTAVE](https://www.gnu.org/software/octave/index).

## Further work
This version of LynFarm is not final and further work is continuously being done to improve the performance.
Future versions are expected to be released.
