% -------------------------------------------------------------------------
% DESCRIPTION: User's inputs are provided here, which include:
%
% * Properties of the layered system
% * Coordinates of the calculation point
% * The farming vehicle's geometry
% * Patch stresses
% 
% -------------------------------------------------------------------------

close all; clear all; clc

% Creating path to ALVA sub-procedure folder
addpath(genpath(pwd))
addpath(genpath('../../basic'))
addpath(genpath('../Examples'))

% -------------------------------------------------------------------------
% Definition of pavement system
% -------------------------------------------------------------------------

% y-> x ------------------- Pavement surface ---------------------------- %
% |                              layer 1      ^      ^       ^
% v z                                         | z1   |       |
%                                             v      |       |
% ---------------------------------------------------|-------|-------------
%                                layer 2             | z2    |
%                                                    v       |
%------------------------------------------------------------|-------------
%                                                            .
%                                                            .
% -------------------------------------------------------------------------
%                                layer n-1                   .
% -----------------------------------------------------------|-------------
%                                layer n                     | infinity
%                                                            v

% -------------------------------------------------------------------------
% Select response analysis type
% -------------------------------------------------------------------------
analysis = 'Full';
% 1) 'Full'     : Conventional full integration with one-step Richardson
%                 extrapolation (for improved convergence near surface) 
%                 applied for evaluation of displacements and stresses
% 2) 'PolFit'   : Use polynomial fit technique to reduce integral length
%                 according to Andersen et al. (2018) for evaluation of
%                 surface displacements

% -------------------------------------------------------------------------
% Select interface 
% -------------------------------------------------------------------------
bond = 'Bonded';
% 1) 'Bonded'       : Full bonding between layers
% 2) 'Slip'         : Interface bonding factor
% 3) 'Frictionless' : No bonding between layers

% -------------------------------------------------------------------------
% Select additional settings
% -------------------------------------------------------------------------
plot_filling_function = 'False';  % Chose whether to visualize the filling
                                 % function of how micro circles are
                                 % allocated within a single patch

ang_steps = 0.1;      % Step size for graphically illustrating the circles.
                      % Choose between (0.1 - 1.0).
                      % A lower number increases resolution.

heat_map = 'True';    % Chose whether to plot a heatmap of the selected
                      % response. Heatmaps provide inituition on where
                      % maxima/minima are located.

heat_map_resp = 'Sigz';         % Select which response to visualize in the
                                % the heat map if it is set to 'True'.

% 'False': No plotting
% 'True' : Execute the plot

% 'Ux'   : Horizontal displacement in x-direction
% 'Uy'   : Horizontal displacement in y-direction
% 'Uz'   : Horizontal displacement in z-direction

% 'Sigx' : Horizontal stress in x-direction
% 'Sigy' : Horizontal stress in y-direction
% 'Sigz' : Vertical stress in z-direction
% 'Sigxy': Shear stress in the x,y plane
% 'Sigyz': Shear stress in the y,z plane
% 'Sigxz': Shear stress in the x,z plane

% 'Epsx' : Horizontal strain in x-direction
% 'Epsy' : Horizontal strain in y-direction
% 'Epsz' : Vertical strain in z-direction
% 'Epsxy': Shear strain in the x,y plane
% 'Epsyz': Shear strain in the y,z plane
% 'Epsxz': Shear strain in the x,z plane

% -------------------------------------------------------------------------
% Patch numbering system
% -------------------------------------------------------------------------

% ------------------------------ Example: ---------------------------------

%   [x]
%    ▲
%    |
%    |
%    #1  #2      #3  #4  #5  #6  #7  #8        
%    █   █       █   █   █   █   █   █     LEFT
%    #9  #10     #11 #12 #13 #14 #15 #16   SIDE
%    █   █       █   █   █   █   █   █         
%    |  
%    |-----------------------------------► [y] (Forward direction)
%    |
%    #17 #18     #19 #20 #21 #22 #23 #24        
%    █   █       █   █   █   █   █   █     RIGHT
%    #25 #26     #27 #28 #29 #30 #31 #32   SIDE
%    █   █       █   █   █   █   █   █      
%    |
%    |

% -------------------------------------------------------------------------
% Location of evaluation point: [x y z] in [mm]
% -------------------------------------------------------------------------
x = 4000; y = 1380; z = 0;

% The orientation is:
% x = 0 denotes the centerpoint of the leftmost patch.
% y = 0 denotes the vehicle's centerline.

% -------------------------------------------------------------------------
% Numerical parameters
% -------------------------------------------------------------------------
N   = 300;  % No. of Bessel zero points in numerical integration
n   = 30;   % No. of Gauss points between zero points
a   = 10;   % Radius of homogeneous micro circles [mm]
q_m = 0.7;  % Uniform load pressure in [MPa] for homogeneous micro-circles

% -------------------------------------------------------------------------
% Select pavement material properties and layers (minimum two layers required)
% -------------------------------------------------------------------------
zi = [165 (165+330)];     % Individual layer heights, cf. ALVA
E  = [3000 200 40];       % Young's modulus for pavement layers, cf. ALVA
nu = [0.35, 0.35, 0.35];  % Layer Poisson's ratio [-]
kh = [1e9, 1e9];          % Interface bonding/horizontal spring [MPa/mm]

% -------------------------------------------------------------------------
% Define external configuration geometry
% -------------------------------------------------------------------------
dist_conf1 = 3300;  % Distance between the outer edges of the patches on
                    % each side of vehicle for configuration 1 [mm]

dist_conf2 = 3300;  % Distance between the outer edges of the patches on
                    % each side of vehicle for configuration 2 [mm]

dist_conf12 = 3080; % Center-Center distance between last patch in
                    % configuration 1 and first patch in configuration 2 [mm]

theta = 1.0485;     % Angular rotation to the X-axis line valid for all
                    % patches [Rad]

% -------------------------------------------------------------------------
% Define internal configuration geometry
% -------------------------------------------------------------------------

% For configuration 1
no_patches_conf1 = 4;    % No. of patches in configuration (must be even)
pitch_conf1 = 160;       % Pitch between patches [mm]
edge_offset_conf1 = 210; % Offset from outer edge [mm]
cent_offset_conf1 = 130; % Offset from centerline [mm]
L_conf1 = 300;           % Length of patch in this configuration [mm]
W_conf1 = 60;            % Width of patch in this configuration [mm]

% For configuration 2
no_patches_conf2 = 26;   % No. of patches in configuration (must be even)
pitch_conf2 = 160;       % Pitch between patches [mm]
edge_offset_conf2 = 210; % Offset from outer edge [mm]
cent_offset_conf2 = 150; % Offset from centerline [mm]
L_conf2 = 300;           % Length of patch in this configuration [mm]
W_conf2 = 60;            % Width of patch in this configuration [mm]

% -------------------------------------------------------------------------
% Define load distribution 
% -------------------------------------------------------------------------   
% Load vector for all patches in [MPa].
% Length must comply with total number of patches.
% Specific patches can be neglected by inserting 0.

q_tot = [0.559 0.559 0.397 0.397 0.055 0.055 0.231 0.231 0.212 0.231 0.231...
         0.055 0.055 0.397 0.397 0.559 0.559 0.397 0.397 0.055 0.055 0.231...
     	 0.231 0.212 0.231 0.231 0.055 0.055 0.397 0.397 0.559 0.559 0.397...
     	 0.397 0.055 0.055 0.231 0.231 0.212 0.231 0.231 0.055 0.055 0.397...
     	 0.397 0.559 0.559 0.397 0.397 0.055 0.055 0.231 0.231 0.212 0.231...
     	 0.231 0.055 0.055 0.397 0.397];

% -------------------------------------------------------------------------
% Define grid discretization
% (The lower the discretization, the higher performance is achieved)
% -------------------------------------------------------------------------
inc = 10; % Distance between mesh points in grid.
          % Minimum allowed increment size must be equal to the rounding
          % precision based on the input geometry.  

% -------------------------------------------------------------------------
% Activating main.m file to initiate LynFarm
% -------------------------------------------------------------------------
run('main.m')
