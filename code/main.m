% LynFarm Ver. 1.0 (July 2022)
% -------------------------------------------------------------------------

% DESCRIPTION: Main executional file taking in the user's inputs and
% initiating a sequence of sub-procedures in the following order:

% (1) Loading user inputs
% (2) Executing filling function
% (3) Calculating pavement responses for a track
% (4) Providing results in form of text (and via an optional plot)
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Performing vital checks of user inputs for the program to run correctly
% -------------------------------------------------------------------------

% In the case that only one configuration is selected, it is by default set
% to Configuration 1. Rest is set to 0.
if no_patches_conf1 == 0 && no_patches_conf2 > 0
    no_patches_conf1 = no_patches_conf2;
    pitch_conf1 = pitch_conf2;
    edge_offset_conf1 = edge_offset_conf2;
    cent_offset_conf1 = cent_offset_conf2;
    L_conf1 = L_conf2;
    W_conf1 = W_conf2;

    no_patches_conf2 = 0; pitch_conf2 = 0; edge_offset_conf2 = 0;
    cent_offset_conf2 = 0; L_conf2 = 0; W_conf2 = 0;
end

% Making sure that the finite grid discretization matches the rounding
% precision of the user's input geometry values
list_ext = [dist_conf1, dist_conf2, dist_conf12];
list_conf1 = [pitch_conf1 edge_offset_conf1 cent_offset_conf1 L_conf1 W_conf1];
list_conf2 = [pitch_conf2 edge_offset_conf2 cent_offset_conf2 L_conf2 W_conf2];

rnd_prec = mod([list_ext, list_conf1, list_conf2]/inc,1);

for i = 1:length(rnd_prec)
    if rnd_prec(i) > 0
        fprintf(['\nGrid discretization does not comply with the\n' ...
                 'rounding precision of input geometry.\n'])
        fprintf(['\nTry minimizing the increment size between mesh points\n' ...
                 'OR using a rougher rounding precision.\n'])
    return;
    end
end

% Rounding pitch distances in case their halves are of higher rounding
% precision than the increment size. This operation is viable for the
% discrete translations where every second row of patches are offset by
% a pitch/2 in x.
if mod((pitch_conf1/2)/inc,1) > 0
    pitch_conf1 = round(pitch_conf1/2/inc,0)*2*inc;
elseif mod((pitch_conf2/2)/inc,1) > 0
    pitch_conf2 = round(pitch_conf2/2/inc,0)*2*inc;
end

% Checking if patches are defined correctly
if no_patches_conf1+no_patches_conf2 == 0
    fprintf(['\nThere are currently 0 patches being modelled.\n'...
        'Please revise.\n'])
    return;
elseif no_patches_conf1 < 0
    fprintf(['\nNegative amount of patches have been defined for\n'...
        'configuration 1. Please revise.\n'])
    return;
elseif no_patches_conf2 < 0
    fprintf(['\nNegative amount of patches have been defined for\n'...
        'configuration 2. Please revise.\n'])
    return;
end

% Checking if user input for number of patches is even
if mod(no_patches_conf1,2) > 0
    fprintf(['\nNumber of patches for configuration 1 is odd.\n' ...
             'Please correct to an even amount.\n'])
    return;
elseif mod(no_patches_conf2,2) > 0
    fprintf(['\nNumber of patches for configuration 2 is odd.\n' ...
             'Please correct to an even amount.\n'])
    return;
end

% Checking whether load vector contains values for all patches
if length(q_tot) < (no_patches_conf1+no_patches_conf2)*2
    fprintf(['\nLoad vector, q_tot, is missing values for some patches.\n'...
             'Please revise.\n'])
    return;
elseif length(q_tot) > (no_patches_conf1+no_patches_conf2)*2
    fprintf(['\nLoad vector, q_tot, contains too many values than for the\n'...
             'specified amount of patches.\n'])
    return;
end

% Checking if patches are wide and long enough to contain atleast a single
% micro circle
d = 2*a;
if no_patches_conf1 > 0
    if L_conf1 < d || W_conf1 < d
        fprintf(['\nPatch dimensions in configuration 1 are too small\n' ...
            'to fit atleast one micro circle. Choose a smaller circle radius.\n'])
        return;
    end
end

if no_patches_conf2 > 0
    if L_conf2 < d || W_conf2 < d
        fprintf(['\nPatch dimensions in configuration 2 are too small\n' ...
            'to fit atleast one micro circle. Choose a smaller circle radius.\n'])
        return;
    end
end

% -------------------------------------------------------------------------
% Geometry calculations for the configurations
% -------------------------------------------------------------------------
no_patches = length(q_tot); % Total no. of patches for the vehicle

if no_patches_conf1 > 0
    % Defining centre coordinates for patches in configuration 1
    y_geo_conf1 = repelem(round((dist_conf1/2-L_conf1/2)/inc)*inc, no_patches_conf1/2);
    x_geo_conf1 = [0, linspace(pitch_conf1,(no_patches_conf1/2-1)*pitch_conf1,...
                   no_patches_conf1/2-1)];
    x_geo_conf1_end = x_geo_conf1(end);
else
    y_geo_conf1 = [];
    x_geo_conf1 = [];
    dist_conf12 = 0;
    x_geo_conf1_end = 0;
end

if no_patches_conf2 > 0
    % Defining centre coordinates for patches in configuration 2
    y_geo_conf2 = repelem(round((dist_conf2/2-L_conf2/2)/inc)*inc, no_patches_conf2/2);
    x_geo_conf2 = [dist_conf12+x_geo_conf1_end,...
               linspace(pitch_conf2+dist_conf12+x_geo_conf1_end,...
              (no_patches_conf2/2-1)*pitch_conf2+dist_conf12+x_geo_conf1_end,...
               no_patches_conf2/2-1)]-round(pitch_conf2/2/inc)*inc;
else
    y_geo_conf2 = [];
    x_geo_conf2 = [];
end

% Combining geometry data for one row of patches from both configurations
y_geo = [y_geo_conf1, y_geo_conf2];
x_geo = [x_geo_conf1, x_geo_conf2];

if no_patches_conf1 == 0 && no_patches_conf2 > 0
    no_patches_conf1 = no_patches_conf2;
    pitch_conf1 = pitch_conf2;
    edge_offset_conf1 = edge_offset_conf2;
    cent_offset_conf1 = cent_offset_conf2;
    L_conf1 = L_conf2;
    W_conf1 = W_conf2;

    no_patches_conf2 = 0;
    pitch_conf2 = 0;
    edge_offset_conf2 = 0;
    cent_offset_conf2 = 0;
    L_conf2 = 0;
    W_conf2 = 0;
end 

% -------------------------------------------------------------------------
% Defining local grid size based on vehicle geometry and patch locations
% -------------------------------------------------------------------------
% Rounding up to ensure that the local grid encapsules the whole vehicle
grid_max_round = ceil(max([x_geo_conf1,x_geo_conf2])/100)*100;
grid_size = grid_max_round*2; % [mm] 

% Checking if the coordinates of evaluation points are within the grid
if x > max(x_geo)+grid_size/2 || x < min(x_geo)-grid_size/2
     fprintf(['\nx-coordinate for evaluation point is beyond the calculation grid\n' ...
            'Please select another x-coordinate.\n'])
     return;

elseif y > max(y_geo)+grid_size/2 || y < min(y_geo)-grid_size/2
    fprintf(['\ny-coordinate for evaluation point is beyond the calculation grid\n' ...
            'Please select another y-coordinate.\n'])
     return;
end

% -------------------------------------------------------------------------
% Calculation points for Alva to be used for the interpolation scheme
% -------------------------------------------------------------------------

% Generating a 100 point vector in the x-direction with offset distances
% w.r.t. to the radius, a, from the micro-circle's center.
% The first third of points are placed within the micro-circle's boundary,
% while the rest are beyond.
% Points are specifically assigned to be more densely spaced apart close to
% the micro-circle's edge - and further away from each other at greater
% offset distances.
% The first offset distance should be 0, although, the calculation engine
% cannot handle a true 0. Instead, a numerical 0 is defined (i.e. a very
% small number).

xp =...
[10^(-10);0.0001;0.01;0.10;0.15;0.20;0.25;0.30;0.35;0.40;0.45;0.50;0.55;...
0.60;0.65;0.70;0.75;0.80;0.82;0.84;0.86;0.88;0.90;0.91;0.92;0.93;0.94;...
0.95;0.96;0.97;0.98;0.99;1.00;1.01;1.02;1.03;1.04;1.05;1.06;1.07;1.08;...
1.09;1.10;1.12;1.14;1.16;1.18;1.2;1.3;1.4;1.5;1.6;1.7;1.8;1.9;2.0;2.2;...
2.4;2.6;2.8;3.0;3.2;3.4;3.6;3.8;4.0;4.2;4.4;4.6;4.8;5;6;7;8;9;10;12;14;...
16;18;20;22;24;26;28;30;35;40;45;50;60;70;80;90;100;200;500;1000;2000;10000]*a;

% Points in the y-direction are all set to 0, while points in the
% z-direction are held at a constant depth.
yp = zeros(length(xp),1);    
zp = repmat(z,length(xp),1);

% -------------------------------------------------------------------------
% Creating ALVA struct file and initiating ALVA sub-routine
% -------------------------------------------------------------------------

% The ALVA sub-routine is initiated for only one micro circle,
% which will be the sample representative for all micro circles.
alva.analysis = analysis;
alva.bond = bond;
alva.N = N;
alva.n = n;
alva.nu = nu;
alva.kh = kh;
alva.zi = zi;
alva.E = E;
alva.a = a;
alva.q = q_m;
alva.xp = xp;
alva.Xl = [0.0 0.0]; % Positioning the sample micro circle in orego
alva.Xd = [xp,yp,zp]; % Constructing coordinate matrix for ALVA

% Launching sub-routine
alva = init_LET(alva);

% -------------------------------------------------------------------------
% Call filling function
% -------------------------------------------------------------------------
% Placing the sample patch in orego
y_geo_sp = 0; x_geo_sp = 0;    

% Creating struct variable for this function
FF.a = a;
FF.x = x; 
FF.y = y;
FF.y_geo_sp = y_geo_sp;
FF.x_geo_sp = x_geo_sp;
FF.theta = theta;
FF.plot_filling_function = plot_filling_function;
FF.ang_steps = ang_steps;

% Finding screen resolution for plotting
set(0,'units','pixels')
pix_ss = get(0,'screensize');

wind_w = 700; % Width of plot window
wind_h = 400; % Height of plot window

% Determining if patch sizes are equal or differ for configurations 1 & 2
if no_patches_conf1 > 0
    L_conf12 = L_conf1;
elseif no_patches_conf2 > 0
    L_conf12 = L_conf2;
end

% Number of patches in both configurations
no_patches_conf12 = no_patches_conf1+no_patches_conf2;

if L_conf1 == L_conf2 && W_conf1 == W_conf2

    % Executing filling function for configurations 1 & 2
    FF.L_b = L_conf1;
    FF.W_b = W_conf1;
    FF.cnt = '1&2';
    FF.plot_window = [(pix_ss(3)-wind_w)/2 (pix_ss(4)-wind_h)/2,...
        wind_w wind_h]; % Size of plot window
    [Mcirc_conf12, Rec_area_conf12, Circ_area_conf12, Fill_ratio, Ncirc] = filling_function(FF);

    % Calculating area ratio (patch/circles)
    A_diff = (Rec_area_conf12-Circ_area_conf12)/Circ_area_conf12;

    % Assigning circle positions and radii to arrays
    Xl_load_conf12 = zeros(size(Mcirc_conf12{1}));
    Xl_load_conf12(:,:) = Mcirc_conf12{1};
    Xl_load_conf12 = reshape(permute(Xl_load_conf12, [1,3,2]), [], 2);
    a_load_conf12  = repmat(a,size(Xl_load_conf12,1),1);

    alva.Mcirc_conf12 = Mcirc_conf12;
    alva.Xl_load_conf12 = Xl_load_conf12;
    alva.a_load_conf12 = a_load_conf12;

elseif  L_conf1 ~= L_conf2 || W_conf1 ~= W_conf2

    if no_patches_conf1 > 0
        % Executing filling function for configuration 1
        FF.L_b = L_conf1;
        FF.W_b = W_conf1;
        FF.cnt = '1';
        FF.plot_window = [pix_ss(3)/2-wind_w (pix_ss(4)-wind_h)/2,...
            wind_w wind_h]; % Size of plot window
        [Mcirc_conf1, Rec_area_conf1, Circ_area_conf1, Fill_ratio, Ncirc] = filling_function(FF);

        % Display statement
        if plot_filling_function == "True"
            fprintf('\n***\nSample patch for configuration 1\n');
            fprintf( ['\nSpace fill ratio of circles/rectangle: %.3f\n', ...
                'No. of micro circles: %.0f\n'], ...
                Fill_ratio, Ncirc);
        else

        end

        % Calculating area ratio (patch/circles) for configuration 1
        A_diff_conf1 = (Rec_area_conf1-Circ_area_conf1)/Circ_area_conf1;

        % Assigning circle positions and radii to arrays
        Xl_load_conf1 = zeros(size(Mcirc_conf1{1}));
        Xl_load_conf1(:,:) = Mcirc_conf1{1};
        Xl_load_conf1 = reshape(permute(Xl_load_conf1, [1,3,2]), [], 2);
        a_load_conf1  = repmat(a,size(Xl_load_conf1,1),1);

        alva.Mcirc_conf1 = Mcirc_conf1;
        alva.Xl_load_conf1 = Xl_load_conf1;
        alva.a_load_conf1 = a_load_conf1;

    elseif no_patches_conf2 > 0
        % Executing filling function for configuration 2
        FF.L_b = L_conf2;
        FF.W_b = W_conf2;
        FF.cnt = '2';
        FF.plot_window = [pix_ss(3)/2 (pix_ss(4)-wind_h)/2,...
            wind_w wind_h]; % Size of plot window
        [Mcirc_conf2, Rec_area_conf2, Circ_area_conf2, Fill_ratio, Ncirc] = filling_function(FF);

        % Display statement
        if plot_filling_function == "True"
            fprintf('\n***\nSample patch for configurations 2\n');
            fprintf( ['\nSpace fill ratio of circles/rectangle: %.3f\n', ...
                'No. of micro circles: %.0f\n'], ...
                Fill_ratio, Ncirc);
        else

        end

        % Calculating area ratio (patch/circles) for configuration 2
        A_diff_conf2 = (Rec_area_conf2-Circ_area_conf2)/Circ_area_conf2;

        % Assigning circle positions and radii to arrays
        Xl_load_conf2 = zeros(size(Mcirc_conf2{1}));
        Xl_load_conf2(:,:) = Mcirc_conf2{1};
        Xl_load_conf2 = reshape(permute(Xl_load_conf2, [1,3,2]), [], 2);
        a_load_conf2  = repmat(a,size(Xl_load_conf2,1),1);

        alva.Mcirc_conf2 = Mcirc_conf2;
        alva.Xl_load_conf2 = Xl_load_conf2;
        alva.a_load_conf2 = a_load_conf2;
    end
    % Indexing correct area ratios among configurations 1 and 2
    if no_patches_conf1 == 0
        A_diff_idx = [repelem(A_diff_conf2,no_patches_conf2/2)];

    elseif no_patches_conf2 == 0
        A_diff_idx = [repelem(A_diff_conf1,no_patches_conf1/2)];

    elseif no_patches_conf1 > 0 && no_patches_conf2 > 0
        A_diff_idx = [repelem(A_diff_conf1,no_patches_conf1/2),...
                      repelem(A_diff_conf2,no_patches_conf1/2)];
    end

    A_diff = repmat(A_diff_idx,1,4);
end

% -------------------------------------------------------------------------
% Preparatory calculations
% -------------------------------------------------------------------------

% Load vector where the missing empty area among the micro circles inside
% the patches are compensated.
q_all = q_tot.*(1+A_diff);

% Size of square grid
Length = grid_size;
Width = Length;

% No. of steps in grid
xm = Length/inc; % In x [mm]
ym = Width/inc;  % In y [mm]

% Increment size
dx = Length/xm; % In x [mm]
dy = Width/ym;  % In y [mm]

% Discretization in 2D
xd = -Length/2:inc:Length/2; % In x [mm]
yd = -Width/2:inc:Width/2;   % In y [mm]

% -------------------------------------------------------------------------
% Calling sub-script for carrying out pavement analysis 
% -------------------------------------------------------------------------
% Assigning last variables to ALVA struct
alva.dx = dx; alva.dy = dy; alva.xd = xd; alva.yd = yd; alva.q_all = q_all;
alva.z = z; alva.theta = theta; alva.x_geo = x_geo; alva.y_geo = y_geo;
alva.L_conf1 = L_conf1; alva.L_conf2 = L_conf2; alva.W_conf1 = W_conf1;
alva.W_conf2 = W_conf2; alva.no_patches = no_patches; alva.dist_conf1 = dist_conf1;
alva.dist_conf2 = dist_conf2; alva.dist_conf12 = dist_conf12;
alva.pitch_conf1 = pitch_conf1; alva.pitch_conf2 = pitch_conf2;
alva.edge_offset_conf1 = edge_offset_conf1; alva.edge_offset_conf2 = edge_offset_conf2;
alva.cent_offset_conf1 = cent_offset_conf1; alva.cent_offset_conf2 = cent_offset_conf2;
alva.no_patches_conf1 = no_patches_conf1; alva.no_patches_conf2 = no_patches_conf2;
alva.no_patches_conf12 = no_patches_conf12; alva.inc = inc;

% Storing results in its own structure array
[sigx, sigy, sigz, sigxy, sigyz, sigxz, epsx, epsy, epsz, epsxy, epsyz,...
 epsxz, ux, uy, uz, X_glob, Y_glob, x_start, y_start, x_stop, y_stop] = pavement_response(alva);

res{1}  = ux;
res{2}  = uy;
res{3}  = uz;
res{4}  = sigx;
res{5}  = sigy;
res{6}  = sigz;
res{7}  = sigxy;
res{8}  = sigyz;
res{9}  = sigxz;
res{10} = epsx;
res{11} = epsy;
res{12} = epsz;
res{13} = epsxy;
res{14} = epsyz;
res{15} = epsxz;

% -------------------------------------------------------------------------
% (1) Response results for user chosen coordinates
% -------------------------------------------------------------------------

% Locating the user chosen coordinates as position in the global grid.
% In case the coordinates are not a node point in the grid, the closest
% value is selected.
[x_idx,x_pos] = min(abs(X_glob(1,:)-x));
[y_idx,y_pos] = min(abs(Y_glob(:,1)-y));

% Print statement output
fprintf('\n***\nResponse results for:\n (x = %.0f; y = %.0f; z = %.0f)\n',...
         x, y, z);

fprintf(['\nHorizontal displacement in x: %.2f mm\n',...
         'Horizontal displacement in y: %.2f mm\n',...
         'Vertical displacement in z: %.2f mm\n',...
         'Horizontal stress in x-direction: %.2f MPa\n',...
         'Horizontal stress in y-direction: %.2f MPa\n',...
         'Vertical stress in z-direction: %.2f MPa\n',...
         'Shear stress in the x,y plane: %.2f MPa\n',...
         'Shear stress in the y,z plane: %.2f MPa\n',...
         'Shear stress in the x,z plane: %.2f MPa\n',...
         'Horizontal strain in x-direction: %.0f microstrain\n',...
         'Horizontal strain in y-direction: %.0f microstrain\n',...
         'Vertical strain in z-direction: %.0f microstrain\n',...
         'Shear strain in the x,y plane: %.0f microstrain\n',...
         'Shear strain in the y,z plane: %.0f microstrain\n',...
         'Shear strain in the x,z plane: %.0f microstrain\n***\n'],...
         res{1}(y_pos,x_pos), res{2}(y_pos,x_pos), res{3}(y_pos,x_pos),...
         res{4}(y_pos,x_pos), res{5}(y_pos,x_pos), res{6}(y_pos,x_pos),...
         res{7}(y_pos,x_pos), res{8}(y_pos,x_pos), res{9}(y_pos,x_pos),...
         res{10}(y_pos,x_pos), res{11}(y_pos,x_pos), res{12}(y_pos,x_pos),...
         res{13}(y_pos,x_pos), res{14}(y_pos,x_pos), res{15}(y_pos,x_pos));

% -------------------------------------------------------------------------
% (2) Heatmap plot of user chosen response
% -------------------------------------------------------------------------

if heat_map == 'True'
    % Extracting the correct response matrix for plot
    list_resp = {'Sigx','Sigy','Sigz','Sigxy','Sigyz','Sigxz',...
        'Epsx','Epsy','Epsz','Epsxy','Epsyz','Epsxz',...
        'Ux','Uy','Uz'};
    idx_resp = find(strcmp(list_resp, heat_map_resp));

    %%% inkluder angivelse af punkt i koordinaterne som User har valgt
    max_resp = max(max(res{idx_resp}(:,:)))+0.1;
    figure('Position', [200 200 1200 600])
    heatmap = surface(X_glob,Y_glob,res{idx_resp}(:,:), 'EdgeAlpha',0.01);
    heatmap.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hold on
    plot3(x, y, max_resp,'ko','MarkerFaceColor',[1 0 0],'linewidth',3,'MarkerSize',10)
    set(gca,'FontSize',11)
    title('Heat map of the defined farming vehicle','Interpreter', 'Latex','FontSize',15)
    subtitle(['Selected response: ' heat_map_resp],'Interpreter', 'Latex','FontSize',14)
    xlabel('(x) [mm]','Interpreter', 'Latex', 'FontSize', 14)
    ylabel('(y) [mm]','Interpreter', 'Latex', 'FontSize', 14)
    xlim([X_glob(1,max(x_start))-L_conf12 X_glob(1,min(x_stop))+L_conf12])
    ylim([Y_glob(max(y_start),1) Y_glob(min(y_stop),1)])
    legend({'Location of evaluation point'},'Interpreter', 'Latex','FontSize',13)
    ax = gca;
    hcb = colorbar;
    if ismember(idx_resp, [4,5,6,7,8,9])
        title(hcb,'[MPa]','Interpreter', 'Latex','FontSize',13);
    elseif ismember(idx_resp, [10,11,12,13,14,15])
        title(hcb,'[$\mu m/m$]','Interpreter','Latex','FontSize',13);
    elseif ismember(idx_resp, [1,2,3])
        title(hcb,'[mm]','Interpreter', 'Latex','FontSize',13);
    end
    hold off
end
