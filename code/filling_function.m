% -------------------------------------------------------------------------
% The filling function creates a square boundary with same dimensions as
% the user input for a patch in the configurations. Micro circles of same
% size with radius, a, are then assigned inside the boundary without any
% overlapping.
% -------------------------------------------------------------------------

function [Mcirc, Rec_area, Circ_area, Fill_ratio, Ncirc] = filling_function(FF)

% Defining 2D rotation matrix
Rx_rot = [cos(FF.theta) -sin(FF.theta); sin(FF.theta) cos(FF.theta)];

% -------------------------------------------------------------------------
% Creating a squared boundary which is filled with micro-circles
% -------------------------------------------------------------------------

for i = 1:length(FF.x_geo_sp)
    patch_no = i;
    
    % Patch construction
    Rec_length = FF.L_b;
    Rec_width = FF.W_b;
    xCenter = FF.x_geo_sp;
    yCenter = FF.y_geo_sp;
    xRight = xCenter + Rec_length/2;
    xLeft = xCenter - Rec_length/2;
    yTop = yCenter + Rec_width/2;
    yBottom = yCenter - Rec_width/2;
    
    % Coordinates for rectangular outline in matrix
    rec_coords = [xLeft xRight xRight xLeft; yBottom yBottom yTop yTop]';
    
    % Defining steps for multiplication of multiple circles
    dm = 2*FF.a;                   % Step size (No overlap at 2*r)
    step_x = length(Rec_length)*dm;  % Step length in x
    step_y = length(Rec_width)*dm;   % Step length in y
    
    % Creating centre coordinates for circles within rectangle
    x_center = xLeft+FF.a:step_x:xRight-FF.a;
    y_center = yBottom+FF.a:step_y:yTop-FF.a;
    
    % Calculating the fill of space ratio between circles and the rectangle
    Rec_area = Rec_length*Rec_width;
    Circ_area = size(y_center,2)*size(x_center,2)*FF.a^2*pi;
    Fill_ratio = Circ_area/Rec_area;
    Ncirc = length(x_center)*length(y_center);
    
    % Gathering all circle coodinates in matrix
    circ_coords = [];
    for j = 1:length(x_center)
        pair = [x_center(j)*ones(length(y_center),1) y_center'];
        circ_coords = [circ_coords; pair];
    end
    
    % 1) Translating circle coordinates by the center position of rectangle
    % 2) Rotating circle coordinates by rotation matrix
    % 3) Translating circle coordinates back to original position
    
    center = repmat([xCenter; yCenter], 1, length(circ_coords));
    rot_circ_coords = Rx_rot*(circ_coords'-center) + center;
    x_rot = rot_circ_coords(1,:);
    y_rot = rot_circ_coords(2,:);
    
    if FF.plot_filling_function == "True"
        
        %%% plotting the graphics
        figure('Position', [FF.plot_window]);
        hold all

        % Rectangle
        f = [1 2 3 4];
        Rec = patch('Faces',f,'Vertices',rec_coords,'EdgeColor','black', ...
            'FaceColor', 'none', 'LineWidth',1);
        rotate(Rec,[0 0 1],rad2deg(FF.theta),[xCenter yCenter 0]) % Rot. parameters

        % Circles
        r = FF.a;         % Radius
        ang=0:FF.ang_steps:2*pi;     % Steps to draw angular shape
        xu=r*cos(ang);      % Outline for x-coordinates
        yu=r*sin(ang);      % Outline for y-coordinates
        plot(x_rot'+xu,y_rot'+yu, 'b.', 'Markersize', 0.01);
        plot(x_rot', y_rot', '.r', 'MarkerSize', 10);
        axis equal;
        grid on;
        ax = gca;
        ax.XAxisLocation = 'origin';
        ax.YAxisLocation = 'origin';
        ax.XRuler.Axle.LineWidth = 1;
        ax.YRuler.Axle.LineWidth = 1;
        ax.GridColor = [0 0 0];
        ax.GridAlpha = 0.4;

        if FF.cnt == '1&2'
            title('Micro circles'' layout','FontSize',14)
        elseif FF.cnt == '1'
            title('Micro circles'' layout in configuration 1','FontSize',14)
        elseif FF.cnt == '2'
            title('Micro circles'' layout in configuration 2','FontSize',14)
        end

        % x-axis settings
        xticks([-xRight:50:xRight]);
        xlabel('(x) [mm]','Interpreter', 'Latex', 'FontSize', 12)
        xlim([-xRight*1.7 xRight*1.7])

        % y-axis seetings
        ylabel('(y) [mm]','Interpreter', 'Latex', 'FontSize', 12)
        ylim([-200 200])
    end

    % Gathering micro circle coordinates in global matrix
    if Ncirc == 1
        Mcirc{i} = [x_rot(1)' y_rot(1)'];
    elseif Ncirc > 1
        Mcirc{i} = [x_rot' y_rot'];
    end
    
end