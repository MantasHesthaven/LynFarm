% -------------------------------------------------------------------------
% This file takes results for one micro circle from the ALVA sub-routine
% and via an interpolation scheme applies it to all patches in the grid.
% By utilizing the principle of superposition, the final responses are in
% the end summed to construct the full picture of the response space for
% the chosen vehicle.
% -------------------------------------------------------------------------

function [sigx_glob, sigy_glob, sigz_glob, sigxy_glob, sigyz_glob, sigxz_glob,...
          epsx_glob, epsy_glob, epsz_glob, epsxy_glob, epsyz_glob, epsxz_glob,...
          ux_glob, uy_glob, uz_glob, X_glob, Y_glob, x_start, y_start, x_stop, y_stop] = pavement_response(alva)

dx = alva.dx; dy = alva.dy; xd = alva.xd; yd = alva.yd; q_all = alva.q_all; z = alva.z;
theta = alva.theta; x_geo = alva.x_geo; y_geo = alva.y_geo;
L_conf1 = alva.L_conf1; L_conf2 = alva.L_conf2; W_conf1 = alva.W_conf1; W_conf2 = alva.W_conf2;
no_patches = alva.no_patches; dist_conf1 = alva.dist_conf1;
dist_conf2 = alva.dist_conf2; dist_conf12 = alva.dist_conf12;
pitch_conf1 = alva.pitch_conf1; pitch_conf2 = alva.pitch_conf2;
edge_offset_conf1 = alva.edge_offset_conf1; edge_offset_conf2 = alva.edge_offset_conf2;
cent_offset_conf1 = alva.cent_offset_conf1; cent_offset_conf2 = alva.cent_offset_conf2;
no_patches_conf1 = alva.no_patches_conf1; no_patches_conf2 = alva.no_patches_conf2;
no_patches_conf12 = alva.no_patches_conf12; inc = alva.inc;

% -------------------------------------------------------------------------
% Gathering needed responses from ALVA sub-routine
% -------------------------------------------------------------------------
M(:,1) = alva.sigx;
M(:,2) = alva.sigy;
M(:,3) = alva.sigz;
M(:,4) = alva.sigxz;
M(:,5) = alva.ux;
M(:,6) = alva.uy;
M(:,7) = alva.uz;

% Gathering all response results in one array
V = [M(1:end,1:end); zeros(1,size(M,2))];

% As only the x-axis is considered, the points are already in radial
% distances. The response must be 0 for a great distance. The last position
% in the array, which is 0, is therefore located at e.g. x*10^3.
x0 = [alva.xp; 10^3*max(alva.xp)];

% Creating grid of coordinates
[X,Y] = meshgrid(alva.xd,alva.yd);

% Creating global coordinates with individual load positions - either for
% configurations 1 and 2, or a common one if configurations are alike

if L_conf1 == L_conf2 && W_conf1 == W_conf2
    Mcirc_conf12 = alva.Mcirc_conf12;
    Xl_load_conf12 = alva.Xl_load_conf12;
    a_load_conf12 = alva.a_load_conf12;

    n_confs = 1;
    R_conf12 = cell(1,size(Xl_load_conf12,1));
    X_patch_conf12 = cell(1,size(Xl_load_conf12,1));
    Y_patch_conf12 = cell(1,size(Xl_load_conf12,1));

    for i = 1:size(Xl_load_conf12,1)
        X_patch_conf12{i} = X - Xl_load_conf12(i,1);
        Y_patch_conf12{i} = Y - Xl_load_conf12(i,2);

        % Calculating the radial distances to each point in grid for interpolation
        R_conf12{i} = sqrt(X_patch_conf12{i}.^2+Y_patch_conf12{i}.^2);

        % If there is a radial distance in the grid where: r = 0, then it is
        % replaced by the first point in x0 because the interpolation will not
        % work otherwise with a distance of 0.
        if nnz(~R_conf12{i}) > 0
            [row,col] = find(R_conf12{i}==0);
            R_conf12{i}(row,col) = x0(1);
        end
    end

elseif  L_conf1 ~= L_conf2 || W_conf1 ~= W_conf2
    if no_patches_conf1 > 0 && no_patches_conf2 > 0
        n_confs = 2;
    else
        n_confs = 1;
    end

    if no_patches_conf1 > 0
        Mcirc_conf1 = alva.Mcirc_conf1;
        Xl_load_conf1 = alva.Xl_load_conf1;
        a_load_conf1 = alva.a_load_conf1;

        R_conf1 = cell(1,size(Xl_load_conf1,1));
        X_patch_conf1 = cell(1,size(Xl_load_conf1,1));
        Y_patch_conf1 = cell(1,size(Xl_load_conf1,1));

        for i = 1:size(Xl_load_conf1,1)
            X_patch_conf1{i} = X - Xl_load_conf1(i,1);
            Y_patch_conf1{i} = Y - Xl_load_conf1(i,2);

            % Calculating the radial distances to each point in grid for interpolation
            R_conf1{i} = sqrt(X_patch_conf1{i}.^2+Y_patch_conf1{i}.^2);

            if nnz(~R_conf1{i}) > 0
                [row_conf1,col_conf1] = find(R_conf1{i}==0);
                R_conf1{i}(row_conf1,col_conf1) = x0(1);
            end
        end

    elseif no_patches_conf2 > 0
        Mcirc_conf2 = alva.Mcirc_conf2;
        Xl_load_conf2 = alva.Xl_load_conf2;
        a_load_conf2 = alva.a_load_conf2;

        R_conf2 = cell(1,size(Xl_load_conf2,1));
        X_patch_conf2 = cell(1,size(Xl_load_conf2,1));
        Y_patch_conf2 = cell(1,size(Xl_load_conf2,1));

        for i = 1:size(Xl_load_conf2,1)
            X_patch_conf2{i} = X - Xl_load_conf2(i,1);
            Y_patch_conf2{i} = Y - Xl_load_conf2(i,2);

            % Calculating the radial distances to each point in grid for interpolation
            R_conf2{i} = sqrt(X_patch_conf2{i}.^2+Y_patch_conf2{i}.^2);

            if nnz(~R_conf2{i}) > 0
                [row_conf2,col_conf2] = find(R_conf2{i}==0);
                R_conf2{i}(row_conf2,col_conf2) = x0(1);
            end
        end
    end
end

% -------------------------------------------------------------------------
% Interpolation algorithm in 2D
% -------------------------------------------------------------------------

% Defining whether to interpolate once if configurations 1 and 2 are alike,
% or twice if they are different.

clear n
zd = z; % Possibility for future modelling in 3D by including a vector,
        % zd, in the z-direction. For now, it is a fixed value.

for n = 1:n_confs
    clear Resx Resy Resz Resxz Resux Resuy Resuz Xl_load

    if L_conf1 == L_conf2 && W_conf1 == W_conf2
        Xl_load = Xl_load_conf12;
        R = R_conf12;
        a_load = a_load_conf12;
        fprintf('\nInterpolation and Transformation Time for Configurations 1 & 2\n')
    elseif L_conf1 ~= L_conf2 || W_conf1 ~= W_conf2
        if n == 1
            if no_patches_conf1 > 0
                Xl_load = Xl_load_conf1;
                R = R_conf1;
                a_load = a_load_conf1;
                fprintf('\nInterpolation and Transformation Time for Configuration 1\n')
            elseif no_patches_conf1 == 0
                Xl_load = Xl_load_conf2;
                R = R_conf2;
                a_load = a_load_conf2;
                fprintf('\nInterpolation and Transformation Time for Configuration 2\n')
            end
        elseif n == 2
            Xl_load = Xl_load_conf2;
            R = R_conf2;
            a_load = a_load_conf2;
            fprintf('\nInterpolation and Transformation Time for Configuration 2\n')
        end
    end

    tic
    % Pre-allocating array sizes used in interpolation
    int_pre = zeros(length(xd),length(yd),size(Xl_load,1),length(zd));
    Resx = int_pre; Resy = int_pre; Resz = int_pre; Resxz = int_pre;
    Resux = int_pre; Resuy = int_pre; Resuz = int_pre;

    for h = 1:length(zd) % Expansion in the third dimension is currently unused.
        z = zd(h);       % But can be implemented through for-looping in z.

        for k = 1:size(Xl_load,1)
            Resx(:,:,k,h) = interp1(log10(x0),V(:,1),log10(R{k}(1:length(xd), 1:length(yd))),'linear');
            Resy(:,:,k,h) = interp1(log10(x0),V(:,2),log10(R{k}(1:length(xd), 1:length(yd))),'linear');
            Resz(:,:,k,h) = interp1(log10(x0),V(:,3),log10(R{k}(1:length(xd), 1:length(yd))),'linear');
            Resxz(:,:,k,h) = interp1(log10(x0),V(:,4),log10(R{k}(1:length(xd), 1:length(yd))),'linear');
            Resux(:,:,k,h) = interp1(log10(x0),V(:,5),log10(R{k}(1:length(xd), 1:length(yd))),'linear');
            Resuy(:,:,k,h) = interp1(log10(x0),V(:,6),log10(R{k}(1:length(xd), 1:length(yd))),'linear');
            Resuz(:,:,k,h) = interp1(log10(x0),V(:,7),log10(R{k}(1:length(xd), 1:length(yd))),'linear');
        end
    end

    % ---------------------------------------------------------------------
    % Tranformation to Cartersian coordinates
    % ---------------------------------------------------------------------
    % Making vectors vertical
    if size(xd,1) < size(xd,2)
        xd = xd';
    end

    if size(yd,1) < size(yd,2)
        yd = yd';
    end

    if size(zd,1) < size(zd,2)
        zd = zd';
    end

    if size(alva.zi,1) > size(alva.zi,2)
        alva.zi = alva.zi';
    end

    %%% Rotation in 3-D
    % Pre-allocating empty arrays
    uxs = zeros(size(Resx,1),size(Resx,2),size(Xl_load,1),length(zd));
    uys = zeros(size(Resy,1),size(Resy,2),size(Xl_load,1),length(zd));
    uzs = zeros(size(Resz,1),size(Resz,2),size(Xl_load,1),length(zd));

    sigxs = zeros(size(Resx,1),size(Resx,2),size(Xl_load,1),length(zd));
    sigys = zeros(size(Resy,1),size(Resy,2),size(Xl_load,1),length(zd));
    sigzs = zeros(size(Resz,1),size(Resz,2),size(Xl_load,1),length(zd));
    sigxys = zeros(size(Resx,1),size(Resx,2),size(Xl_load,1),length(zd));
    sigyzs = zeros(size(Resy,1),size(Resy,2),size(Xl_load,1),length(zd));
    sigxzs = zeros(size(Resz,1),size(Resz,2),size(Xl_load,1),length(zd));

    epsxs = zeros(size(Resx,1),size(Resx,2),size(Xl_load,1),length(zd));
    epsys = zeros(size(Resy,1),size(Resy,2),size(Xl_load,1),length(zd));
    epszs = zeros(size(Resz,1),size(Resz,2),size(Xl_load,1),length(zd));
    epsxys = zeros(size(Resx,1),size(Resx,2),size(Xl_load,1),length(zd));
    epsyzs = zeros(size(Resy,1),size(Resy,2),size(Xl_load,1),length(zd));
    epsxzs = zeros(size(Resz,1),size(Resz,2),size(Xl_load,1),length(zd));

    for h = 1:length(zd)

        % Gathering deformation points
        Xd = [xd,yd,repelem(zd(h),[length(xd)])'];
        H = alva.zi(end);
        xdef = size(Xd,1);

        for k = 1:size(Xl_load,1)
            alva_Xl = [Xl_load(k,1),Xl_load(k,2)];
            xl = size(alva_Xl,1);
            Xd = repmat(Xd,[xl,1]);
            Xl = repmat(alva_Xl,[xdef,1]);

            if length(a_load) < xl
                a_load = a_load(1)*ones(xl,1);
            end

            drep = repmat(1:xdef,xl,1);
            drep = drep(:)' + repmat(0:xdef:xdef*(xl-1),1,xdef);
            Xd = Xd(drep,:);

            % Distances from load to points in grid
            r  = sqrt((Xd(:,1)-Xl(:,1)).^2 + (Xd(:,2)-Xl(:,2)).^2);

            % Transformation matrix
            T = sparse(3*length(r),3*length(r),(3*3-4)*length(r));

            COS = (Xd(:,1)-Xl(:,1))./r;
            SIN = (Xd(:,2)-Xl(:,2))./r;

            COS(~r) = 1;
            SIN(~r) = 0;

            % Constructing T matrix
            T(1:3*3*length(r)+3:end)               =  COS;
            T(2:3*3*length(r)+3:end)               = -SIN;
            T(3*length(r)+1:3*3*length(r)+3:end)   =  SIN;
            T(3*length(r)+2:3*3*length(r)+3:end)   =  COS;
            T(2*3*length(r)+3:3*3*length(r)+3:end) =  1;

            for j = 1:size(Resx,2)

                % DISPLACEMENTS
                u_rtz = zeros(length(r)*3,1);
                u_rtz(1:3:end) = Resx(:,j,k,h);
                u_rtz(2:3:end) = 0; % Displacements transverse to the radius (u_theta = 0)
                u_rtz(3:3:end) = Resz(:,j,k,h);

                % Transformation
                u_trans = T'*u_rtz;

                % Organize matrix Txyz for adding displacements together
                unos                = zeros(3*xdef*xl,1);
                unos(1:xl)          = 1:9*xdef:9*xl*xdef;
                unos(xl+1:2*xl)     = unos(1:xl)+3*xdef+1;
                unos(2*xl+1:3*xl)   = unos(1:xl)+2*3*xdef+2;

                for g=1:xdef-1
                    unos(3*xl+1+3*xl*(g-1):3*xl+3*xl*g) = unos(1:3*xl)+g*(3*xdef*3*xl+3);
                end

                % Matrix adding contributions together from different loads
                Txyz       = spalloc(3*xdef,3*length(r),length(unos));
                Txyz(unos) = 1;

                % Final response vector
                u1           = Txyz*u_trans;
                uxs(:,j,k,h)    = u1(1:3:end-2);
                uys(:,j,k,h)    = u1(2:3:end-1);
                uzs(:,j,k,h)    = u1(3:3:end);

                % STRESSES
                sig_rtz = sparse(length(r)*3,length(r)*3,(3*3-4)*length(r));
                sig_rtz(1:3*3*length(r)+3:end)               = Resx(:,j,k,h);
                sig_rtz(3*length(r)+2:3*3*length(r)+3:end)   = Resy(:,j,k,h);
                sig_rtz(2*3*length(r)+3:3*3*length(r)+3:end) = Resz(:,j,k,h);
                sig_rtz(3:3*3*length(r)+3:end)               = Resxz(:,j,k,h);
                sig_rtz(2*3*length(r)+1:3*3*length(r)+3:end) = Resxz(:,j,k,h);

                % Transformation
                sig_trans = T'*sig_rtz*T;

                % Transform 3x3 matrices into 6x1 vectors
                dos            = zeros(6*xdef*xl,1);
                dos(1:6:end-5) = 1:9*xdef*xl+3:9*xdef*xl*length(r);
                dos(2:6:end-4) = dos(1:6:end-5) + 3*xdef*xl+1;
                dos(3:6:end-3) = dos(1:6:end-5) + 2*3*xdef*xl+2;
                dos(4:6:end-2) = dos(1:6:end-5) + 1;
                dos(5:6:end-1) = dos(2:6:end-4) + 1;
                dos(6:6:end)   = dos(1:6:end-5) + 2;

                sigv = sig_trans(dos);

                % Organize matrix Mxyz
                tres                = zeros(6*xdef*xl,1);
                tres(1:xl)          = 1:9*4*xdef:9*4*xl*xdef;
                tres(xl+1:2*xl)     = tres(1:xl)+6*xdef+1;
                tres(2*xl+1:3*xl)   = tres(1:xl)+2*6*xdef+2;
                tres(3*xl+1:4*xl)   = tres(1:xl)+3*6*xdef+3;
                tres(4*xl+1:5*xl)   = tres(1:xl)+4*6*xdef+4;
                tres(5*xl+1:6*xl)   = tres(1:xl)+5*6*xdef+5;

                for i=1:xdef-1
                    tres(6*xl+1+6*xl*(i-1):6*xl+6*xl*i)=tres(1:6*xl)+i*(6*xdef*6*xl+6);
                end

                % Matrix adding contributions together from different loads
                Mxyz = spalloc(6*xdef,6*length(r),length(tres));
                Mxyz(tres) = 1;

                % Final response vector
                sig1_final = Mxyz*full(sigv);
                sigxs(:,j,k,h)  = sig1_final(1:6:end-5);
                sigys(:,j,k,h)  = sig1_final(2:6:end-4);
                sigzs(:,j,k,h)  = sig1_final(3:6:end-3);
                sigxys(:,j,k,h) = sig1_final(4:6:end-2);
                sigyzs(:,j,k,h) = sig1_final(5:6:end-1);
                sigxzs(:,j,k,h) = sig1_final(6:6:end);

                % STRAINS
                if exist('sig1_final')
                    rho = r/H;
                    rho(rho==0) = 1e-20;
                    zl = Xd(:,3);
                    % Layer number where the point we use is given
                    layer_no = sum(repmat(zl,1,length(alva.zi)) > repmat(alva.zi,length(rho),1),2)'+1;
                    Nu = alva.nu(layer_no)';
                    Ee = alva.E(layer_no)';
                    Eel = Ee(1:xl:end);
                    Nul = Nu(1:xl:end);

                    % Final response vector
                    epsxs(:,j,k,h) = 1./Eel.*(sigxs(:,j,k,h)-Nul.*(sigys(:,j,k,h)+sigzs(:,j,k,h)));
                    epsys(:,j,k,h) = 1./Eel.*(sigys(:,j,k,h)-Nul.*(sigzs(:,j,k,h)+sigxs(:,j,k,h)));
                    epszs(:,j,k,h) = 1./Eel.*(sigzs(:,j,k,h)-Nul.*(sigxs(:,j,k,h)+sigys(:,j,k,h)));
                    epsxys(:,j,k,h) = (1+Nul)./Eel.*sigxys(:,j,k,h);
                    epsyzs(:,j,k,h) = (1+Nul)./Eel.*sigyzs(:,j,k,h);
                    epsxzs(:,j,k,h) = (1+Nul)./Eel.*sigxzs(:,j,k,h);
                end
            end
        end
    end
    toc

    % -------------------------------------------------------------------------
    % Superpositioning
    % -------------------------------------------------------------------------

    % Performing superposition for all micro circles in the sample patch
    % i.e., for a single patch.

    % Patch in one side of the track
    ux{n} = squeeze(sum(uxs,3)); % [mm]
    uy{n} = squeeze(sum(uys,3)); % [mm]
    uz{n} = squeeze(sum(uzs,3)); % [mm]

    sigx{n} = squeeze(sum(sigxs,3));   % [MPa]
    sigy{n} = squeeze(sum(sigys,3));   % [MPa]
    sigz{n} = squeeze(sum(sigzs,3));   % [MPa]
    sigxy{n} = squeeze(sum(sigxys,3)); % [MPa]
    sigyz{n} = squeeze(sum(sigyzs,3)); % [MPa]
    sigxz{n} = squeeze(sum(sigxzs,3)); % [MPa]

    epsx{n} = squeeze(sum(epsxs,3)).*1e6;   % [micro strain]
    epsy{n} = squeeze(sum(epsys,3)).*1e6;   % [micro strain]
    epsz{n} = squeeze(sum(epszs,3)).*1e6;   % [micro strain]
    epsxy{n} = squeeze(sum(epsxys,3)).*1e6; % [micro strain]
    epsyz{n} = squeeze(sum(epsyzs,3)).*1e6; % [micro strain]
    epsxz{n} = squeeze(sum(epsxzs,3)).*1e6; % [micro strain]

    % Patch on the other side of the track - mirrored

    % Pre-allocating arrays
    ux_fl{n} = zeros(size(Resux,1),size(Resux,2),length(zd));
    uy_fl{n} = zeros(size(Resuy,1),size(Resuy,2),length(zd));
    uz_fl{n} = zeros(size(Resuz,1),size(Resuz,2),length(zd));

    sigx_fl{n} = zeros(size(Resx,1),size(Resx,2),length(zd));
    sigy_fl{n} = zeros(size(Resy,1),size(Resy,2),length(zd));
    sigz_fl{n} = zeros(size(Resz,1),size(Resz,2),length(zd));
    sigxy_fl{n} = zeros(size(Resx,1),size(Resx,2),length(zd));
    sigyz_fl{n} = zeros(size(Resy,1),size(Resy,2),length(zd));
    sigxz_fl{n} = zeros(size(Resz,1),size(Resz,2),length(zd));

    epsx_fl{n} = zeros(size(Resx,1),size(Resx,2),length(zd));
    epsy_fl{n} = zeros(size(Resy,1),size(Resy,2),length(zd));
    epsz_fl{n} = zeros(size(Resz,1),size(Resz,2),length(zd));
    epsxy_fl{n} = zeros(size(Resx,1),size(Resx,2),length(zd));
    epsyz_fl{n} = zeros(size(Resy,1),size(Resy,2),length(zd));
    epsxz_fl{n} = zeros(size(Resz,1),size(Resz,2),length(zd));

    for h = 1:length(zd)
        % [mm]
        ux_fl{n}(:,:,h) = flip(ux{n}(:,:,h),1); % [mm]
        uy_fl{n}(:,:,h) = flip(uy{n}(:,:,h),1); % [mm]
        uz_fl{n}(:,:,h) = flip(uz{n}(:,:,h),1); % [mm]
    end

    for  h = 1:length(zd)
        % [MPa]
        sigx_fl{n}(:,:,h) = flip(sigx{n}(:,:,h),1);   % [MPa]
        sigy_fl{n}(:,:,h) = flip(sigy{n}(:,:,h),1);   % [MPa]
        sigz_fl{n}(:,:,h) = flip(sigz{n}(:,:,h),1);   % [MPa]
        sigxy_fl{n}(:,:,h) = flip(sigxy{n}(:,:,h),1); % [MPa]
        sigyz_fl{n}(:,:,h) = flip(sigyz{n}(:,:,h),1); % [MPa]
        sigxz_fl{n}(:,:,h) = flip(sigxz{n}(:,:,h),1); % [MPa]

        % [micro strain]
        epsx_fl{n}(:,:,h) = flip(epsx{n}(:,:,h),1);   % [micro strain]
        epsy_fl{n}(:,:,h) = flip(epsy{n}(:,:,h),1);   % [micro strain]
        epsz_fl{n}(:,:,h) = flip(epsz{n}(:,:,h),1);   % [micro strain]
        epsxy_fl{n}(:,:,h) = flip(epsxy{n}(:,:,h),1); % [micro strain]
        epsyz_fl{n}(:,:,h) = flip(epsyz{n}(:,:,h),1); % [micro strain]
        epsxz_fl{n}(:,:,h) = flip(epsxz{n}(:,:,h),1); % [micro strain]
    end
end

% -------------------------------------------------------------------------
% Creating a global grid and the full track pattern
% -------------------------------------------------------------------------

% (1a) In this sub-procedure, the micro circles' center coordinates in the
% sample patch are copied and translated in the two dimensions, x & y.

% (1b) Firstly, this is done for a single row of patches in x.

% (1c) Secondly, the single row is copied, mirrored and translated in y
% to represent the second row of patches. (1b) and (1c) combined comprises
% a half vehicle, i.e., one full track and e.g. a full flotation tyre.

% (1d) Lastly, the array of (1c) is again copied and translated in y
% to represent the other half of the vehicle. After this, micro circles
% have been distributed to cover patches for a full two tracks
% (and/or flotation tyres).

% (2a) Once the mapping of micro circles is done, it is calculated where to
% insert all the patches in a global grid.
% The global grid is discretized with the increment size that user chose,
% therefore, all patches must be shifted into the correct positions.

% (2b) When the whole layout has been established in (2a), all the
% responses found for each micro circle in this global grid are then summed
% by the principle of superposition to create a full 2D response plate
% of the vehicle at the user-specific depth, z.

clear Xl_load

% Calculating the offset distances for a rotated patch
offset_conf1 = round((cos(pi/2-theta)*(edge_offset_conf1/cos(pi/2-theta)-L_conf1/2))/inc)*inc;
offset_cent_conf1 = round((cos(pi/2-theta)*(cent_offset_conf1/cos(pi/2-theta)-L_conf1/2))/inc)*inc;
offset_conf2 = round((cos(pi/2-theta)*(edge_offset_conf2/cos(pi/2-theta)-L_conf2/2))/inc)*inc;
offset_cent_conf2 = round((cos(pi/2-theta)*(cent_offset_conf2/cos(pi/2-theta)-L_conf2/2))/inc)*inc;

% Assembling coordinate list for the upper half of the tracked vehicle,
% which comprises the configurations
if L_conf1 == L_conf2 && W_conf1 == W_conf2
    Xl_load(:,:) = Mcirc_conf12{1};
    for i = 1:no_patches_conf12/2
        Xl_load_glob_sn{i} = [Xl_load(:,1)+x_geo(i), Xl_load(:,2)+y_geo(i)];
    end
elseif L_conf1 ~= L_conf2 || W_conf1 ~= W_conf2
    if no_patches_conf1 > 0 && no_patches_conf2 > 0
        Mcirc = cat(1, Mcirc_conf1{1}, Mcirc_conf2{1});
        Xl_load(:,:) = Mcirc;
        Xl_load_conf1(:,:) = Xl_load(1:size(Mcirc_conf1{1},1), 1:2);
        Xl_load_conf2(:,:) = Xl_load(size(Xl_load_conf1,1)+1:end, 1:2);

        for i = 1:no_patches_conf1/2
            Xl_load_glob_sn{i} = [Xl_load_conf1(:,1)+x_geo(i), Xl_load_conf1(:,2)+y_geo(i)];
        end

        for i = no_patches_conf1/2+1:no_patches_conf12/2
            Xl_load_glob_sn{i} = [Xl_load_conf2(:,1)+x_geo(i), Xl_load_conf2(:,2)+y_geo(i)];
        end

    elseif no_patches_conf2 == 0
        Xl_load(:,:) = cat(1, Mcirc_conf1{1});
        Xl_load_conf1(:,:) = Xl_load(1:size(Mcirc_conf1{1},1), 1:2);
        for i = 1:no_patches_conf1/2
            Xl_load_glob_sn{i} = [Xl_load_conf1(:,1)+x_geo(i), Xl_load_conf1(:,2)+y_geo(i)];
        end

    elseif no_patches_conf1 == 0
        Xl_load(:,:) = cat(1, Mcirc_conf2{1});
        Xl_load_conf2(:,:) = Xl_load(1:size(Mcirc_conf2{1},1), 1:2);
        for i = no_patches_conf1/2+1:no_patches_conf12/2
            Xl_load_glob_sn{i} = [Xl_load_conf2(:,1)+x_geo(i), Xl_load_conf2(:,2)+y_geo(i)];
        end
    end
end

Xl_load_glob_sn = vertcat(Xl_load_glob_sn{:});

% Assembling coordinate list for the other half of the tracked vehicle
Xl_load_glob_bn = Xl_load_glob_sn;

% Counting amount of micro circles in each configuration
if L_conf1 == L_conf2 && W_conf1 == W_conf2
    circ_conf12 = no_patches_conf12*length(Mcirc_conf12{1});
    circ_conf12_1 = no_patches_conf1*length(Mcirc_conf12{1});
    circ_conf12_2 = no_patches_conf2*length(Mcirc_conf12{1});
elseif L_conf1 ~= L_conf2 || W_conf1 ~= W_conf2
    if no_patches_conf1 > 0 && no_patches_conf2 > 0
        circ_conf1 = no_patches_conf1*length(Mcirc_conf1{1});
        circ_conf2 = no_patches_conf2*length(Mcirc_conf2{1});
        circ_conf12 = circ_conf1+circ_conf2;
    elseif no_patches_conf2 == 0
        circ_conf1 = no_patches_conf1*length(Mcirc_conf1{1});
        circ_conf12 = circ_conf1;
    elseif no_patches_conf1 == 0
        circ_conf2 = no_patches_conf2*length(Mcirc_conf2{1});
        circ_conf12 = circ_conf2;
    end

end

% Translating micro circles in x & y dimensions
if L_conf1 == L_conf2 && W_conf1 == W_conf2
    for i = 1:circ_conf12_1/2
        Xl_load_glob_bn(i,2) = Xl_load_glob_sn(i,2)-(L_conf1+2*offset_cent_conf1);
        Xl_load_glob_bn(i,1) = Xl_load_glob_bn(i,1)+pitch_conf1/2;
    end
    for i = circ_conf12_1/2+1:circ_conf12_2/2
        Xl_load_glob_bn(i,2) = Xl_load_glob_sn(i,2)-(L_conf2+2*offset_cent_conf2);
        Xl_load_glob_bn(i,1) = Xl_load_glob_bn(i,1)+pitch_conf2/2;
    end

elseif L_conf1 ~= L_conf2 || W_conf1 ~= W_conf2
    if no_patches_conf1 > 0 && no_patches_conf2 > 0
        for i = 1:circ_conf1/2
            Xl_load_glob_bn(i,2) = Xl_load_glob_sn(i,2)-(L_conf1+2*offset_cent_conf1);
            Xl_load_glob_bn(i,1) = Xl_load_glob_bn(i,1)+pitch_conf1/2;
        end
        for i = circ_conf1/2+1:circ_conf2/2
            Xl_load_glob_bn(i,2) = Xl_load_glob_sn(i,2)-(L_conf2+2*offset_cent_conf2);
            Xl_load_glob_bn(i,1) = Xl_load_glob_bn(i,1)+pitch_conf2/2;
        end

    elseif no_patches_conf2 == 0
        for i = 1:circ_conf1/2
            Xl_load_glob_bn(i,2) = Xl_load_glob_sn(i,2)-(L_conf1+2*offset_cent_conf1);
            Xl_load_glob_bn(i,1) = Xl_load_glob_bn(i,1)+pitch_conf1/2;
        end

    elseif no_patches_conf1 == 0
        for i = 1:circ_conf2/2
            Xl_load_glob_bn(i,2) = Xl_load_glob_sn(i,2)-(L_conf2+2*offset_cent_conf2);
            Xl_load_glob_bn(i,1) = Xl_load_glob_bn(i,1)+pitch_conf2/2;
        end
    end
end

% Coordinate list for micro circles for both upper and lower patch rows
Xl_load_glob_bn(:,2) = flip(Xl_load_glob_bn(:,2));
Xl_load_glob_half = [Xl_load_glob_sn; Xl_load_glob_bn];
Xl_load_glob = [Xl_load_glob_half; Xl_load_glob_half];

% Constructing second half of tracked vehicle by translation in the y dimension
if L_conf1 == L_conf2 && W_conf1 == W_conf2
    for i = [size(Xl_load_glob_half,1)+1:size(Xl_load_glob_half,1)+circ_conf12_1/2,...
            size(Xl_load_glob_half,1)+circ_conf12/2+1:size(Xl_load_glob_half,1)+circ_conf12/2+circ_conf12_1/2]
        Xl_load_glob(i,2) = Xl_load_glob(i,2)-(dist_conf1-(2*L_conf1+2*offset_cent_conf1+offset_conf1));
    end
    for i = [size(Xl_load_glob_half,1)+circ_conf12_1/2+1:size(Xl_load_glob_half,1)+circ_conf12/2,...
            size(Xl_load_glob_half,1)+circ_conf12/2+circ_conf12_1/2+1:Xl_load_glob]
        Xl_load_glob(i,2) = Xl_load_glob(i,2)-(dist_conf2-(2*L_conf2+2*offset_cent_conf2+offset_conf2));
    end

elseif L_conf1 ~= L_conf2 || W_conf1 ~= W_conf2
    if no_patches_conf1 > 0 && no_patches_conf2 > 0
        for i = [size(Xl_load_glob_half,1)+1:size(Xl_load_glob_half,1)+circ_conf1/2,...
                size(Xl_load_glob_half,1)+circ_conf12/2+1:size(Xl_load_glob_half,1)+circ_conf12/2+circ_conf1/2]
            Xl_load_glob(i,2) = Xl_load_glob(i,2)-(dist_conf1-(2*L_conf1+2*offset_cent_conf1+offset_conf1));
        end
        for i = [size(Xl_load_glob_half,1)+circ_conf1/2+1:size(Xl_load_glob_half,1)+circ_conf12/2,...
                size(Xl_load_glob_half,1)+circ_conf12/2+circ_conf1/2+1:Xl_load_glob]
            Xl_load_glob(i,2) = Xl_load_glob(i,2)-(dist_conf2-(2*L_conf2+2*offset_cent_conf2+offset_conf2));
        end

    elseif no_patches_conf2 == 0
        for i = [size(Xl_load_glob_half,1)+1:size(Xl_load_glob_half,1)+circ_conf1/2,...
                size(Xl_load_glob_half,1)+circ_conf12/2+1:size(Xl_load_glob_half,1)+circ_conf12/2+circ_conf1/2]
            Xl_load_glob(i,2) = Xl_load_glob(i,2)-(dist_conf1-(2*L_conf1+2*offset_cent_conf1+offset_conf1));
        end

    elseif no_patches_conf1 == 0
        for i = [size(Xl_load_glob_half,1)+1:size(Xl_load_glob_half,1)+circ_conf2/2,...
                size(Xl_load_glob_half,1)+circ_conf12/2+1:size(Xl_load_glob_half,1)+circ_conf12/2+circ_conf2/2]
            Xl_load_glob(i,2) = Xl_load_glob(i,2)-(dist_conf2-(2*L_conf2+2*offset_cent_conf2+offset_conf2));
        end
    end
end

% Locating the maximum coordinates (endings) in both X & Y directions

if no_patches_conf1 > 0 && no_patches_conf2 > 0
    rot_end1 = cos(theta)*L_conf1; % Because the patches are rotated, their ends
    rot_end2 = cos(theta)*L_conf2; % protrude further than the initial grid boundary
elseif no_patches_conf2 == 0
    rot_end1 = cos(theta)*L_conf1;
    rot_end2 = rot_end1; 
elseif no_patches_conf1 == 0
    rot_end1 = cos(theta)*L_conf2;
    rot_end2 = rot_end1; 
end
x_ends = [min(Xl_load_glob(:,1))+min(xd)-rot_end1, max(Xl_load_glob(:,1))+max(xd)+rot_end2];


offset_conf = max(offset_conf1, offset_conf2);
if abs(min(Xl_load_glob(:,2))+min(yd)) > max(Xl_load_glob(:,2))+max(yd)
    y_ends = [min(Xl_load_glob(:,2))+min(yd)-offset_conf, -(min(Xl_load_glob(:,2))+min(yd)-offset_conf)];
else
    y_ends = [-(max(Xl_load_glob(:,2))+max(yd)+offset_conf), max(Xl_load_glob(:,2))+max(yd)+offset_conf];
end
    
% Discretizing the global grid
if min(xd) < min(x_ends)
    xd_glob = xd(1):dx:ceil(max(x_ends)/dx)*dx;
elseif min(xd) > min(x_ends)
    xd_glob = floor(min(x_ends)/dx)*dx:dx:ceil(max(x_ends)/dx)*dx;
end

if min(yd) < min(y_ends)
    if max(yd) > max(y_ends)
        yd_glob = yd;
    else
        yd_glob = yd(1):dy:ceil(max(y_ends)/dy)*dy;
    end
elseif min(yd) > min(y_ends)
    if max(yd) > max(y_ends)
        yd_glob = floor(min(y_ends)/dy)*dy:dy:yd(end);
    else
        yd_glob = floor(min(y_ends)/dy)*dy:dy:ceil(max(y_ends)/dy)*dy;
    end
end

% Meshing global coordinates
[X_glob,Y_glob] = meshgrid(xd_glob,yd_glob);


[row_glob,  col_glob]  = size(X_glob); % Size of global grid
[row_input, col_input] = size(X); % Size of response grid used as input

% Positioning response matrices into global matrix
y_start = zeros(1,no_patches);
y_stop = zeros(1,no_patches);
x_start = zeros(1,no_patches);
x_stop = zeros(1,no_patches);

if no_patches_conf1 > 0 && no_patches_conf2 > 0
y_geo_glob = [y_geo, y_geo(1:no_patches_conf1/2)-(L_conf1+2*offset_cent_conf1),...
              y_geo(no_patches_conf1/2+1:end)-(L_conf2+2*offset_cent_conf2),...
             -y_geo(1:no_patches_conf1/2)+(L_conf1+2*offset_cent_conf1),...
             -y_geo(no_patches_conf1/2+1:end)+(L_conf2+2*offset_cent_conf2), -y_geo];

elseif no_patches_conf2 == 0
y_geo_glob = [y_geo, y_geo(1:no_patches_conf1/2)-(L_conf1+2*offset_cent_conf1),...
             -y_geo(1:no_patches_conf1/2)+(L_conf1+2*offset_cent_conf1), -y_geo];

elseif no_patches_conf1 == 0
y_geo_glob = [y_geo, y_geo(1:no_patches_conf2/2)-(L_conf2+2*offset_cent_conf2),...
             -y_geo(1:no_patches_conf2/2)+(L_conf2+2*offset_cent_conf2), -y_geo];
end
    
if pitch_conf1 == pitch_conf2
    if no_patches_conf1 == 0
        x_geo_glob = [x_geo, x_geo+pitch_conf2/2, x_geo+pitch_conf2/2, x_geo];

    else
        x_geo_glob = [x_geo, x_geo+pitch_conf1/2, x_geo+pitch_conf1/2, x_geo];
    end

elseif pitch_conf1 ~= pitch_conf2
    if no_patches_conf1 > 0 && no_patches_conf2 > 0
        x_geo_glob = [x_geo, x_geo(1:no_patches_conf1/2)+pitch_conf1/2,...
                      x_geo(no_patches_conf1/2+1:end)+pitch_conf2/2,...
                      x_geo(1:no_patches_conf1/2)+pitch_conf1/2,...
                      x_geo(no_patches_conf1/2+1:end)+pitch_conf2/2,...
                      x_geo];

    elseif no_patches_conf2 == 0
        x_geo_glob = [x_geo, x_geo(1:no_patches_conf1/2)+pitch_conf1/2,...
                      x_geo(1:no_patches_conf1/2)+pitch_conf1/2,...
                      x_geo];

    elseif no_patches_conf1 == 0
        x_geo_glob = [x_geo, x_geo(1:no_patches_conf2/2)+pitch_conf2/2,...
                      x_geo(1:no_patches_conf2/2)+pitch_conf2/2,...
                      x_geo];
    end
end
    
for i = 1:no_patches
    y_start(1,i) = find(Y(1,1)+y_geo_glob(i)   == yd_glob); % Start row
    y_stop(1,i)  = find(Y(end,1)+y_geo_glob(i) == yd_glob); % End Row
    x_start(1,i) = find(X(1,1)+x_geo_glob(i)   == xd_glob); % Start column
    x_stop(1,i)  = find(X(1,end)+x_geo_glob(i) == xd_glob); % End column
end

% Assembling the global grid and superpositioning all patches
for n = 1:n_confs
    if L_conf1 == L_conf2 && W_conf1 == W_conf2
        tic
        np_conf = no_patches_conf12;
        fprintf('\nAssembling Full Grid of Responses for Configurations 1 & 2\n')

        run('grid_assembly.m') 
        toc

    elseif L_conf1 ~= L_conf2 || W_conf1 ~= W_conf2
        if n == 1
            if no_patches_conf1 > 0
                tic
                np_conf = no_patches_conf1;
                fprintf('\nAssembling Full Grid of Responses for Configuration 1\n')

            elseif no_patches_conf1 == 0
                tic
                np_conf = no_patches_conf2;
                fprintf('\nAssembling Full Grid of Responses for Configuration 2\n')
            end

            run('grid_assembly.m')

            ux_glob_n1 = ux_glob;
            uy_glob_n1 = uy_glob;
            uz_glob_n1 = uz_glob;

            sigx_glob_n1  = sigx_glob;
            sigy_glob_n1  = sigy_glob;
            sigz_glob_n1  = sigz_glob;
            sigxy_glob_n1 = sigxy_glob;
            sigyz_glob_n1 = sigyz_glob;
            sigxz_glob_n1 = sigxz_glob;

            epsx_glob_n1  = epsx_glob;
            epsy_glob_n1  = epsy_glob;
            epsz_glob_n1  = epsz_glob;
            epsxy_glob_n1 = epsxy_glob;
            epsyz_glob_n1 = epsyz_glob;
            epsxz_glob_n1 = epsxz_glob;
            toc

        elseif n == 2
            tic
            np_conf = no_patches_conf2;
            fprintf('\nAssembling Full Grid of Responses for Configuration 2\n')
            run('grid_assembly.m')

            ux_glob_n2 = ux_glob;
            uy_glob_n2 = uy_glob;
            uz_glob_n2 = uz_glob;

            sigx_glob_n2  = sigx_glob;
            sigy_glob_n2  = sigy_glob;
            sigz_glob_n2  = sigz_glob;
            sigxy_glob_n2 = sigxy_glob;
            sigyz_glob_n2 = sigyz_glob;
            sigxz_glob_n2 = sigxz_glob;

            epsx_glob_n2  = epsx_glob;
            epsy_glob_n2  = epsy_glob;
            epsz_glob_n2  = epsz_glob;
            epsxy_glob_n2 = epsxy_glob;
            epsyz_glob_n2 = epsyz_glob;
            epsxz_glob_n2 = epsxz_glob;
            toc
        
            % Superpositioning the final grid for configurations 1 & 2
            ux_glob = ux_glob_n1 + ux_glob_n2;
            uy_glob = uy_glob_n1 + uy_glob_n2;
            uz_glob = uz_glob_n1 + uz_glob_n2;

            sigx_glob  = sigx_glob_n1+sigx_glob_n2;
            sigy_glob  = sigy_glob_n1+sigy_glob_n2;
            sigz_glob  = sigz_glob_n1+sigz_glob_n2;
            sigxy_glob = sigxy_glob_n1+sigxy_glob_n2;
            sigyz_glob = sigyz_glob_n1+sigyz_glob_n2;
            sigxz_glob = sigxz_glob_n1+sigxz_glob_n2;

            epsx_glob  = epsx_glob_n1+epsx_glob_n2;
            epsy_glob  = epsy_glob_n1+epsy_glob_n2;
            epsz_glob  = epsz_glob_n1+epsz_glob_n2;
            epsxy_glob = epsxy_glob_n1+epsxy_glob_n2;
            epsyz_glob = epsyz_glob_n1+epsyz_glob_n2;
            epsxz_glob = epsxz_glob_n1+epsxz_glob_n2;
        end
    end
end
