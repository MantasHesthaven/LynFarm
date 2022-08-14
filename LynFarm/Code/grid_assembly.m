%% Displacements

% Assigning empty matrices
M_ux = zeros(row_glob, col_glob, no_patches, length(zd));
M_uy = zeros(row_glob, col_glob, no_patches, length(zd));
M_uz = zeros(row_glob, col_glob, no_patches, length(zd));

q = q_all;
count_col = 0; count_row = 0;
np = no_patches_conf12;
np2 = np*1.5;

for h = 1:length(zd)
    if n == 1
        nx = 0;
    elseif n == 2
        nx = no_patches_conf1/2;
    end
    for i = nx+1:np_conf/2+nx
        % Checking if global matrix can fit
        if x_stop(i) <= col_glob
            if y_stop(i) <= row_glob

                % If fits, paste displacements
                M_ux(y_start(i):y_stop(i), x_start(i):x_stop(i),i,h) = q(i)/alva.q*ux{n}(:,:,h);
                M_ux(y_start(i+np):y_stop(i+np), x_start(i+np):x_stop(i+np),i+np,h) = q(i+np)/alva.q*ux{n}(:,:,h);
                M_ux(y_start(i+np/2):y_stop(i+np/2), x_start(i+np/2):x_stop(i+np/2),i+np/2,h) = q(i+np/2)/alva.q*ux_fl{n}(:,:,h);
                M_ux(y_start(i+np2):y_stop(i+np2), x_start(i+np2):x_stop(i+np2),i+np2,h) = q(i+np2)/alva.q*ux_fl{n}(:,:,h);

                M_uy(y_start(i):y_stop(i), x_start(i):x_stop(i),i,h) = q(i)/alva.q*uy{n}(:,:,h);
                M_uy(y_start(i+np):y_stop(i+np), x_start(i+np):x_stop(i+np),i+np,h) = q(i+np)/alva.q*uy{n}(:,:,h);
                M_uy(y_start(i+np/2):y_stop(i+np/2), x_start(i+np/2):x_stop(i+np/2),i+np/2,h) = q(i+np/2)/alva.q*uy_fl{n}(:,:,h);
                M_uy(y_start(i+np2):y_stop(i+np2), x_start(i+np2):x_stop(i+np2),i+np2,h) = q(i+np2)/alva.q*uy_fl{n}(:,:,h);

                M_uz(y_start(i):y_stop(i), x_start(i):x_stop(i),i,h) = q(i)/alva.q*uz{n}(:,:,h);
                M_uz(y_start(i+np):y_stop(i+np), x_start(i+np):x_stop(i+np),i+np,h) = q(i+np)/alva.q*uz{n}(:,:,h);
                M_uz(y_start(i+np/2):y_stop(i+np/2), x_start(i+np/2):x_stop(i+np/2),i+np/2,h) = q(i+np/2)/alva.q*uz_fl{n}(:,:,h);
                M_uz(y_start(i+np2):y_stop(i+np2), x_start(i+np2):x_stop(i+np2),i+np2,h) = q(i+np2)/alva.q*uz_fl{n}(:,:,h);

            else
                % Error if no fit
                count_row = i;
                err = "rows";
            end
        else
            % Error if no fit
            count_col = i;
            err = "cols";
        end
    end
end
% Superposition in global grid
ux_glob = squeeze(sum(M_ux,3));
uy_glob = squeeze(sum(M_uy,3));
uz_glob = squeeze(sum(M_uz,3));

%% Stresses and Strains

% Assigning empty matrices
M_sigx = zeros(row_glob, col_glob, no_patches, length(zd));
M_sigy = zeros(row_glob, col_glob, no_patches, length(zd));
M_sigz = zeros(row_glob, col_glob, no_patches, length(zd));
M_sigxy = zeros(row_glob, col_glob, no_patches, length(zd));
M_sigyz = zeros(row_glob, col_glob, no_patches, length(zd));
M_sigxz = zeros(row_glob, col_glob, no_patches, length(zd));

M_epsx = zeros(row_glob, col_glob, no_patches, length(zd));
M_epsy = zeros(row_glob, col_glob, no_patches, length(zd));
M_epsz = zeros(row_glob, col_glob, no_patches, length(zd));
M_epsxy = zeros(row_glob, col_glob, no_patches, length(zd));
M_epsyz = zeros(row_glob, col_glob, no_patches, length(zd));
M_epsxz = zeros(row_glob, col_glob, no_patches, length(zd));

count_col = 0; count_row = 0;

for h = 1:length(zd)
    if n == 1
        nx = 0;
    elseif n == 2
        nx = no_patches_conf1/2;
    end
    for i = nx+1:np_conf/2+nx
        % Checking if global matrix can fit
        if x_stop(i) <= col_glob
            if y_stop(i) <= row_glob

                % If fits, paste stresses
                M_sigx(y_start(i):y_stop(i), x_start(i):x_stop(i),i,h) = q(i)/alva.q*sigx{n}(:,:,h);
                M_sigx(y_start(i+np):y_stop(i+np), x_start(i+np):x_stop(i+np),i+np,h) = q(i+np)/alva.q*sigx{n}(:,:,h);
                M_sigx(y_start(i+np/2):y_stop(i+np/2), x_start(i+np/2):x_stop(i+np/2),i+np/2,h) = q(i+np/2)/alva.q*sigx_fl{n}(:,:,h);
                M_sigx(y_start(i+np2):y_stop(i+np2), x_start(i+np2):x_stop(i+np2),i+np2,h) = q(i+np2)/alva.q*sigx_fl{n}(:,:,h);

                M_sigy(y_start(i):y_stop(i), x_start(i):x_stop(i),i,h) = q(i)/alva.q*sigy{n}(:,:,h);
                M_sigy(y_start(i+np):y_stop(i+np), x_start(i+np):x_stop(i+np),i+np,h) = q(i+np)/alva.q*sigy{n}(:,:,h);
                M_sigy(y_start(i+np/2):y_stop(i+np/2), x_start(i+np/2):x_stop(i+np/2),i+np/2,h) = q(i+np/2)/alva.q*sigy_fl{n}(:,:,h);
                M_sigy(y_start(i+np2):y_stop(i+np2), x_start(i+np2):x_stop(i+np2),i+np2,h) = q(i+np2)/alva.q*sigy_fl{n}(:,:,h);

                M_sigz(y_start(i):y_stop(i), x_start(i):x_stop(i),i,h) = q(i)/alva.q*sigz{n}(:,:,h);
                M_sigz(y_start(i+np):y_stop(i+np), x_start(i+np):x_stop(i+np),i+np,h) = q(i+np)/alva.q*sigz{n}(:,:,h);
                M_sigz(y_start(i+np/2):y_stop(i+np/2), x_start(i+np/2):x_stop(i+np/2),i+np/2,h) = q(i+np/2)/alva.q*sigz_fl{n}(:,:,h);
                M_sigz(y_start(i+np2):y_stop(i+np2), x_start(i+np2):x_stop(i+np2),i+np2,h) = q(i+np2)/alva.q*sigz_fl{n}(:,:,h);

                M_sigxy(y_start(i):y_stop(i), x_start(i):x_stop(i),i,h) = q(i)/alva.q*sigxy{n}(:,:,h);
                M_sigxy(y_start(i+np):y_stop(i+np), x_start(i+np):x_stop(i+np),i+np,h) = q(i+np)/alva.q*sigxy{n}(:,:,h);
                M_sigxy(y_start(i+np/2):y_stop(i+np/2), x_start(i+np/2):x_stop(i+np/2),i+np/2,h) = q(i+np/2)/alva.q*sigxy_fl{n}(:,:,h);
                M_sigxy(y_start(i+np2):y_stop(i+np2), x_start(i+np2):x_stop(i+np2),i+np2,h) = q(i+np2)/alva.q*sigxy_fl{n}(:,:,h);

                M_sigyz(y_start(i):y_stop(i), x_start(i):x_stop(i),i,h) = q(i)/alva.q*sigyz{n}(:,:,h);
                M_sigyz(y_start(i+np):y_stop(i+np), x_start(i+np):x_stop(i+np),i+np,h) = q(i+np)/alva.q*sigyz{n}(:,:,h);
                M_sigyz(y_start(i+np/2):y_stop(i+np/2), x_start(i+np/2):x_stop(i+np/2),i+np/2,h) = q(i+np/2)/alva.q*sigyz_fl{n}(:,:,h);
                M_sigyz(y_start(i+np2):y_stop(i+np2), x_start(i+np2):x_stop(i+np2),i+np2,h) = q(i+np2)/alva.q*sigyz_fl{n}(:,:,h);

                M_sigxz(y_start(i):y_stop(i), x_start(i):x_stop(i),i,h) = q(i)/alva.q*sigxz{n}(:,:,h);
                M_sigxz(y_start(i+np):y_stop(i+np), x_start(i+np):x_stop(i+np),i+np,h) = q(i+np)/alva.q*sigxz{n}(:,:,h);
                M_sigxz(y_start(i+np/2):y_stop(i+np/2), x_start(i+np/2):x_stop(i+np/2),i+np/2,h) = q(i+np/2)/alva.q*sigxz_fl{n}(:,:,h);
                M_sigxz(y_start(i+np2):y_stop(i+np2), x_start(i+np2):x_stop(i+np2),i+np2,h) = q(i+np2)/alva.q*sigxz_fl{n}(:,:,h);

                % --------------- Paste Strains ---------------
                M_epsx(y_start(i):y_stop(i), x_start(i):x_stop(i),i,h) = q(i)/alva.q*epsx{n}(:,:,h);
                M_epsx(y_start(i+np):y_stop(i+np), x_start(i+np):x_stop(i+np),i+np,h) = q(i+np)/alva.q*epsx{n}(:,:,h);
                M_epsx(y_start(i+np/2):y_stop(i+np/2), x_start(i+np/2):x_stop(i+np/2),i+np/2,h) = q(i+np/2)/alva.q*epsx_fl{n}(:,:,h);
                M_epsx(y_start(i+np2):y_stop(i+np2), x_start(i+np2):x_stop(i+np2),i+np2,h) = q(i+np2)/alva.q*epsx_fl{n}(:,:,h);

                M_epsy(y_start(i):y_stop(i), x_start(i):x_stop(i),i,h) = q(i)/alva.q*epsy{n}(:,:,h);
                M_epsy(y_start(i+np):y_stop(i+np), x_start(i+np):x_stop(i+np),i+np,h) = q(i+np)/alva.q*epsy{n}(:,:,h);
                M_epsy(y_start(i+np/2):y_stop(i+np/2), x_start(i+np/2):x_stop(i+np/2),i+np/2,h) = q(i+np/2)/alva.q*epsy_fl{n}(:,:,h);
                M_epsy(y_start(i+np2):y_stop(i+np2), x_start(i+np2):x_stop(i+np2),i+np2,h) = q(i+np2)/alva.q*epsy_fl{n}(:,:,h);

                M_epsz(y_start(i):y_stop(i), x_start(i):x_stop(i),i,h) = q(i)/alva.q*epsz{n}(:,:,h);
                M_epsz(y_start(i+np):y_stop(i+np), x_start(i+np):x_stop(i+np),i+np,h) = q(i+np)/alva.q*epsz{n}(:,:,h);
                M_epsz(y_start(i+np/2):y_stop(i+np/2), x_start(i+np/2):x_stop(i+np/2),i+np/2,h) = q(i+np/2)/alva.q*epsz_fl{n}(:,:,h);
                M_epsz(y_start(i+np2):y_stop(i+np2), x_start(i+np2):x_stop(i+np2),i+np2,h) = q(i+np2)/alva.q*epsz_fl{n}(:,:,h);

                M_epsxy(y_start(i):y_stop(i), x_start(i):x_stop(i),i,h) = q(i)/alva.q*epsxy{n}(:,:,h);
                M_epsxy(y_start(i+np):y_stop(i+np), x_start(i+np):x_stop(i+np),i+np,h) = q(i+np)/alva.q*epsxy{n}(:,:,h);
                M_epsxy(y_start(i+np/2):y_stop(i+np/2), x_start(i+np/2):x_stop(i+np/2),i+np/2,h) = q(i+np/2)/alva.q*epsxy_fl{n}(:,:,h);
                M_epsxy(y_start(i+np2):y_stop(i+np2), x_start(i+np2):x_stop(i+np2),i+np2,h) = q(i+np2)/alva.q*epsxy_fl{n}(:,:,h);

                M_epsyz(y_start(i):y_stop(i), x_start(i):x_stop(i),i,h) = q(i)/alva.q*epsyz{n}(:,:,h);
                M_epsyz(y_start(i+np):y_stop(i+np), x_start(i+np):x_stop(i+np),i+np,h) = q(i+np)/alva.q*epsyz{n}(:,:,h);
                M_epsyz(y_start(i+np/2):y_stop(i+np/2), x_start(i+np/2):x_stop(i+np/2),i+np/2,h) = q(i+np/2)/alva.q*epsyz_fl{n}(:,:,h);
                M_epsyz(y_start(i+np2):y_stop(i+np2), x_start(i+np2):x_stop(i+np2),i+np2,h) = q(i+np2)/alva.q*epsyz_fl{n}(:,:,h);

                M_epsxz(y_start(i):y_stop(i), x_start(i):x_stop(i),i,h) = q(i)/alva.q*epsxz{n}(:,:,h);
                M_epsxz(y_start(i+np):y_stop(i+np), x_start(i+np):x_stop(i+np),i+np,h) = q(i+np)/alva.q*epsxz{n}(:,:,h);
                M_epsxz(y_start(i+np/2):y_stop(i+np/2), x_start(i+np/2):x_stop(i+np/2),i+np/2,h) = q(i+np/2)/alva.q*epsxz_fl{n}(:,:,h);
                M_epsxz(y_start(i+np2):y_stop(i+np2), x_start(i+np2):x_stop(i+np2),i+np2,h) = q(i+np2)/alva.q*epsxz_fl{n}(:,:,h);

            else
                % Error if no fit
                count_row = i;
                err = "rows";
            end
        else
            % Error if no fit
            count_col = i;
            err = "cols";
        end
    end
end
% Superposition in global grid
% Stresses
sigx_glob = squeeze(sum(M_sigx,3));
sigy_glob = squeeze(sum(M_sigy,3));
sigz_glob = squeeze(sum(M_sigz,3));
sigxy_glob = squeeze(sum(M_sigxy,3));
sigyz_glob = squeeze(sum(M_sigyz,3));
sigxz_glob = squeeze(sum(M_sigxz,3));

epsx_glob = squeeze(sum(M_epsx,3));
epsy_glob = squeeze(sum(M_epsy,3));
epsz_glob = squeeze(sum(M_epsz,3));
epsxy_glob = squeeze(sum(M_epsxy,3));
epsyz_glob = squeeze(sum(M_epsyz,3));
epsxz_glob = squeeze(sum(M_epsxz,3));


% Check for error criteria
if exist('err')
    if err == "rows"
        fprintf('\nPatch (%.0f) does not fit in global matrix. Rows extended.\n', count_row);
    elseif err == "cols"
        fprintf('\nPatch (%.0f) does not fit in global matrix. Columns extended.\n', count_col);
    end
end
