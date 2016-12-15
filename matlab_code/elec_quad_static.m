%% Calculate electric field for static electric dipole

% Determine directory where script is located 
% (not necessarily current matlab directory)
% to determine where to save results
[scriptdir,~,~] = fileparts(mfilename('fullpath'));
resultsdir = [scriptdir,'/../results/'];
imgdir = [resultsdir,'static_img/'];

% Add required packages to matlab path
% 'genpath' adds subfolders as well
addpath(genpath([scriptdir,'/../required_packages/']))

% nxm Cartesian grid
nn = 1024;
mm = 1024;

% Cartesian mesh
xmax = nn/1000;
xmin = -xmax;
ymax = mm/1000;
ymin = -ymax;
x_arr = linspace(xmin,xmax,nn);
y_arr = linspace(ymin,ymax,mm);
[x_mesh,y_mesh] = meshgrid(x_arr,y_arr);

% Polar mesh
[th_mesh,r_mesh] = cart2pol(x_mesh,y_mesh);

% Parameters
QQ = 1e-3;
kk = 1e1;

% Create noise canvas for LIC
MM = randn([nn,mm]);

% Magnitude matrix
% Dimensions: ii grid, jj grid
mag_mat = zeros(nn,mm);

% Vector matrix
% Dimensions: ii grid, jj grid, comp # (x,y)
vec_mat = zeros(nn,mm,2);

%%%%%%%%%%%%%%%%%%%%%
%% Calculate Field %%
%%%%%%%%%%%%%%%%%%%%%

% Loop through space
for ii=1:nn
    for jj=1:mm
        % Shorthand for r and theta for this grid point
        rr = r_mesh(ii,jj);
        th = th_mesh(ii,jj);

        % Spherical components of electric field
        Er = 6*QQ/rr^4*(3/2*cos(th)^2-.5);
        Eth = 6*QQ/rr^4*cos(th)*sin(th);

        % Combine into one vector
        vec_sph = [Er;Eth];

        % Convert to cartesian
        vecs_cart(1) = vec_sph(1)*sin(th) + vec_sph(2)*cos(th);
        vecs_cart(2) = vec_sph(1)*cos(th) - vec_sph(2)*sin(th);
        
        % Calculate vector norm
        mag = sqrt(vecs_cart(1)^2+vecs_cart(2)^2);
        
        % Save to appropriate matrices
        mag_mat(ii,jj) = mag;
        vec_mat(ii,jj,:) = vecs_cart(:);
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Line Integral Convolution %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Normalize vector field
vec_mat = perform_vf_normalization(vec_mat);

% LIC Parameters
options.bound = 'sym';
options.histogram = 'linear'; % keep contrast fixed
options.verb = 1;
options.dt = 1.5; % time steping
options.flow_correction = 1;
options.isoriented=0;
options.niter_lic = 3; % several iterations gives better results
options.M0 = MM;
LL = 50; % "Smear length"

% Perform LIC
lic_out = perform_lic(vec_mat, LL, options);

% Use constant-brightness colormap
cmap = cmocean('phase');

% Use log scale for colors
mag_mat = log(abs(mag_mat));

% Determine colors over magnitude range
mmax = max(mag_mat(:));
mmin = min(mag_mat(:));

% Handle out of bounds colors by clipping them
mag_mat(mag_mat<mmin) = mmin;
mag_mat(mag_mat>mmax) = mmax;

% Determine colors from magnitude
colors = floor((mag_mat(:)-mmin)/(mmax-mmin).*size(cmap,1));

% Set vectors with the smallest magnitude to the first color
% (as opposed to the zeroth color)
colors(colors==0) = 1;

% Multiply b/w img by colors & reshape to n x m x 3
img = lic_out.*reshape(cmap(colors,:),[size(lic_out) 3]);

% Plot final colored image
%imshow(img)
%drawnow

% Save image to file
imwrite(img,[imgdir,'elec_static.png'])

disp 'Finished!'

 
