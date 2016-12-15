%% Calculate electric field for oscillating electric dipole

% Add required packages to matlab path
% 'genpath' adds subfolders as well
addpath(genpath('required_packages/cmocean'))
addpath(genpath('required_packages/toolbox_image'))

% Determine directory where script is located 
% (not necessarily current matlab directory)
% to determine where to save results
[scriptdir,~,~] = fileparts(mfilename('fullpath'));
resultsdir = [scriptdir,'/../results/elec_img'];

% nxm Cartesian grid
nn = 1000;
mm = 1000;

% Cartesian mesh
start = 1;
xg = linspace(-start,start,nn);
yg = linspace(-start,start,mm);
[xx,yy] = meshgrid(xg,yg);

% Polar mesh
[tt,rr] = cart2pol(xx,yy);

% Parameters
Q = 1e-3;
w = 2*pi;
k = 1e1;

% Create noise canvas for LIC
M = randn([nn,mm]);

% Loop through time
nsteps = 50;
tmax = 1;
dt = tmax/nsteps;
tstep = 1;
for t=0:dt:tmax
    % Magnitude matrix
    % Dimensions: ii grid, jj grid
    mag_mat = zeros(nn,mm);
    
    % Vector matrix
    % Dimensions: ii grid, jj grid, comp # (x,y)
    vec_mat = zeros(nn,mm,2);
    
    % Print current timestep
    fprintf('tstep = %d/%d\n',tstep,nsteps)
    
    % Loop through space
    for ii=1:nn
        for jj=1:mm
            % Shorthand for r and theta for this grid point
            r = rr(ii,jj);
            th = tt(ii,jj);

            % Spherical components of electric field
            Er = -2*Q*k^2/r^2 * (cos(k*r - w*t)*(1-3/(k*r)^2)... 
                - 3/(k*r)*sin(k*r-w*t)) * (3/2*cos(th)^2-.5);
            Eth = -Q*k^2/r^2 * (sin(k*r - w*t)*(k*r - 6/(k*r)) ...
                + (3-6/(k*r)^2)*cos(k*r-w*t))*cos(th)*sin(th);

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
    options.isoriented=1;
    options.niter_lic = 2; % several iterations gives better results
    options.M0 = M;
    L = 50; % "Smear length"

    % Perform LIC
    lic_out = perform_lic(vec_mat, L, options);

    % Use constant-brightness colormap
    cmap = cmocean('phase');
    
    % Use log scale for colors
    mag_mat = log(abs(mag_mat));
    
    % Determine colors over magnitude range only one the first timestep
    if tstep == 1
        mmax = max(mag_mat(:));
        mmin = min(mag_mat(:));
    end
    
    % Determine colors from magnitude
    colors = floor((mag_mat(:)-mmin)/(mmax-mmin).*size(cmap,1));
    
    % Set vectors with the smallest magnitude to the darkest color
    % as opposed to the zeroth color
    colors(colors==0) = 1;
    
    % Multiply b/w img by colors & reshape to n x m x 3
    img = lic_out.*reshape(cmap(colors,:),[size(lic_out) 3]);

    % Plot final colored image
    imshow(img)
    drawnow
    
    % Write image to file
    imwrite(img,sprintf([resultsdir,'/elec_quad_osc_%03d.png'],tstep))
    
    % Increment timestep
    tstep = tstep + 1;
    
    % Only run one timestep
    break
end

