%% Calculate gravitoelectric field for static gravitational quadrupole

% Determine directory where script is located 
% (not necessarily current matlab directory)
% to determine where to save results
[scriptdir,~,~] = fileparts(mfilename('fullpath'));
resultsdir = [scriptdir,'/../results/'];

% Add required packages to matlab path
% 'genpath' adds subfolders as well
addpath(genpath([scriptdir,'/../required_packages/']))

% Suffixes for positive and negative eigenvalues
labels = {'pos','neg'};

% Determine image directory
imgdir = [resultsdir,sprintf('static_img/')];

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
GG = 1e-2;
QQ = 1e-2;
kk = 1e1;
ww = 2*pi;

% Create noise canvas for LIC
MM = randn([nn,mm]);

% Eigenvalue matrix
% Dimensions: ii grid, jj grid, eval # (pos,neg)
eval_mat = zeros(nn,mm,2);

% Eigenvector matrix
% Dimensions: ii grid, jj grid, comp # (x,y), evec # (pos,neg)
evec_mat = zeros(nn,mm,2,2);

% Loop through x,y grid
for ii=1:nn
    for jj=1:mm
        % Loop through positive/negative evals
        for pp=1:2
            % Shorthand for r and th for this grid point
            rr = r_mesh(ii,jj);
            th = th_mesh(ii,jj);

            % Potential derivatives for GE tensor
            Err = -12*GG*QQ./rr.^5 .* (3*cos(th).^2 - 1);
            Ert = -24*GG*QQ./rr.^5 .* sin(th)*cos(th);
            Ett = 3*GG*QQ./rr.^5 .* (6*cos(th).^2 - sin(th).^2 - 2);
            
            % Equivalent 2x2 tensor
            GE = [Err, Ert; ...
                  Ert, Ett];

            % Calculate eigenvalues
            [evecs,evals] = eig(GE);

            % Extract eigenvalues from diagonal matrix
            evals = diag(evals);

            % Convert eigenvectors back to Cartesian
            evecs_cart = zeros([2,2]);

            % Vector change of coordinates sph -> cart given by
            % v_cart = (v_r*sin(th)+v_th*cos(th))*x_hat
            %        + (v_r*cos(th)-v_th*sin(th))*y_hat
            % where v_r = evecs(1,:) and v_th = evecs(2,:)
            evecs_cart(1,:) = evecs(1,:)*sin(th) + evecs(2,:)*cos(th);
            evecs_cart(2,:) = evecs(1,:)*cos(th) - evecs(2,:)*sin(th);

            % Save to appropriate matrices
            eval_mat(ii,jj,:) = evals;
            % Since direction of eigenvalues is arbitrary, we need to apply
            % some transformations before plotting to ensure continuity.
            evec_mat(ii,jj,:,1) = evecs_cart(:,1)...
                .*sign(evecs_cart(1,1))*sign(x_mesh(ii,jj));
            evec_mat(ii,jj,:,2) = evecs_cart(:,2)...
                .*sign(evecs_cart(2,2))*sign(y_mesh(ii,jj));
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Line Integral Convolution %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Normalize vector field
evec_mat = perform_vf_normalization(evec_mat);

% LIC Parameters
options.bound = 'sym';
options.histogram = 'linear'; % keep contrast fixed
options.verb = 1;
options.dt = 1.5; % time steping
options.flow_correction = 3;
options.isoriented=0;
options.niter_lic = 3; % several iterations gives better results
options.M0 = MM;
LL = 50; % "Smear length"

% Use constant-brightness colormap
cmap = cmocean('phase');

for pp = 1:2
    % Perform LIC
    lic_out = perform_lic(evec_mat(:,:,:,pp), LL, options);
    
    % Use log scale for colors
    eval_mat = log(abs(eval_mat));

    % Determine colors over magnitude range only on the first timestep
    mmax = max(eval_mat(:));
    mmin = min(eval_mat(:));
    
    % Handle out of bounds colors by clipping them
    eval_mat(eval_mat<mmin) = mmin;
    eval_mat(eval_mat>mmax) = mmax;
    
    % Determine colors from magnitude
    colors = floor((eval_mat(:,:,pp) - mmin)...
        / (mmax-mmin).*size(cmap,1));
    
    % Set vectors with the smallest magnitude to the first color
    % (as opposed to the zeroth color)
    colors(colors==0) = 1;
    
    % Colors to column vector
    colors_col_vec = reshape(colors,numel(colors),1);
    
    % Multiply b/w img by colors & reshape to n x m x 3
    img = lic_out.*reshape(cmap(colors_col_vec,:),[size(lic_out) 3]);

    imwrite(img,[imgdir,sprintf('grav_%s_static.png',labels{pp})])
end

disp 'Finished!'

