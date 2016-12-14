%% Calculate eigenfield for gravitostatic tensor for gravitational quadrupole

% nxm Cartesian grid
nn = 1000;
mm = 1000;
start = 1;
xg = linspace(-start,start,nn);
yg = linspace(-start,start,mm);

% Cartesian mesh
[xx,yy] = meshgrid(xg,yg);

% Polar mesh
[tt,rr] = cart2pol(xx,yy);

% Parameters
G = 1;
Q = 1;

% Eigenvalue matrix
% Dimensions: ii grid, jj grid,
eval_mat = zeros(nn,mm);

% Eigenvector matrix
% Dimensions: ii grid, jj grid, comp # (x,y)
evec_mat = zeros(nn,mm,2);

% Gravitoelectric Tensor
% Err Ert Erp
% Etr Ett Etp
% Epr Ept Epp

% Define components
%Err = -12*G*Q./rr.^5 .* (3*cos(tt).^2 - 1);
%Ert = -24*G*Q./rr.^5 .* sin(tt)*cos(tt);
%Ett = 3*G*Q./rr.^5 .* (6*cos(tt).^2 - sin(tt).^2 - 2);
% Epp = 3*G*Q./rr.^5 .* (5*cos(tt).^2 - 1);

% Loop through x,y grid
for ii=1:nn
    for jj=1:mm
        % Shorthand for r and theta for this grid point
        r1 = rr(ii,jj);
        t1 = tt(ii,jj);
        
        Er = 6*Q/r1^4*(3/2*cos(t1)^2-.5);
        Et = 6*Q/r1^4*cos(t1)*sin(t1);
        
        evecs = [Er;Et];
        
%         % Potential derivatives for GE tensor
%         Err = -12*G*Q./r1.^5 .* (3*cos(t1).^2 - 1);
%         Ert = -24*G*Q./r1.^5 .* sin(t1)*cos(t1);
%         Ett = 3*G*Q./r1.^5 .* (6*cos(t1).^2 - sin(t1).^2 - 2);
%         
%         % Equivalent 2x2 tensor
% %         GE = [Err(ii,jj), Ert(ii,jj); ...
% %               Ert(ii,jj), Ett(ii,jj)];
%         
%         % Equivalent 2x2 tensor
%         GE = [Err, Ert; ...
%               Ert, Ett];
%         
%         % Calculate eigenvalues
%         [evecs,evals] = eig(GE);
%         
%         % Extract eigenvalues from diagonal matrix
%         evals = diag(evals);
%         
        % Convert eigenvectors back to Cartesian
        evecs_cart = zeros([2,1]);
        
        % Vector change of coordinates sph -> cart given by
        % v_cart = (v_r*sin(th)+v_th*cos(th))*x_hat
        %        + (v_r*cos(th)-v_th*sin(th))*y_hat
        % where v_r = evecs(1,:) and v_th = evecs(2,:)
        evecs_cart(1) = evecs(1)*sin(t1) + evecs(2)*cos(t1);
        evecs_cart(2) = evecs(1)*cos(t1) - evecs(2)*sin(t1);
        
        eval = norm(evecs_cart);
        
        % Save "evals" (magnitude) and "evecs" (vector)
<<<<<<< HEAD
        eval_mat(ii,jj) = eval;
=======
        eval_mat(ii,jj) = evals;
>>>>>>> grad_lab
        evec_mat(ii,jj,:) = evecs_cart(:);
        
%         % Save to appropriate matrices
%         eval_mat(ii,jj,:) = evals;
%         % Since direction of eigenvalues is arbitrary, we need to apply
%         % some transformations before plotting to ensure continuity.
%         evec_mat(ii,jj,:,1) = evecs_cart(:,1)...
%             .*sign(evecs_cart(1,1))*sign(xx(ii,jj));
%         evec_mat(ii,jj,:,2) = evecs_cart(:,2)...
%             .*sign(evecs_cart(2,2))*sign(yy(ii,jj));
    end
end

%% Plot using quiver

% figure(1)
% clf
% skip = 5;
% x_sec = xx(1:skip:end,1:skip:end);
% y_sec = yy(1:skip:end,1:skip:end);
% v_sec = vv(1:skip:end,1:skip:end,:);
% 
% quiver(x_sec,y_sec,v_sec(:,:,1),v_sec(:,:,2))
% xlim([-5,5]);
% ylim([-5,5]);


%% LIC

% Plot tiles
titles = {'Static Electric Quadrupole Electric Field'};

figure(1)
clf 
% Store first eigenvector
vv = evec_mat(:,:,:);
ll = eval_mat(:,:);

% Normalize eigenvector array
vv = perform_vf_normalization(vv);

% Create noise canvas
M = randn([nn,mm]);

% parameters for the LIC
options.bound = 'sym';
options.histogram = 'linear'; % keep contrast fixed
options.verb = 1;
options.dt = 1.5; % time steping
% size of the features
options.flow_correction = 1;
%options.isoriented=1;
options.niter_lic = 5; % several iterations gives better results
options.M0 = M;
L = 60; % "Smear length"

% Perform LIC
lic_out = perform_lic(vv, L, options);

% Color output from eigenvalues
cmap = cmocean('phase');
ll = log(abs(ll));
colors = floor((ll(:)-min(ll(:)))/(max(ll(:))-min(ll(:))).*size(cmap,1));
colors(colors==0) = 1;
img = lic_out.*reshape(cmap(colors,:),[size(lic_out) 3]);

% Plot final colored image
%imshow(img)
%title('Static Electric Quadrupole Electric Field')
<<<<<<< HEAD
imwrite(img,'elec_quad_static.png')

=======
imwrite,'elec_quad_static.png')
>>>>>>> grad_lab
