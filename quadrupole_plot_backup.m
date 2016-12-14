% Plot gravitational quadrupole using quiver
% n=201;
% m=401;
% tx = linspace(-5,5,n);
% ty = linspace(-5,5,m);

% nxm grid
nn = 500;
mm = 500;
xg = linspace(-5,5,nn);
yg = linspace(-5,5,mm);

% Cartesian mesh
[xx,yy] = meshgrid(xg,yg);

% Polar mesh
[tt,rr] = cart2pol(xx,yy);

% Parameters
G = 1;
Q = 1;

% Polar vector field
%Err = 12*G*Q./rr.^5 .* (3 .* cos(tt).^2 - 1);
%Ert = -24*G*Q./rr.^5 .* (sin(tt) .* cos(tt));

% Eigenvalue matrix
% Dimensions: ii grid, jj grid, eval # (pos,neg)
eval_mat = zeros(nn,mm,2);

% Eigenvector matrix
% Dimensions: ii grid, jj grid, comp # (x,y), evec # (pos,neg)
evec_mat = zeros(nn,mm,2,2);

% Gravitoelectric Tensor
% Err Ert Erp
% Etr Ett Etp
% Epr Ept Epp

% Define components
Err = -12*G*Q./rr.^5 .* (3*cos(tt).^2 - 1);
Ert = -24*G*Q./rr.^5 .* sin(tt)*cos(tt);
Ett = 3*G*Q./rr.^5 .* (6*cos(tt).^2 - sin(tt).^2 - 2);
Epp = 3*G*Q./rr.^5 .* (5*cos(tt).^2 - 1);

% Loop through x,y grid
for ii=1:nn
    for jj=1:mm
        % Define tensor
        % Full tensor
%         GE = [Err(ii,jj), Ert(ii,jj),          0; ...
%               Ert(ii,jj), Ett(ii,jj),          0; ...
%                        0,          0, Epp(ii,jj)];
        
        % Equivalent 2x2 tensor
        GE = [Err(ii,jj), Ert(ii,jj); ...
              Ert(ii,jj), Ett(ii,jj)];
        
        % Calculate eigenvalues
        [evecs,evals] = eig(GE);
        
        % Extract eigenvalues from diagonal matrix
        evals = diag(evals);
        
        % Save to appropriate matrices
        eval_mat(ii,jj,:) = evals;
        evec_mat(ii,jj,:,:) = evecs;
    end
end

% Convert to cartesian vector field
%Fx = Err.*cos(Ert);
%Fy = Ert.*sin(Ert);

% Normalize all vectors
% N = sqrt(Fx.^2 + Fy.^2);
% Prevent NANs
% for ii=1:nn
%     for jj=1:mm
%         if N(ii,jj) > 1
%            N(ii,jj) = 0;
%         end
%     end
% end

%Fx = Fx/N;
%Fy = Fy/N;

% vv = zeros(nn,mm,2);
% % Create nxnx2 array
% for kk = 1:nn
% for jj = 1:nn
%     vv(kk,jj,1) = Fx(kk,jj);
%     vv(kk,jj,2) = Fy(kk,jj);
% end
% end

vv = evec_mat(:,:,:,1);

%figure(2)

% Normalize
vv = perform_vf_normalization(vv);

%% Plot using quiver
figure(1)
clf
skip = 50;
x_sec = xx(1:skip:end,1:skip:end);
y_sec = yy(1:skip:end,1:skip:end);
v_sec = vv(1:skip:end,1:skip:end,:);

size(x_sec)
size(y_sec)
size(v_sec)
quiver(x_sec,y_sec,v_sec(:,:,1),v_sec(:,:,2))
xlim([-5,5]);
ylim([-5,5]);


%% LIC
figure(2)
clf
M = randn([nn,mm]);
size(vv)
size(M)
% parameters for the LIC
options.histogram = 'linear'; % keep contrast fixed
options.verb = 0;
options.dt = 1.5; % time steping
% size of the features
options.flow_correction = 1;
options.niter_lic = 2; % several iterations gives better results
options.M0 = M;
L = 10; % "Smear length"
lic_out = perform_lic(vv, L, options);
imshow(lic_out)