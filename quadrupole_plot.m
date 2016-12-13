% Plot gravitational quadrupole using quiver
% n=201;
% m=401;
% tx = linspace(-5,5,n);
% ty = linspace(-5,5,m);

% nxm grid
nn = 100;
mm = 100;
xg = linspace(-5,5,nn);
yg = linspace(-5,5,mm);

% Cartesian mesh
[xx,yy] = meshgrid(xg,yg);

size(xx)

% Polar mesh
[rr,tt] = cart2pol(xx,yy);

size(rr)

% Parameters
G = 1;
Q = 1;

% Polar vector field
Err = 12*G*Q./rr.^5 .* (3 .* cos(tt).^2 - 1);
Ert = -24*G*Q./rr.^5 .* (sin(tt) .* cos(tt));

% Convert to cartesian vector field
Fx = Err.*cos(Ert);
Fy = Ert.*sin(Ert);

size(Fx)
% Normalize all vectors
N = sqrt(Fx.^2 + Fy.^2);
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

vv = zeros(nn,mm,2);
% Create nxnx2 array
for kk = 1:nn
for jj = 1:nn
    vv(kk,jj,1) = Fx(kk,jj);
    vv(kk,jj,2) = Fy(kk,jj);
end
end

%figure(2)

% Normalize
vv = perform_vf_normalization(vv);

fprintf('test1\n')

%Plot using quiver
figure(1)
skip = 1;
x_sec = x(1:skip:end,1:skip:end);
y_sec = y(1:skip:end,1:skip:end);
v_sec = vv(1:skip:end,1:skip:end,:);
quiver(x_sec,y_sec,v_sec(:,:,1),v_sec(:,:,2))
xlim([-5,5]);
ylim([-5,5]);

fprintf('test\n')

%% LIC
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