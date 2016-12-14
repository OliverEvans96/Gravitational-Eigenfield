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
G = 1e-2;
Q = 1e-2;
k = 1e1;
w = 2*pi;

% Create noise canvas
M = randn([nn,mm]);

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
%Err = -12*G*Q./rr.^5 .* (3*cos(tt).^2 - 1);
%Ert = -24*G*Q./rr.^5 .* sin(tt)*cos(tt);
%Ett = 3*G*Q./rr.^5 .* (6*cos(tt).^2 - sin(tt).^2 - 2);
% Epp = 3*G*Q./rr.^5 .* (5*cos(tt).^2 - 1);

% Loop through time
tstep = 1;
nframes = 50;
tmax=1;
dt = tmax/nframes;
for t=0:dt:tmax
    % Loop through x,y grid
    for ii=1:nn
        for jj=1:mm
            % Shorthand for r and th for this grid point
            r = rr(ii,jj);
            th = tt(ii,jj);

            % Potential derivatives for GE tensor
            %Err = -12*G*Q./r.^5 .* (3*cos(th).^2 - 1);
            %Ert = -24*G*Q./r.^5 .* sin(th)*cos(th);
            %Ett = 3*G*Q./r.^5 .* (6*cos(th).^2 - sin(th).^2 - 2);

            Err = -4*G*Q*k^2/r^3 * ((-1 + 3/(k*r)^2)*cos(k*r-w*t) ...
                + 3*sin(k*r-w*t)/(k*r)) * (3*cos(th)^2 - 1);
            Ert = -4*G*Q*k^2/r^3 * ((6/(k*r)^2-3)*cos(k*r-w*t) ...
                - (k*r-6/(k*r))*sin(k*r-w*t)) * sin(th)*cos(th);
            Ett = -G*Q*k^2*(-(-2*k*r + 3/(k*r))*sin(k*r - w*t) ...
                + (-k^2*r^2 + 3 - 3/(k^2*r^2))*cos(k*r - w*t))*sin(th)^2/r^3;

            % Epp = 2*G*Q*k^2*((-1 + 3/(k^2*r^2))*cos(k*r - w*t) + 3*sin(k*r - w*t)/(k*r))*(3*cos(th)^2 - 1)/r^3;

            % Equivalent 2x2 tensor
    %         GE = [Err(ii,jj), Ert(ii,jj); ...
    %               Ert(ii,jj), Ett(ii,jj)];

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
                .*sign(evecs_cart(1,1))*sign(xx(ii,jj));
            evec_mat(ii,jj,:,2) = evecs_cart(:,2)...
                .*sign(evecs_cart(2,2))*sign(yy(ii,jj));
            %evecs_cart(:,:)
            %sign(evecs_cart(1,:))
            %evecs_cart.*sign(evecs_cart(:,1))
            %break
        end
        %break
    end
    
    %%% LIC

    % Plot tiles
    titles = {'Gravitational Quadrupole Tendex Field: Positive Eigenvalue',...
        'Gravitational Quadrupole Tendex Field: Negative Eigenvalue'};
    labels = {'pos','neg'};

    for kk=1:2
        figure(kk)
        clf 
        % Store first eigenvector
        vv = evec_mat(:,:,:,kk);
        ll = eval_mat(:,:,kk);

        % Normalize eigenvector array
        vv = perform_vf_normalization(vv);

        % parameters for the LIC
        options.bound = 'sym';
        options.histogram = 'linear'; % keep contrast fixed
        options.verb = 1;
        options.dt = 1.5; % time steping
        % size of the features
        options.flow_correction = 1;
        %options.isoriented=1;
        options.niter_lic = 2; % several iterations gives better results
        options.M0 = M;
        L = 30; % "Smear length"

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
        %drawnow
        %title(titles{kk})
        imwrite(img,sprintf('img/grav_quad_osc_%s_%03d.png',labels{kk},tstep))
    end
    
    tstep = tstep + 1;
end

