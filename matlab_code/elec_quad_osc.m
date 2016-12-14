%% Calculate electric field for oscillating electric dipole

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
Q = 1e-3;
w = 2*pi;
k = 1e1;

% Eigenvalue matrix
% Dimensions: ii grid, jj grid, eval # (pos,neg)
eval_mat = zeros(nn,mm);

% Eigenvector matrix
% Dimensions: ii grid, jj grid, comp # (x,y), evec # (pos,neg)
evec_mat = zeros(nn,mm,2);

% Create noise canvas
M = randn([nn,mm]);

% Loop through x,y grid
nsteps = 50;
tmax = 1;
dt = tmax/nsteps;
tstep = 1;
for t=0:dt:tmax
    % Eigenvalue matrix
    % Dimensions: ii grid, jj grid, eval # (pos,neg)
    eval_mat = zeros(nn,mm);
    
    % Eigenvector matrix
    % Dimensions: ii grid, jj grid, comp # (x,y), evec # (pos,neg)
    evec_mat = zeros(nn,mm,2);
    fprintf('\nt=%f\n',t)
    for ii=1:nn
        %fprintf('ii=%d\n',ii)
        for jj=1:mm
            % Shorthand for r and theta for this grid point
            r = rr(ii,jj);
            th = tt(ii,jj);

            % Potential derivatives for GE tensor
            %Err = -12*G*Q./r.^5 .* (3*cos(t).^2 - 1);
            %Ert = -24*G*Q./r.^5 .* sin(t)*cos(t);
            %Ett = 3*G*Q./r.^5 .* (6*cos(t).^2 - sin(t).^2 - 2);

            Er = -2*Q*k^2/r^2 * (cos(k*r - w*t)*(1-3/(k*r)^2)... 
                - 3/(k*r)*sin(k*r-w*t)) * (3/2*cos(th)^2-.5);
            Eth = -Q*k^2/r^2 * (sin(k*r - w*t)*(k*r - 6/(k*r)) ...
                + (3-6/(k*r)^2)*cos(k*r-w*t))*cos(th)*sin(th);
            
            if ii==190 && jj == 73
                fprintf('Eth = %f',Eth)
            end

            evecs = [Er;Eth];

            evecs_cart(1) = evecs(1)*sin(th) + evecs(2)*cos(th);
            evecs_cart(2) = evecs(1)*cos(th) - evecs(2)*sin(th);
            
            evals = sqrt(evecs_cart(1)^2+evecs_cart(2)^2);

            if ii==10  && jj == 89
                fprintf('\nEvec_cart = %f,%f\n',evecs_cart)
                fprintf('\nunit = %f,%f\n',evecs_cart/evals);
            end
            
            % Save to appropriate matrices
            eval_mat(ii,jj) = evals;
            evec_mat(ii,jj,:) = evecs_cart(:);

            % Since direction of eigenvalues is arbitrary, we need to apply
            % some transformations before plotting to ensure continuity.
            %evec_mat(ii,jj,:,1) = evecs_cart(:,1)...
            %    .*sign(evecs_cart(1,1))*sign(xx(ii,jj));
            %evec_mat(ii,jj,:,2) = evecs_cart(:,2)...
            %    .*sign(evecs_cart(2,2))*sign(yy(ii,jj));
            %evecs_cart(:,:)
            %sign(evecs_cart(1,:))
            %evecs_cart.*sign(evecs_cart(:,1))
            %break
        end
        %break
    end


    %%% LIC

    % Plot tiles
    figure(1)
    clf 
    % Store first eigenvector
    vv = evec_mat(:,:,:);
    ll = eval_mat(:,:);

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
    options.niter_lic = 3; % several iterations gives better results
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
    drawnow
    title('Electric Field of Oscillating Electric Dipole')
    imwrite(img,sprintf('img/elec_quad_osc_%03d.png',tstep))
    
    tstep = tstep + 1;
end

