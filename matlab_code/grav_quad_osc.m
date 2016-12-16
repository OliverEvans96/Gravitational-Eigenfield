%% Calculate gravitoelectric field for oscillating gravitational quadrupole

% Determine directory where script is located 
% (not necessarily current matlab directory)
% to determine where to save results
[scriptdir,~,~] = fileparts(mfilename('fullpath'));
resultsdir = [scriptdir,'/../results/'];
videodir = [resultsdir,'movies/'];

% Add required packages to matlab path
% 'genpath' adds subfolders as well
addpath(genpath([scriptdir,'/../required_packages/']))

% Suffixes for positive and negative eigenvalues
labels = {'pos','neg'};

% Video framerate
fps = 30;

imgdir = cell(2);
videoname = cell(2);
videoname_loop = cell(2);
videopath = cell(2);
videopath_loop = cell(2);
video = cell(2); %% NOT CHANGED ON MASTER %%
% Open video files
for pp=1:2
    imgdir{pp} = [resultsdir,sprintf('grav_%s_img/',labels{pp})];
	videoname{pp} = sprintf('grav_%s.avi',labels{pp});
	videoname_loop{pp} = sprintf('grav_%s_loop.avi',labels{pp});
	videopath{pp} = [videodir,videoname{pp}];
	videopath_loop{pp} = [videodir,videoname_loop{pp}];
	video{pp} = VideoWriter(videopath{pp},'Uncompressed AVI');
	video{pp}.FrameRate = fps;
	open(video{pp})
end

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

% Allocate time arrays
% Loop should be over one period & take 4 seconds
nsteps = 4*fps;
tmax = 1;
dt = tmax/nsteps;
tstep = 1;
dt = tmax/nsteps;
for t=0:dt:tmax-dt
	% Eigenvalue matrix
	% Dimensions: ii grid, jj grid, eval # (pos,neg)
	eval_mat = zeros(nn,mm,2);

	% Eigenvector matrix
	% Dimensions: ii grid, jj grid, comp # (x,y), evec # (pos,neg)
	evec_mat = zeros(nn,mm,2,2);

    fprintf('tstep = %d/%d\n',tstep,nsteps)

    % Loop through x,y grid
    for ii=1:nn
        for jj=1:mm
			% Loop through positive/negative evals
				for pp=1:2
				% Shorthand for r and th for this grid point
				rr = r_mesh(ii,jj);
				th = th_mesh(ii,jj);

				% Potential derivatives for GE tensor
				%Err = -12*G*Q./r.^5 .* (3*cos(th).^2 - 1);
				%Ert = -24*G*Q./r.^5 .* sin(th)*cos(th);
				%Ett = 3*G*Q./r.^5 .* (6*cos(th).^2 - sin(th).^2 - 2);

				Err = -4*GG*QQ*kk^2/rr^3 * ((-1 + 3/(kk*rr)^2)*cos(kk*rr-ww*t) ...
					+ 3*sin(kk*rr-ww*t)/(kk*rr)) * (3*cos(th)^2 - 1);
				Ert = -4*GG*QQ*kk^2/rr^3 * ((6/(kk*rr)^2-3)*cos(kk*rr-ww*t) ...
					- (kk*rr-6/(kk*rr))*sin(kk*rr-ww*t)) * sin(th)*cos(th);
				Ett = -GG*QQ*kk^2*(-(-2*kk*rr + 3/(kk*rr))*sin(kk*rr - ww*t) ...
					+ (-kk^2*rr^2 + 3 - 3/(kk^2*rr^2))*cos(kk*rr - ww*t))*sin(th)^2/rr^3;

				% Epp = 2*G*Q*k^2*((-1 + 3/(k^2*r^2))*cos(k*r - w*t) + 3*sin(k*r - w*t)/(k*r))*(3*cos(th)^2 - 1)/r^3;

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
		if tstep == 1
			mmax = max(eval_mat(:));
			mmin = min(eval_mat(:));
		end
		
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

		% Save data to image & video
		writeVideo(video{pp},img);
        imwrite(img,[imgdir{pp},sprintf('grav_%s_%03d.png',labels{pp},tstep)])
    end
    
    tstep = tstep + 1;
end

% Save & close video file
for pp=1:2
    close(video{pp})
end

%%%%%%%%%%%%%%%%%%%
%% Create Videos %%
%%%%%%%%%%%%%%%%%%%

% Write new video file with several loops
disp 'Writing data to movie'

for pp = 1:2
    % Read video data
    % More info: see videoreader.readframe doc page
    % videoData will probably by ~150MB
    videoIn = VideoReader(videopath{pp});
    vidWidth = videoIn.Width;
    vidHeight = videoIn.Height;
    videoData = struct('cdata',zeros(vidHeight,vidWidth,3,'uint8'),...
        'colormap',[]);

    % Ignore last frame for loop because it's the same as first frame
    for kk = 1:nsteps-1
        videoData(kk).cdata = readFrame(videoIn);
    end

    % Open output video file
    videoOut = VideoWriter(videopath_loop{pp});
    videoOut.FrameRate = fps;
    open(videoOut)

    % Write to file, looping several times
    nloops = 20;
    for ll=1:nloops
        % Loop through frames in each loop
        for kk=1:nsteps-1
            writeVideo(videoOut,videoData(kk).cdata)
        end
    end

    % Save and close once again
    close(videoOut)
end

disp 'Finished!'

