%% Calculate electric field for oscillating electric dipole

% Determine directory where script is located 
% (not necessarily current matlab directory)
% to determine where to save results
[scriptdir,~,~] = fileparts(mfilename('fullpath'));
resultsdir = [scriptdir,'/../results/'];
imgdir = [resultsdir,'elec_img/'];
moviedir = [resultsdir,'movies/'];

% Add required packages to matlab path
% 'genpath' adds subfolders as well
addpath(genpath([scriptdir,'/../required_packages/']))


% Create video file to write
videoname = 'elec_osc.avi';
videoname_loop = 'elec_osc_loop.avi';
videopath = [moviedir,videoname];
videopath_loop = [moviedir,videoname_loop];
video = VideoWriter(videopath,'Uncompressed AVI');
fps = 30;
video.FrameRate = fps;
open(video)

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
ww = 2*pi;
kk = 1e1;

% Create noise canvas for LIC
MM = randn([nn,mm]);

% Allocate time arrays
% Loop should be over one period & take 4 seconds
nsteps = 4*fps;
tmax = 1;
dt = tmax/nsteps;
tstep = 1;

% Loop through time
for tt=0:dt:tmax
    % Magnitude matrix
    % Dimensions: ii grid, jj grid
    mag_mat = zeros(nn,mm);
    
    % Vector matrix
    % Dimensions: ii grid, jj grid, comp # (x,y)
    vec_mat = zeros(nn,mm,2);
    
    % Print current timestep
    fprintf('tstep = %d/%d\n',tstep,nsteps)
    
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
            Er = -2*QQ*kk^2/rr^2 * (cos(kk*rr - ww*tt)*(1-3/(kk*rr)^2)... 
                - 3/(kk*rr)*sin(kk*rr-ww*tt)) * (3/2*cos(th)^2-.5);
            Eth = -QQ*kk^2/rr^2 * (sin(kk*rr - ww*tt)*(kk*rr - 6/(kk*rr)) ...
                + (3-6/(kk*rr)^2)*cos(kk*rr-ww*tt))*cos(th)*sin(th);

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
    
    % Determine colors over magnitude range only on the first timestep
    if tstep == 1
        mmax = max(mag_mat(:));
        mmin = min(mag_mat(:));
    end
    
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
    
    % Write data image & video
    imwrite(img,sprintf([imgdir,'elec_quad_osc_%03d.png'],tstep));
    writeVideo(video,img);
    
    % Increment timestep
    tstep = tstep + 1;
    
end

% Save & close video file
close(video)

%%%%%%%%%%%%%%%%%%%
%% Create Videos %%
%%%%%%%%%%%%%%%%%%%

% Write new video file with several loops
disp 'Writing data to movie'

% Read data
% More info: see videoreader.readframe doc page
% videoData will probably by ~150MB
videoIn = VideoReader(videopath);
vidWidth = videoIn.Width;
vidHeight = videoIn.Height;
videoData = struct('cdata',zeros(vidHeight,vidWidth,3,'uint8'),...
    'colormap',[]);

% Ignore last frame for loop because it's the same as first frame
for kk = 1:nsteps-1
    videoData(kk).cdata = readFrame(videoIn);
end

% Open output video file
videoOut = VideoWriter(videopath_loop);
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

disp 'Finished!'

