n = 201;    % size of the image
% sigma = 60; % regularity of the vector field
% options.bound = 'sym'; % boundary handling
% v = perform_blurring(randn(n,n,2), sigma, options);
% v = perform_vf_normalization(v);

nn  = 201;
xx = linspace(0,10,nn);
vv = zeros(nn,nn,2);
v = zeros(size(vv));


[x,y] = meshgrid(xx,xx);
Fx = x;
Fy = y;
for kk = 1:nn
for jj = 1:nn
    vv(kk,jj,1) = Fx(kk,jj);
    vv(kk,jj,2) = Fy(kk,jj);
end
end
%v = vv/norm(vv);
v = perform_vf_normalization(vv);


M = randn(n);
% parameters for the LIC
options.histogram = 'linear'; % keep contrast fixed
options.verb = 0;
options.dt = 1.5; % time steping
% size of the features
options.flow_correction = 1;
options.niter_lic = 2; % several iterations gives better results
% iterated lic
Mlist = {};
%for i=1:4
%    options.M0 = M;
%    Mlist{end+1} = perform_lic(v, i, options);
%end
out_image = perform_lic(v,5,options);
out_cut = out_image(1:20,1:20);
out_rot = imrotate(out_image,90);

% display
figure(1)
clf; %imageplot(Mlist,'',2,2);
imshow(out_rot)
axis equal


figure(2); clf;
% quiver(x,y,vv(:,:,1),vv(:,:,2))
skip = 5;
x_sec = x(1:skip:end,1:skip:end);
y_sec = y(1:skip:end,1:skip:end);
v_sec = vv(1:skip:end,1:skip:end,:);
quiver(x_sec,y_sec,v_sec(:,:,1),v_sec(:,:,2))
axis equal
