% imshow test
x = 0:0.01:1;
y = 0:0.01:1;

[X,Y] = meshgrid(x,y);

Z = Y;
C = X;

% Colormap
A = cmocean('ice');
F = floor(C(:).*size(A,1));
F(F==0) = 1;
img = Z.*reshape(A(F,:),[size(Z) 3]);
figure(1)
imshow(img)