% imshow test
x = 0:0.01:1;
y = 0:0.01:1;

[X,Y] = meshgrid(x,y);

Z = Y;

figure(1)
imshow(Z)