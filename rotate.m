I = imread("manor.png");
I_gray = rgb2gray(I);

x0 = 512;
y0 = 512;
scale1 = 1.1;
scale2 = 0.9;
theta1 = -10;
theta2 = -350;
width = 256;
height = 256;       %size of interst region

R1 = imageRotation(I_gray, x0, y0, theta1, scale1);
R2 = imageRotation(I_gray, x0, y0, theta2, scale2);

figure();
imshow(R1, [0 255]);
hold on;
plot(y0, x0, 'ro');
hold off;
pause(1);

figure();
imshow(R2, [0 255]);
hold on;
plot(y0, x0, 'ro');
hold off;


