function [I_rotated] = imageRotation(I, x0, y0, theta, scale)
[m, n, ~] = size(I);
I = padarray(I, [m, n], 0);
I_copy = zeros(3*m, 3*n);
x0new = x0 + m;
y0new = y0 + n;
[m_new, n_new] = size(I);
xCenter = int32(floor(m_new / 2 + 1));
yCenter = int32(floor(n_new / 2 + 1));     %central coordiante of the padded image

x_shift = xCenter - x0new;
y_shift = yCenter - y0new;     %shift value

for i = 1 + m + x_shift : 2*m + x_shift
    for j = 1 + n + y_shift : 2*n + y_shift        
        I_copy(i, j) = I(i - x_shift, j - y_shift);        
    end
end
I_copy = imrotate(I_copy, theta, 'bicubic', 'crop');           %rotate the image
I_copy = imresize(I_copy, scale);
I_copy = padarray(I_copy, [m, n], 0);
[M, N, ~] = size(I_copy);
xCenter_new = int32(floor(M / 2 + 1));
yCenter_new = int32(floor(N / 2 + 1));
I_copy = imcrop(I_copy, [yCenter_new - (y0-1),xCenter_new - (x0-1), n-1, m-1]);

I_rotated = I_copy;


        



