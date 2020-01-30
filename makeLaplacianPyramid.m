I = imread("manor.png");
I_gray = rgb2gray(I);

I_new = I_gray;

G_P = cell(1, 7);
L_P = cell(1, 6);
filter = fspecial('gaussian', [11, 11], 1);
G_P{1} = conv2(I_gray, filter, 'same');
filter = fspecial('gaussian', [11, 11], 2);
figure();
imshow(G_P{1}, [0 255]);
for level = 2:7
    temp = conv2(G_P{level - 1}, filter, 'same');
    G_P{level} = imresize(temp, 0.5);
    L_P{level - 1} = G_P{level - 1} - imresize(G_P{level}, 2);
    figure();
    imshow(L_P{level - 1}, [min(min(L_P{level - 1})), max(max(L_P{level - 1}))]);
end