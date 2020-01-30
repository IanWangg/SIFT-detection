I = imread("manor.png");
I_gray = rgb2gray(I);

I_new = I_gray;

G_P = cell(1, 7);
G_P{1} = imgaussfilt(I_new, 1);
figure();
imshow(G_P{1});
for level = 2:7
    temp = imgaussfilt(G_P{level - 1}, 2);
    G_P{level} = imresize(temp, 0.5);
    figure();
    imshow(G_P{level});
end

