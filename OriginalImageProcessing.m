I = imread("manor.png");
I_gray = rgb2gray(I);

I_new = I_gray;

G_P = cell(1, 7);
L_P = cell(1, 6);
filter = fspecial('gaussian', [11, 11], 1);
G_P{1} = conv2(I_gray, filter, 'same');
filter = fspecial('gaussian', [11, 11], 2);
for level = 2:7
    temp = conv2(G_P{level-1}, filter, 'same');
    G_P{level} = imresize(temp, 0.5);
    L_P{level-1} = G_P{level-1} - imresize(G_P{level}, 2);
end


%%now we are trying to find sift points in level 1, 2, 3 and 4 since we
%%only have 6 layers in total

%%
%find sift keypoints, plot them, and find the neighborhoods(batch) of them
%batch size is 17-by-17
figure(1);
imshow(I);
hold on;
threshold = [25, 25, 25, 25, 25];
Keypoints = cell(1, 5000);
keypoint = [];
contrast = [];
SIFTvectors = [];
orientation = cell(1, 5000);
magnitude = cell(1, 5000);
weighted = cell(1, 5000);
original = cell(1, 5000);
num = 0; %number of keypoints

for level = (1 + 1):(5 + 1)
    Pre_level = imresize(L_P{level-1}, 1/2);
    if (level ~= 6)
        Next_level = imresize(L_P{level+1}, 2);
    end
    [m, n] = size(L_P{level});
    for i = 3:(m - 2)
        for j = 3:(n - 2)
            point = L_P{level}(i, j);
            pre_scale = Pre_level(i-1:i+1, j-1:j+1);
            pre_min = min(min(pre_scale));
            pre_max = max(max(pre_scale));
            if (level ~= 6)
                next_scale = Next_level(i-1:i+1, j-1:j+1);
                next_min = min(min(next_scale));
                next_max = max(max(next_scale));
            end
            cur_scale = L_P{level}(i - 1:i + 1, j - 1:j + 1);
            cur_min = min(min(cur_scale));
            cur_max = max(max(cur_scale));
            iffind = 0;
            if (level ~= 6)
                if (point == cur_max)
                    if (abs(point) > threshold(level-1) && point - next_max > 0 && point - pre_max > 0)
                        iffind = 1;
                        num = num + 1;
                        keypoint = [i, j, 2^(level - 1)];
                        Keypoints{num} = keypoint;
                        if (level == 2)
                            plot(j*2^(level - 1), i*2^(level - 1), "bo", 'MarkerSize', 4);
                        elseif (level == 3)
                            plot(j*2^(level - 1), i*2^(level - 1), "go", 'MarkerSize', 8);
                        elseif (level == 4)
                            plot(j*2^(level - 1), i*2^(level - 1), "yo", 'MarkerSize', 16);
                        else
                            plot(j*2^(level - 1), i*2^(level - 1), "mo", 'MarkerSize', 32);
                        end
                    end
                elseif (point == cur_min)
                    if (pre_min - point > 0 && next_min - point > 0 && abs(point) > threshold(level-1))
                        iffind = 1;
                        num = num + 1;
                        keypoint = [i, j, 2^(level - 1)];
                        Keypoints{num} = keypoint;
                        if (level == 2)
                            plot(j*2^(level - 1), i*2^(level - 1), "bo", 'MarkerSize', 4);
                        elseif (level == 3)
                            plot(j*2^(level - 1), i*2^(level - 1), "go", 'MarkerSize', 8);
                        elseif (level == 4)
                            plot(j*2^(level - 1), i*2^(level - 1), "yo", 'MarkerSize', 16);
                        else
                            plot(j*2^(level - 1), i*2^(level - 1), "mo", 'MarkerSize', 32);
                        end
                    end
                else
                    continue;
                end
            else
                if (point == cur_max)
                    if (point - pre_max > 0 && abs(point) > threshold(level-1))
                        iffind = 1;
                        num = num + 1;
                        keypoint = [i, j, 2^(level - 1)];
                        Keypoints{num} = keypoint;
                        plot(j*2^(level - 1), i*2^(level - 1), "ro", 'MarkerSize', 64);
                    end
                elseif (point == cur_min)
                    if (pre_min - point > 0 && abs(point) > threshold(level-1))
                        iffind = 1;
                        num = num + 1;
                        keypoint = [i, j, 2^(level - 1)];
                        Keypoints{num} = keypoint;
                        plot(j*2^(level - 1), i*2^(level - 1), "ro", 'MarkerSize', 64);
                    end
                else
                    continue;
                end
            end
            if (iffind == 1)
                [o, m, w, or] = Batches(keypoint, G_P{log2(keypoint(3))+1});
                orientation{num} = o;
                magnitude{num} = m;
                weighted{num} = w;
                original{num} = or;
            end
        end
    end
end
hold off
disp('finish finding location of keypoints')
disp('the neighbourhood size is: 3-by-3')
disp('the threshold is:')
disp(threshold)

%now we plot magnitude of gradient, orientation of gradient and weighted
%magnitude of gradient of the neighbourhood(17*17) of one keypoint

%%
%show the batch of a selected keypoint
selectedPointIndex = 20;
figure();
imshow(magnitude{selectedPointIndex}, [min(min(magnitude{selectedPointIndex})), max(max(magnitude{selectedPointIndex}))], 'InitialMagnification', 'fit');
figure();
imshow(original{selectedPointIndex},[0 255] ,'InitialMagnification', 'fit');
hold on
x = 1:17;
y = 1:17;
u = orientation{selectedPointIndex};
v = magnitude{selectedPointIndex};
quiver(x, y, u, v)
hold off
figure()
imshow(weighted{selectedPointIndex}, [min(min(weighted{selectedPointIndex})), max(max(weighted{selectedPointIndex}))], 'InitialMagnification', 'fit');

%%
%calculate the histogram

hist = cell(1, num);
finalKeypointTuple_all = cell(1, num); %39-tuples of keypoints
for n = 1:num
    hist{n} = zeros(2, 36);
    hist{n}(2, 1:36) = 0:35;

    for i = 1:17
        for j = 1:17
            range = orientation{n}(i, j);
            if (range >= 0)
                range = floor(range/10);
                hist{n}(1, range + 1) = hist{n}(1, range + 1) + weighted{n}(i, j);
            else
                range = range + 360;
                range = floor(range/10);
                hist{n}(1, range + 1) = hist{n}(1, range + 1) + weighted{n}(i, j);
            end
        end
    end
    maximum = max(hist{n}(1, 1:36));
    idx = find(hist{n}(1, 1:36) == maximum);
    finalKeypointTuple_all{n} = [Keypoints{n}(1),Keypoints{n}(2), Keypoints{n}(3)];

    for degree = idx:36
        finalKeypointTuple_all{n} = [finalKeypointTuple_all{n}, hist{n}(1, degree)];
    end

    for degree = 1:idx - 1
        finalKeypointTuple_all{n} = [finalKeypointTuple_all{n}, hist{n}(1, degree)];
    end
end

%%
%create the rotated image
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
figure();
imshow((R2), [0 255]);
pause(3);

%%
%now we try to find the sift points in the interest region
%For R1
figure();
R1 = imcrop(R1, [y0 - floor(width / 2), x0 - floor(height / 2), width - 1, height - 1]);
imshow(R1, [0 255]);
hold on
[Key_1, num1] = findSIFT(R1);
hold off
pause(1);
%%
%For R2
figure();
R2 = imcrop(R2, [y0 - floor(width / 2), x0 - floor(height / 2), width - 1, height - 1]);
imshow(R2, [0 255]);
hold on
[Key_2, num2] = findSIFT(R2);
hold off
pause(1);
%%
%For Original
figure();
I_part = imcrop(I_gray, [y0 - floor(width / 2), x0 - floor(height / 2), width - 1, height - 1]);
imshow(I_part, [0 255]);
hold on
[Key_part, num_part] = findSIFT(R2);
hold off

%%
%match sift points in a specific region
%the selected point is (200, 512)

pause(1);
catImage = cat(2, R1, I_part, R2);
figure();
imshow(catImage, [0, 255]);
hold on;

for i = 1:num_part
    midpoint = Key_part{i};
    for s = 1:num1
       
            leftpoint = Key_1{s};
            if (Bhattacharya_coefficient(midpoint, leftpoint, 0.9) == 1)
                plot(midpoint(2)*midpoint(3) + 256, midpoint(1)*midpoint(3), 'yo');
                plot(leftpoint(2)*leftpoint(3), leftpoint(1)*leftpoint(3), 'mo');
                plot([midpoint(2)*midpoint(3) + 256, leftpoint(2)*leftpoint(3)], [midpoint(1)*midpoint(3), leftpoint(1)*leftpoint(3)], 'g');
                
                break;
            end
        
    end
    for t = 1:num2
        
            rightpoint = Key_2{t};
            if (Bhattacharya_coefficient(midpoint, rightpoint, 0.9) == 1)
                plot(midpoint(2)*midpoint(3) + 256, midpoint(1)*midpoint(3), 'yo');
                plot(rightpoint(2)*rightpoint(3) + 512, rightpoint(1)*rightpoint(3), 'ro');
                plot([midpoint(2)*midpoint(3) + 256, rightpoint(2)*rightpoint(3) + 512], [midpoint(1)*midpoint(3), rightpoint(1)*rightpoint(3)], 'g');
                
                break;
            end
       
    end
    
end
hold off





    


    
