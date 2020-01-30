function [finalKeypointTuple_all, num] = findSIFT(I)
I_gray = I;
G_P = cell(1, 7);
L_P = cell(1, 6);
filter = fspecial('gaussian', [11, 11], 1);
G_P{1} = conv2(I_gray, filter, 'same');
filter = fspecial('gaussian', [11, 11], 2);
for level = 2:7
    temp = conv2(G_P{level - 1}, filter, 'same');
    G_P{level} = imresize(temp, 0.5);
    L_P{level - 1} = G_P{level - 1} - imresize(G_P{level}, 2);
end



%%now we are trying to find sift points in level 1, 2, 3 and 4 since we
%%only have 6 layers in total

%%
%find sift keypoints, plot them, and find the neighborhoods(batch) of them
%batch size is 17-by-17


threshold = [25 25 25 25 25];
Keypoints = cell(1, 5000);
keypoint = [];
contrast = [];
SIFTvectors = [];
orientation = cell(1, 5000);
magnitude = cell(1, 5000);
weighted = cell(1, 5000);
original = cell(1, 5000);
num = 0; %number of keypoints
level = 2;
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
                    if (pre_min - point > 0 &&  abs(point) > threshold(level-1))
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
disp('finish finding location of keypoints')
disp('the neighbourhood size is: 3-by-3')
disp('the threshold is:')
disp(threshold)

hist = cell(1, num);
finalKeypointTuple_all = cell(1, num);      %39-tuples of keypoints
for n = 1:num
    hist{n} = zeros(2, 36);
    hist{n}(2,1:36) = 0:35;
    
    for i = 1:17
        for j = 1:17
            range = orientation{n}(i, j);
            if (range >= 0) 
                range = floor(range/10);
                hist{n}(1, range+1) = hist{n}(1, range+1) + weighted{n}(i, j);
            else
                range = range + 360;
                range = floor(range/10);
                hist{n}(1, range+1) = hist{n}(1, range+1) + weighted{n}(i, j);
            end
        end
    end
    maximum = max(hist{n}(1,1:36));
    idx = find(hist{n}(1,1:36)==maximum);
    finalKeypointTuple_all{n} = [Keypoints{n}(1:2) Keypoints{n}(3)];
    
    for degree = idx:36
        finalKeypointTuple_all{n} = [finalKeypointTuple_all{n} hist{n}(1, degree)];
    end
    
    for degree = 1:idx-1
        finalKeypointTuple_all{n} = [finalKeypointTuple_all{n} hist{n}(1, degree)];
    end
end