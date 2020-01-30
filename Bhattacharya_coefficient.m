function [isPair] = Bhattacharya_coefficient(point1, point2, threshold)
B_co = 0;
sum1 = 0;
sum2 = 0;
for i = 4:39
    sum1 = sum1 + point1(i);
    sum2 = sum2 + point2(i);
end

point1 = point1 / sum1;
point2 = point2 / sum2;


for i = 4:39
    B_co = B_co + (point1(i)*point2(i))^(1/2);
end

if (B_co >= threshold)
    isPair = 1;
else
    isPair = 0;
end
        