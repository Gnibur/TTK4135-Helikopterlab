clc;

load('travel.mat');
plot(ans(1, :), ans(2, :) - (ans(2,200) - 180) * ans(2, :) ./ ans(2, :));
plot(ans(1, :), ans(2, :));