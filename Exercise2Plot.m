clc;

load('travel.mat');
figure(1);
%plot(ans(1, :), ans(2, :) - (ans(2,200) - 180) * ans(2, :) ./ ans(2, :));
plot(ans(1, :), ans(2, :));

load('pitch.mat');
figure(2);
%plot(ans(1, :), ans(2, :) - (ans(2,200) - 180) * ans(2, :) ./ ans(2, :));
plot(ans(1, :), ans(2, :));