%% Exercise 3: LQ controller

Exercise2;

Q = diag([1 0 0 0]);
R = .5;

[K,S,e] = dlqr(A,B,Q,R);