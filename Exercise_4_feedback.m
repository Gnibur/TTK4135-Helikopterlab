%% Exercise 3: LQ controller

Exercise4;

Q = diag([1 0 0 0 0 0]);
R = diag([1 1]);

[K,S,e] = dlqr(A,B,Q,R);