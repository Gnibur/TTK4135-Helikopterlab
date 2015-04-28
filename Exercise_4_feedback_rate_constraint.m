%% Exercise 4: LQ controller

Exercise4_rate_constraint;

Q = diag([1 0 0 0 0 0]);
R = diag([1 1]);

[K,S,e] = dlqr(A,B,Q,R);