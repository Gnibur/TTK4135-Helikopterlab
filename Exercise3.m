%% Exercise 3: LQ controller

hints_problem2;


Q = diag([1 0 0 0]);
R = 0.1;



[K,S,e] = dlqr(A,B,Q,R);





