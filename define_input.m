function [A,lhs,rhs,lb,ub,f,cluster_id_from,cluster_id_to] = define_input()
%Define the problem through A, lhs, rhs, lb, ub, f

A = [4    -1  0.9    0    0    0    0    0    0 ; 
     0     1   -1    0    0    0    0    0    0 ;    
     0   0.9   -1    4    0    0    0    0    0 ;
     0     0    0    5    6    0    0    0    0 ;
     0     0    0    0    7   -1  0.9    0    0 ; 
     0     0    0    0    0   -1    1    0    0 ; 
     0     0    0    0    0  0.9   -1    7    5 ; 
     0     0    0    0    0    0    0    1    9];
 
A=sparse(A);

lhs = [ 0;  0; 0;  0;   0;  0;  0;  0];
rhs = [13; 0; 15; 56;  29;  0; 30; 39];
lb = [0 0 0 0 0 0 0 0 0];
ub = [6 10 10 6 2 8 8 2 7];

f=[-1 -1.5 -2.2 -1.3 -1.4 -1.5 -1.6 -1.7 -1.8];

%annotation für Clustering
cluster_id_from = [10  10  100  100 100   100 1000 1000 1000];
cluster_id_to   = [ 0 100   10    0   0  1000  100    0    0];


end

