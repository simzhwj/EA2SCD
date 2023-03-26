function [H,iternum,obj] = SymmetricNMF(V,H_init,epsilon,itermax)
% SymmetricNMF - The Non-negative Matrix Factorization algorithm.
%   To do the NMF for matrix V. The formula is shown as follows:
%   V(n*m) = H(n*k)*H'(k*n)
%   For example, (n,n,k) = (1000,1000,8) or (112,112,2).
% 
%   Here are some message about NMF with the URL as follows:
%   https://zhuanlan.zhihu.com/p/27460660
%   https://blog.csdn.net/pipisorry/article/details/52098864
%
%   [H] = StandardNMF(V,H_init,epsilon,iternum)
% 
%   Input - 
%   V: the n*m data matrix, m n-vectors arranging by columns;
%   H_init: initialize Coefficient matrix or Dictionary matrix
%   epsilon: error of distance between V and V' of two iterations;
%   itermax: iteration number upper limit.
%
%   Output - 
%   H : the k*n coefficient matrix, each columns are obtained by 
%       projecting each column of V matrix onto W matrix;
%   iternum : the number of interations consumed to meet epsilon;
%   obj : the objective between V and H_new*H_new' in each iter
%   
%   Copyright (c) 2018 CHEN Tianyang
%   more info contact: tychen@whu.edu.cn

%% pre-work
% n is the number of dimension, n is the number of n-vectors, k is the number of basis vectors;
[n] = size(V,1);
k = size(H_init,2);
% initialize H
H_old = H_init;
H_new = zeros(n,k);
% some prepare variables
obj = zeros(itermax,1); % The value of objective function
count = 1; % The number of itertions
obj(count) = sum(sum( (V-H_init*H_init').^2) );
error = realmax; % The difference between obj(count-1) and obj(count);

%% iterate
% iterate as the following formula:
% 1. H_new(i,j) = 0.5*H_old(i,j)*[1+(V*H_old)(i,j)/(H_old*H_old'*H_old)(i,j)]
while error >= epsilon
    count = count+1;
    % update matrix H
    Hcoematx_up = V*H_old + H_old*(H_old')*H_old;
    Hcoematx_dn = H_old*(H_old')*H_old;
    for i=1:n
        for j=1:k
            if Hcoematx_dn(i,j)==0
                H_new(i,j) = H_old(i,j);
            else
                H_new(i,j) = 0.5*H_old(i,j)*Hcoematx_up(i,j)/Hcoematx_dn(i,j);
            end
        end
    end
    % calculate the difference between two iteration approximation matrices 
    obj(count) = sum(sum( (V-H_new*H_new').^2) );
    error = obj(count-1)-obj(count);

    if mod(count,100)==0
        %fprintf('SymmetricNMF：%d轮迭代误差为%d.\n',count,error);
    end
    % The results of this round of iteration
    if count == itermax
        %fprintf('%d轮迭代已毕，误差仍未收敛.\n',itermax);
        break;
    end
    % prepare for the next round of round
    H_old = H_new;
end

%% get result
H = H_new;
iternum = count;
end
%%