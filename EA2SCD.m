function [C_ea2synmf, P_ea2synmf, cost_ea2synmf] = EA2SCD(A, ncls, lambda, maxin, maxout, tol, epsilon)

%--------------------------------------------------------
% inputs:
% A: adjacent matrix
% ncls: no. of cmmunities
% parameters. maxin, maxout, tol, epsilon
%--------------------------------------------------------

%% initialization
C_init = 2 * sqrt(mean(mean(A)) / ncls) * rand(size(A, 1), ncls);
cost_ea2synmf = zeros(maxout,1);
cost_ea2synmf(1) = sum(sum((A - C_init * C_init').^2));

%% SymmetricNMF
C_ea2synmf = SymmetricNMF(A, C_init, epsilon, maxin);
cost_ea2synmf(2) = sum(sum( (A - C_ea2synmf * C_ea2synmf').^2) );

%% EA2SCD
for i = 1 : maxout
    
    P_ea2synmf = max ( ( ones(size(A)).* A - ones(size(A)).*(C_ea2synmf * C_ea2synmf')) /  (lambda-1), - ones(size(A)).*A);
    
    A_new = ones(size(A)).*(A + P_ea2synmf); 

    C_ea2synmf = SymmetricNMF(A_new, C_ea2synmf, epsilon, 10);
    
    A_cur = C_ea2synmf * C_ea2synmf';
    
    cost_ea2synmf(i+2) = sum(sum( (A - A_cur).^2) );
    
    costerror = abs(cost_ea2synmf(i+1)-cost_ea2synmf(i+2));
    
    if costerror < tol
        disp(['EA2SCD Reach inner stopping criterion at out ' num2str(i+2)]);
        break;
    end
end