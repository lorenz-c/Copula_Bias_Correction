function C = ecopula(U)
% Computes the empirical copula process according to the dimension
% rediction approach from Deheuvels (1979) and Genest and Rémillard (2008).

% References: - Genest and Rémillard (2008),
%             - Deheuvels (1979)


% Number of "time-steps"
n = size(U, 1);

% Number of variables
d = size(U, 2);

% Loop over the number of time-steps
for i = 1:n
    % Within this loop, we're evaluating an indicator function, which is 1
    % if all elements in a row of U are less or equal than the elements of 
    % a specific row u_i(:). 

    % Create a matrix where all elements, which are less than or equal the
    % elements of the ith row, are 1.
    ind = bsxfun(@le, U, U(i, :));
    % Sum up over all variables (i.e. the columns)
    ind = sum(ind, 2);
    % Set all rows to zero, where at least one element in ind is zero (i.e.
    % is not less than or equal u_ij). 
    ind(ind < d)  = 0;
    % The remaining elements are set to 1.
    ind(ind == d) = 1;
    % Count the number of 1s and divide the result by the number of
    % time-steps plus 1. 
    C(i, 1) = sum(ind)/(n+1);    
end

    
    