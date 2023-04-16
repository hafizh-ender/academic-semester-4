% LU Decompisition: Choleski for Jacobian Rotation
% by Hafizh Renanto Akhmad (13621060)

% Define the matrices
A = [
      5,-4, 1, 0, 0, 0, 0, 0;
     -4, 6,-4, 1, 0, 0, 0, 0;
      1,-4, 6,-4, 1, 0, 0, 0;
      0, 1,-4, 6,-4, 1, 0, 0;
      0, 0, 1,-4, 6,-4, 1, 0;
      0, 0, 0, 1,-4, 6,-4, 1;
      0, 0, 0, 0, 1,-4, 6,-4;
      0, 0, 0, 0, 0, 1,-4, 7;
]

B = [
      2,-1, 0, 0, 0, 0, 0, 0;
     -1, 2,-1, 0, 0, 0, 0, 0;
      0,-1, 2,-1, 0, 0, 0, 0;
      0, 0,-1, 2,-1, 0, 0, 0;
      0, 0, 0,-1, 2,-1, 0, 0;
      0, 0, 0, 0,-1, 2,-1, 0;
      0, 0, 0, 0, 0,-1, 2,-1;
      0, 0, 0, 0, 0, 0,-1, 2;
]

n = length(B);

% Show L
L = choleski(B, n)

% Show inverse L
inverse_L = L^-1

% Show inverse L transpose
inverse_L_transpose = transpose(L^-1)

% Show H
H = (L^-1) * A * transpose(L^-1)

function L = choleski(matrix, n)

% Initiate lower-triangular matrix
L = zeros(n, n);

% Fill L
for i = 1:n
    for j = 1:(i - 1)
        sum_LL = 0;
    
        for k = 1:(j - 1)
            sum_LL = sum_LL + L(j, k) * L(i, k);
        end
    
        L(i, j) = (matrix(i, j) - sum_LL) / L(j, j);
    end

    sum_L2 = 0;

    for k = 1:(i - 1)
        sum_L2 = sum_L2 + (L(i, k))^2;
    end

    L(i, i) = sqrt(matrix(i,i) - sum_L2);
end

end

function transposed = transpose(matrix)

% Initialize the result
transposed = matrix(:,:);
n = length(transposed);

for i = 1:n
    for j = 1:n
        transposed(i, j) = matrix(j, i);
    end
end

end
