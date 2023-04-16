% LU Decompisition: Choleski for Jacobian Rotation

% Define the matrices\
A = [
     6, -4,  1,  0;
    -4,  6, -4,  1;
     1, -4,  6, -4;
     0,  1, -4,  7
]

B = [
     1, -2,  3, -1;
    -2,  6, -2,  3;
     3, -2,  6, -2;
    -1,  3, -2,  9
]

n = length(A);

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
