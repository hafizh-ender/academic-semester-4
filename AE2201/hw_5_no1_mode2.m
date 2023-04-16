% LU Decompisition: Choleski Program

% Define the matrices
disp("Ax = b")
A = [
     4, -2,  2;
    -2,  2, -4;
     2, -4, 11
]

b = [
    2;
    5;
    4
]

% Forward elminiation
n = length(A);

% Get the lower-triangular matrix
disp("Lower-Triangular Matrix:")
L = choleski(A, n)

% Define the upper triangular matrix
disp("Upper-Triangular Matrix:")
U = transpose(L)

% Substitution for y
y = b(:);

y(1) = y(1) / L(1, 1);

for i = 2:n
    pengurang = 0;

    for j = 1:(i-1)
        pengurang = pengurang + L(i, j) * y(j);
    end

    y(i) = (y(i) - pengurang) / L(i, i);
end

% Tampilkan y
disp("Result of Substitution for y")
y


% Back Substitution for x
x = y(:);

x(n) = x(n) / U(n, n);

for i = (n-1):-1:1
    pengurang = 0;

    for j = n:-1:(i+1)
        pengurang = pengurang + U(i, j) * x(j);
    end
    
    x(i) = (x(i) - pengurang) / U(i, i);
end

% Tampilkan x
disp("Result of Back Substitution for x")
x

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