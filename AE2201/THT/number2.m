%% Problem 2 THT AE2201 2022/2023
% by Hafizh Renanto Akhmad (13621060)

clc, clear

n = input("Input number of partition: ");
et = 1.0e-15;

% Get the A and B matrix for propped cantilever beam problem of
% Ax = lambda Bx
[A, B] = proppedCantileverAB(n);

% Get standardize form Hz = lambda z 
[H, T] = standardizeForm(A, B, n);

% Display H
H

% Get eigenvalues and eigenvector from H
[lambda, Z] = jacobi(H, n, et);

% Get the X and normalize
X = T * Z;

% Normalize x
for j = 1:n
    magX = 0;
    
    for i = 1:n
        magX = magX + (X(i, j))^2;
    end
    
    magX = sqrt(magX);

    for i = 1:n
        X(i, j) = X(i, j) / magX;
    end
end

% Sort eigenvalues and eigenvectors based on eigenvalues
[lambda, X] = sortEigenAscending(lambda, X, n);

% Get the three smalles eigenvalues for three smalles buckling loads
smallestLambda = lambda(1:3);
smallestLoadsCoeff = (n + 1)^2 * smallestLambda;
smallestX = X(:,1:3);

% Get the analytical solution
smallestLoadsAnalyticalCoeff = [20.19; 59.68; 118.9];

% Calculate the error
errors = abs((smallestLoadsCoeff - smallestLoadsAnalyticalCoeff) ./ smallestLoadsAnalyticalCoeff);

% Display the solution)
fprintf("Your three smallest eigenvalues with its eigenvectors, and the smallest buckling loads are\n\n")
smallestEigenvalue = smallestLambda'
smallestEigenvector = smallestX
smallestLoadsCoeff

% Plot the solution
x = (0:(n+1))';

plotX = zeros(n+2, 3);
for i = 1:n+2
    if i == 1 || i == n+2
        plotX(i,1:3) = 0;
    else
        plotX(i,1:3) = smallestX(i - 1,1:3);
    end
end

plot(x,plotX(:,1),"Marker","o");
hold on
grid on
plot(x,plotX(:,2),"Marker","square");
plot(x,plotX(:,3),"Marker","diamond");
xlim([0 n+1])
xlabel("Partition")
ylabel("u")

%% Functions
function [eigenvalues, eigenvectors] = sortEigenAscending(eigenvalues, eigenvectors, n)
for i = 1:n-1
    minIndex = i;
    minVal = eigenvalues(minIndex);

    for j = i+1:n
        if eigenvalues(j) < minVal
            minIndex = j;
            minVal = eigenvalues(minIndex);
        end
    end
    if minIndex ~= i
        eigenvalues = swapRow(eigenvalues, 1, i, minIndex);
        eigenvectors = swapCol(eigenvectors, n, i, minIndex);
    end
end
end

function newMatrix = swapRow(matrix, ncol, row1, row2)
% Store the matrix into new matrix
newMatrix = matrix;
% Apply swapping row
for i = 1:ncol
    newMatrix(row1, i) = matrix(row2, i);
    newMatrix(row2, i) = matrix(row1, i);
end
end

function newMatrix = swapCol(matrix, nrow, col1, col2)
% Store the matrix into new matrix
newMatrix = matrix;
% Apply swapping row
for i = 1:nrow
    newMatrix(i, col1) = matrix(i, col2);
    newMatrix(i, col2) = matrix(i, col1);
end
end

function [A, B] = proppedCantileverAB(n)
A = zeros(n, n);
B = zeros(n, n);

for i = 1:n
    if (i == 1)
        A(i, i) = 5;
    elseif (i == n)
        A(i, i) = 7;
    else
        A(i, i) = 6;
    end

    B(i, i) = 2;

    if (i > 1)
        A(i, i - 1) = -4;
        A(i - 1, i) = -4;
        B(i, i - 1) = -1;
        B(i - 1, i) = -1;
    end

    if (i < n)
        A(i, i + 1) = -4;
        A(i + 1, i) = -4;
        B(i, i + 1) = -1;
        B(i + 1, i) = -1;
    end

    if (i > 2)
        A(i, i - 2) = 1;
        A(i - 2, i) = 1;
    end

    if (i < n - 1)
        A(i, i + 2) = 1;
        A(i + 2, i) = 1;
    end
end
end

function [H, T] = standardizeForm(A, B, n)
% Standardize Ax = lBx to Hz = lz

% Get lower triangular Choleski matrix of B
L = choleski(B, n);

% Get the inverse of L
invL = inverse(L);

% Get the H and T
H = invL * (A * transpose(invL));
T = transpose(invL);
end

function inversed = inverse(matrix, n)
inversed = matrix^-1;
end

function transposed = transpose(matrix, n)
transposed = matrix';
end

function L = choleski(A, n)
% Initiate lower-triangular matrix
L = zeros(n, n);

% Fill L
for i = 1:n
    for j = 1:(i - 1)
        sum_LL = 0;
    
        for k = 1:(j - 1)
            sum_LL = sum_LL + L(j, k) * L(i, k);
        end
    
        L(i, j) = (A(i, j) - sum_LL) / L(j, j);
    end

    sum_L2 = 0;

    for k = 1:(i - 1)
        sum_L2 = sum_L2 + (L(i, k))^2;
    end

    L(i, i) = sqrt(A(i,i) - sum_L2);
end
end

function [eigenvalues, eigenvectors] = jacobi(A, n, eTol)
% Limit the number of rotations
maxRot = 5 * (n^2);

%  Initialize identity matrix for rotation
P = eye(n);

% Loop rotation
for i = 1:maxRot
    [maxEl, maxRow, maxCol] = maxElem(A, n);

    if maxEl < eTol
        eigenvalues = diagonal(A, n);
        eigenvectors = P;

        return
    end

    [A, P] = rotate(A, P, maxRow, maxCol, n);
end
end

function diagonalArray = diagonal(matrix, n)
diagonalArray = zeros(n, 1);

for i = 1:n
    diagonalArray(i) = matrix(i, i);
end
end

function [A, P] = rotate(A, P, maxRow, maxCol, n)
% Target: zeroes A(maxRow, maxCol)

% Calculate difference in diagonal
diagDiff = A(maxCol, maxCol) - A(maxRow, maxRow);

% Calculate tan theta (tanT) or tangent of angle of rotation
if abs(A(maxRow, maxCol)) < abs(diagDiff) * 1.0e-36
    t = A(maxRow, maxCol);
else
    phi = diagDiff / (2 * A(maxRow, maxCol));
    t = 1/(abs(phi) + sqrt(phi^2 + 1));
    
    if phi < 0
        t = -t;
    end
end

% Calculate sine and cosine of theta using trigonometric identity
c = 1/sqrt(t^2 + 1);
s = t * c;

% Calculate tau
tau = s/(1 + c);

% Start the rotation
temp = A(maxRow, maxCol);
A(maxRow, maxCol) = 0;

A(maxRow, maxRow) = A(maxRow, maxRow) - t * temp;
A(maxCol, maxCol) = A(maxCol, maxCol) + t * temp;

for i = 1:maxRow-1
    % Loop for row less than maxRow
    temp = A(i, maxRow);

    A(i, maxRow) = temp - s * (A(i, maxCol) + tau * temp);
    A(i, maxCol) = A(i, maxCol) + s * (temp - tau * A(i, maxCol));
end

for i = maxRow+1:maxCol-1                   
    % Loop for row more than maxRow and less than maxCol
    temp = A(maxRow, i);

    A(maxRow, i) = temp - s * (A(i, maxCol) + tau * A(maxRow, i));
    A(i, maxCol) = A(i, maxCol) + s * (temp - tau * A(i, maxCol));
end

for i = maxCol+1:n
    % Loop for row more than maxCol
    temp = A(maxRow,i);

    A(maxRow, i) = temp - s * (A(maxCol, i) + tau * temp);
    A(maxCol, i) = A(maxCol, i) + s * (temp - tau * A(maxCol,i));
end

for i = 1:n
    % Update result
    temp = P(i, maxRow);
    P(i, maxRow) = temp - s * (P(i, maxCol) + tau * P(i, maxRow));
    P(i, maxCol) = P(i, maxCol) + s * (temp - tau * P(i, maxCol));
end

end

function [maxA, maxRow, maxCol] = maxElem(A, n)
% Initialize max value
maxA = 0;

% Begin looping
for i = 1:n-1
    for j = i+1:n
        if abs(A(i, j)) >= maxA
            maxA = abs(A(i, j));
            maxRow = i;
            maxCol = j;
        end
    end
end

end