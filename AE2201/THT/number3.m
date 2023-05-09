%% Main Program
clc, clear

nim = input("Insert your NIM: ", "s");
charNim = char(nim);

% Parse the NIM into BCs
rightBC = str2double(charNim(4:5));
lowerBC = str2double(charNim(6:7));
leftBC  = str2double(charNim(8));
upperBC = str2double(charNim(1:3));

% rightBC = 50;
% lowerBC = 0;
% leftBC  = 75;
% upperBC = 100;

% Initialize number of point
n = input("Insert number of n (will result in n x n temperature matrix): ");
n = [n 0];

if ~n
    n = 3;
else
    n = n(1);
end

% Initialize stop conditions
stopIter = 1000000000;
stopError = 0.01;

% Intialize overrelaxation
lambda = 1.5;

% Get the temperature matrix
[liebmanT, liebmanE] = liebmann(n, upperBC, lowerBC, leftBC, rightBC, stopIter, stopError, lambda);
[gaussElimT, gaussSeidelT] = SPL(upperBC, lowerBC, leftBC, rightBC, stopIter, stopError, lambda);

% Flip for display
flippedLiebmanT = zeros(n + 2, n + 2);
for i = 1:(n + 2)
    flippedLiebmanT(n + 3 - i, :) = liebmanT(i, :);
end

flippedGaussSeidelT = zeros(5, 5);
for i = 1:5
    flippedGaussSeidelT(6 - i, :) = gaussSeidelT(i, :);
end

flippedGaussElimT = zeros(5, 5);
for i = 1:5
    flippedGaussElimT(6 - i, :) = gaussElimT(i, :);
end

% Plot the temperature matrix
fprintf("\nLiebman: \n")
disp(flippedLiebmanT)
fprintf("\nGauss-Seidel from Problem 1: \n")
disp(flippedGaussSeidelT)
fprintf("\nGauss Elimination from Problem 1: \n")
disp(flippedGaussElimT)

fig1 = figure();
X = 0:n+1;
Y = 0:n+1;
h = heatmap(X,Y,liebmanT);
h.YDisplayData=flip(h.YDisplayData);
h.Title = "Temperature Distribution Heatmap in Heated Pipe (°C)";

fig2 = figure();
X = 0:n+1;
Y = 0:n+1;
s = surf(X,Y,liebmanT);
zlabel('Temperature (°C)');
title(sprintf('Temperature Distribution in Heated Pipe, n = %d', n));

fig3 = figure();
ax = axes(fig3);
X = 0:n+1;
Y = 0:n+1;
b = bar3(Y,liebmanT);
set(ax,'XTickLabel',X);
set(ax,'YTickLabel',Y);
for k = 1:n+1
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
end
zlabel('Temperature (°C)');
title(sprintf('Temperature Distribution in Heated Pipe, n = %d', n));


%% Liebman Method
function [T, e] = liebmann(n, upperBC, lowerBC, leftBC, rightBC, stopIter, stopError, lambda)
% Initialize T
T = zeros(n + 2, n + 2);

% Fill the boundary conditions
for i = 2:(n+1)
    T(1, i)     = lowerBC;
    T(n + 2, i) = upperBC;
    T(i, 1)     = leftBC;
    T(i, n + 2) = rightBC;
end

% Insert corner as null (for pretty plot purpose)
T(1, 1)         = NaN;
T(n + 2, n + 2) = NaN;
T(n + 2, 1)     = NaN;
T(1, n + 2)     = NaN;

% Initialize storeT (for storing previous value)
prevT = T;

% Initialize iter values and error matrix
iter = 0;
e = ones(n, n);
maxError = 1;

% Do the loop
while ((iter < stopIter) && (maxError > stopError))
    % Initialize max error for an iteration
    errorMaxElement = 0;

    % Fill the temperature matrix with Liebmann method and the error matrix
    % with the approximate error
    for j = 2:(n + 1)
        for i = 2:(n + 1)
            T(i, j) = (T(i + 1, j) + T(i - 1, j) + T(i, j + 1) + T(i, j - 1)) / 4;
            T(i, j) = lambda * T(i, j) + (1 - lambda) * prevT(i, j);
            
            e(i - 1, j - 1) = abs((T(i, j) - prevT(i, j)) / T(i, j));
            
            % Replace errorMaxElement if larger error is shown
            if errorMaxElement < e(i - 1, j - 1)
                errorMaxElement = e(i - 1, j - 1);
            end
        end
    end

    % Update storeT
    for i = 2:(n + 1)
        for j = 2:(n + 1)
            prevT(i, j) = T(i, j);
        end
    end

    % Increment iter
    iter = iter + 1;

    % Replace max error with max error in the iteration
    if maxError > errorMaxElement
        maxError = errorMaxElement;
    end
end
end

%% Gausses for SPL.
function [gaussT, gaussSeidelT] = SPL(upperBC, lowerBC, leftBC, rightBC, stopIter, stopError, lambda)
% Buld up the coefficient matrix
T = zeros(9, 9);

for i = 1:9
    T(i, i) = 4;
    if (i < 9)
        if (mod(i, 3))
            T(i, i + 1) = -1;
            T(i + 1, i) = -1;
        end
    end
    if (i < 7)
        T(i, i + 3) = -1;
        T(i + 3, i) = -1;
    end
    if (i > 1)
        if (mod(i,3) ~= 1)
            T(i, i - 1) = -1;
            T(i - 1, i) = -1;
        end
    end
    if (i > 3)
        T(i, i - 3) = -1;
        T(i - 3, i) = -1;
    end
end

% Buld up the constant array by inserting boundary condition
T0 = zeros(9, 1);

for i = 1:9
    if (mod(i, 3) == 1)
        T0(i, 1) = T0(i, 1) + leftBC;
    elseif (mod(i, 3) == 0)
        T0(i, 1) = T0(i, 1) + rightBC;
    end

    if (i < 4)
        T0(i, 1) = T0(i, 1) + lowerBC;
    elseif (i > 6)
        T0(i, 1) = T0(i, 1) + upperBC;
    end
end

% Use Gauss-Seidel from number 1 to solve
[~, ~, gaussSeidelTResult] = gaussSeidel(T, T0, 9, stopIter, stopError, lambda);
[~, ~, gaussTResult] = gaussElimination(T, T0, 9);

% Build up result matrix from array of solutions
gaussT = zeros(5, 5);
gaussSeidelT = zeros(5, 5);

for i = 2:4
    for j = 2:4
       gaussSeidelT(i, 1) = leftBC;
       gaussSeidelT(i, 5) = rightBC;
       gaussSeidelT(1, i) = lowerBC;
       gaussSeidelT(5, i) = upperBC;

       gaussT(i, 1) = leftBC;
       gaussT(i, 5) = rightBC;
       gaussT(1, i) = lowerBC;
       gaussT(5, i) = upperBC;
    end
end

for i = 2:4
    gaussSeidelT(2, i) = gaussSeidelTResult(i - 1, 1); 
    gaussSeidelT(3, i) = gaussSeidelTResult(i + 2, 1); 
    gaussSeidelT(4, i) = gaussSeidelTResult(i + 5, 1); 

    gaussT(2, i) = gaussTResult(i - 1, 1); 
    gaussT(3, i) = gaussTResult(i + 2, 1); 
    gaussT(4, i) = gaussTResult(i + 5, 1); 
end

% Insert corner as null (for pretty plot purpose)
gaussSeidelT(1, 1) = NaN;
gaussSeidelT(5, 5) = NaN;
gaussSeidelT(5, 1) = NaN;
gaussSeidelT(1, 5) = NaN;

gaussT(1, 1) = NaN;
gaussT(5, 5) = NaN;
gaussT(5, 1) = NaN;
gaussT(1, 5) = NaN;
end

%% Gauss Elmination Functions
function [a, b] = swap(a, b)
temp = a;
a = b;
b = temp;
end

function newMatrix = additionEro(matrix, ncol, usedRow, targetRow, multiplier)
% Store the matrix into new matrix
newMatrix = matrix;

% Apply ERO
for i = 1:ncol
    newMatrix(targetRow, i) = matrix(targetRow, i) + matrix(usedRow) * multiplier;
end
end

function newMatrix = swapEro(matrix, ncol, row1, row2)
% Store the matrix into new matrix
newMatrix = matrix;
% Apply swapping row
for i = 1:ncol
    newMatrix(row1, i) = matrix(row2, i);
    newMatrix(row2, i) = matrix(row1, i);
end
end

function newMatrix = multiplyEro(matrix, ncol, row, multiplier)
% Store the matrix into new matrix
newMatrix = matrix;
% Apply multiplying row
for i = 1:ncol
    newMatrix(row, i) = matrix(row, i) * multiplier;
end
end

function pivotMatrix = partialPivoting(matrix, nrow, ncol, currentRow, currentCol)
maxRow = currentRow;
elemetMaxInCol = abs(matrix(maxRow, currentCol));

for i = (currentRow+1):nrow
    if (elemetMaxInCol < abs(matrix(i, currentCol)))
        maxRow = i;
        elementMaxInCol = abs(matrix(maxRow, currentCol));
    end
end

if maxRow ~= currentRow
    pivotMatrix = swapEro(matrix, ncol, maxRow, currentRow);
else
    pivotMatrix = matrix;
end
end

function rowEchelonMatrix = toRowEchelonForm(matrix, nrow, ncol)
barisPivot = 1;
kolomPivot = barisPivot;
pivotedMatrix = matrix;
while(barisPivot < nrow)
    % Compute partial pivoting
    pivotedMatrix = partialPivoting(pivotedMatrix, nrow, ncol, barisPivot, kolomPivot);

    % Check pivot element value
    if (pivotedMatrix(barisPivot, kolomPivot))
        % Store pivot element value
        pivotValue = pivotedMatrix(barisPivot, kolomPivot);

        for i = (barisPivot + 1):nrow
            % Store multiplier for row i
            multiplier = pivotedMatrix(i, kolomPivot) / pivotValue;

            for j = kolomPivot:ncol
                % Do ERO for every column in row i from pivot column to the
                % end
                % pivotedMatrix = additionEro(pivotedMatrix, ncol, barisPivot, i, -multiplier);
                pivotedMatrix(i, j) = pivotedMatrix(i, j) - pivotedMatrix(barisPivot, j) * multiplier;
            end
        end

        barisPivot = barisPivot + 1;
        kolomPivot = barisPivot;
    else
        kolomPivot = kolomPivot + 1;
    end
end

rowEchelonMatrix = pivotedMatrix;
end

function reducedRowEchelonMatrix = toReducedRowEchelonForm(matrix, nrow, ncol)
barisPivot = 1;
kolomPivot = barisPivot;
pivotedMatrix = matrix;
while(barisPivot <= nrow)
    % Compute partial pivoting
    pivotedMatrix = partialPivoting(pivotedMatrix, nrow, ncol, barisPivot, kolomPivot);

    % Check pivot element value
    if (pivotedMatrix(barisPivot, kolomPivot))
        % Store pivot element value
        pivotValue = pivotedMatrix(barisPivot, kolomPivot);

        for i = 1:nrow
            if (i ~= barisPivot)
                % Store multiplier for row i
                multiplier = pivotedMatrix(i, kolomPivot) / pivotValue;
    
                for j = kolomPivot:ncol
                    % Do ERO for every column in row i
                    % pivotedMatrix = additionEro(pivotedMatrix, ncol, barisPivot, i, -multiplier);
                    pivotedMatrix(i, j) = pivotedMatrix(i, j) - pivotedMatrix(barisPivot, j) * multiplier;
                end
            end
        end

        for j = kolomPivot:ncol
            pivotedMatrix(barisPivot, j) = pivotedMatrix(barisPivot, j) / pivotValue;
        end

        barisPivot = barisPivot + 1;
        kolomPivot = barisPivot;
    else
        kolomPivot = kolomPivot + 1;
    end
end

reducedRowEchelonMatrix = pivotedMatrix;
end

function x = backSubstitution(Ab, n)
% Initialize x from b
x = Ab(:, n + 1);

for i = n:-1:1
    pengurang = 0;

    for j = n:-1:(i+1)
        pengurang = pengurang + Ab(i, j) * x(j);
    end
    
    x(i) = (x(i) - pengurang) / Ab(i, i);
end
end

function rank = calculateRank(matrixRowEchelon, nrow, ncol)
% Initialize sum of absolute value of elements in each row
absSum = zeros(nrow, 1);

% Calculate zero rows
for i = 1:nrow
    for j = 1:ncol
        absSum(i, 1) = absSum(i, 1) + abs(matrixRowEchelon(i, j));
    end
end

% Initialize rank
rank = nrow;

% Calculate value in absSum which are zero
for i = 1:nrow
    rank = rank - ~absSum(i);
end
end

function [stateResult, strResult, xResult] = gaussElimination(A, b, n)
% Transform the matrix into row echelon form
rowEchelonForm = toRowEchelonForm([A b], n, n + 1);

% Calculate the rank of the A matrix in row echelon form
rankMatrix = calculateRank(rowEchelonForm(:, 1:n), n, n);

% Validate the result based on the rank
if (rankMatrix == n)
    strResult = "Solution found";
    xResult = backSubstitution(rowEchelonForm, n);
    stateResult = true;
else
    strResult = "Solution not found";
    xResult = 0;
    stateResult = false;
end
end

function [stateResult, strResult, xResult] = gaussJordanElimination(A, b, n)
% Transform the matrix into row echelon form
reducedRowEchelonForm = toReducedRowEchelonForm([A b], n, n + 1);

% Calculate the rank of the A matrix in row echelon form
rankMatrix = calculateRank(reducedRowEchelonForm(:, 1:n), n, n);

% Validate the result based on the rank
if (rankMatrix == n)
    strResult = "Solution found";
    xResult = reducedRowEchelonForm(:, n + 1);
    stateResult = true;
else
    strResult = "Solution not found";
    xResult = 0;
    stateResult = false;
end
end

%% Gauss Seidel
function [stateResult, strResult, xResult] = gaussSeidel(A, b, n, iterLim, es, lambda)
% Check diagonally dominance of the matrix
isDiagonallyDominant = checkDiagonallyDominant(A, n);
% isDiagonallyDominant = true;

% Ask the user to let the system do some row exchange or want to do it
% themself, as the system is not perfect
if ~isDiagonallyDominant
    warning("The matrix is not diagonally dominant, potentially broke the algorithm")
    rearrangeSelf = input(sprintf("Please choose your action:\n1: Reinput rearranged matrix\n2: Let the system rearrange (NOT STABLE)\nOther: Let it be\n\nChoice:"));
else
    rearrangeSelf = 0;
end

% Do some row exchange if not diagonally dominant
% Limited only to 
if rearrangeSelf == 2
    % Check every diagonal element except last one
    for i = 1:(n-1)
        % Check below the diagonal element
        for j = 1:n
            % If the abs value of element below the diagonal element is
            % larger than the abs value of current diagonal element and
            % value of element on our diagonal row but on that diagonal
            % column is not zero, do row swapping
            if (abs(A(j, i)) >= abs(A(i, i)) && A(i, j))
                A = swapEro(A, n, i, j);
                b = swapEro(b, 1, i, j);
            end
        end
    end

    % Check again and give warning based on diagonally dominance
    printChangeSection();
    printMiddle("Row Swapping Result", 80)
    printEquation([A b], n);
elseif rearrangeSelf == 1
    printChangeSection();
    
    inputChoice = input(sprintf("Please input matrix insertion process: \n1: Element based\n2: All at one (only accepting decimals)\n\nChoice: "));
    
    printChangeSection();
    
    if (inputChoice == 1)
        A = inputCoefficientPerElement(n);
        b = inputConstantPerLine(n);
    else
        [A, b] = inputSemua(n);
    end
    
    printChangeSection();
    
    printMiddle("Confirmation for System of Linear Equations", 80)
    printEquation([A b], n);
end

% Initialize array of solution, a helping array, and array of errors
x = zeros(n,1);
xStore = zeros(n,1);
xErrors = ones(n, 1);

% Initialize iteration variable
iter = 0;
done = false;

% Normalize row with diagona
for i = 1:n
    diagVal = A(i, i);

    for j = 1:n
        A(i, j) = A(i, j) / diagVal;
    end

    b(i) = b(i) / diagVal;
end

while(~done)
    % Do the iterative scheme
    for i = 1:n
        sumAx = 0;

        for j = 1:n
            if (i ~= j)
                sumAx = sumAx + A(i, j) * xStore(j);
            end
        end

        xStore(i) = b(i) - sumAx;
        xStore(i) = lambda * xStore(i) + (1 - lambda) * x(i);
    end

    % Increment iteration
    iter = iter + 1;

    % Calculate approximate error
    for i = 1:n
        xErrors(i) = (xStore(i) - x(i)) / xStore(i);
    end

    % Update x
    x = xStore;

    if ((iter >= iterLim) || (max(abs(xErrors)) < es))
        done = true;
    end
end

xResult = x;

if ~(isnan(xResult) + isinf(xResult))
    strResult = "Result obtained";
    stateResult = true;
else
    strResult = "Result not obtained";
    stateResult = false;
end
end

function result = checkDiagonallyDominant(matrix, n)
% Initialize matrix of zeros for summation of absolute value of elements in
% the matrix
sum = zeros(n, 1);

% Iterate to fill the matrix of zeros
for i = 1:n
    for j = 1:n
        if (i ~= j)
            sum(i) = sum(i) + abs(matrix(i, j));
        end
    end
end

% Initialize result
result = true;

for i = 1:n
    if (abs(matrix(i, i)) < sum(i))
        % Set to false if a diagonal element is less than sum of absolute
        % of other element in the row
        result = false;
    end
end
end
