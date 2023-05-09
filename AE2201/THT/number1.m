%% Welcome Page
printWelcome();


%% Input the Equation
printChangeSection();
n = input("Insert number of equations: ");

printChangeSection();

inputChoice = input(sprintf("Please input matrix insertion process: \n1: Row based\n2: All at one (only accepting decimals)\n\nChoice: "));

printChangeSection();

if (inputChoice == 1)
    A = inputCoefficientPerLine(n);
    b = inputConstantPerLine(n);
else
    [A, b] = inputSemua(n);
end

printChangeSection();

printMiddle("Confirmation for System of Linear Equations", 80)
printEquation([A b], n);
%% Choose Solving Method
printChangeSection();

methodChoice = input(sprintf("Please choose solving method:\n1. Gauss Elimination-Back Substitution\n2. Gauss-Seidel\n3. Gauss-Jordan Elimination\n\nPilihan: "));

if (methodChoice == 1)
    [stateResult, strResult, xResult] = gaussElimination(A, b, n);
elseif (methodChoice == 2)
    printChangeSection();
    imax = input("Please input maximum iteration: ");
    es = input("Please input minimum approximate error: ");

    [stateResult, strResult, xResult] = gaussSeidel(A, b, n, imax, es);
else
    [stateResult, strResult, xResult] = gaussJordanElimination(A, b, n);
end

%% Show Result
printChangeSection();

fprintf("%s\n", strResult)
if (stateResult)
    printResult(xResult, n)
end

%% Testing function

%% Input Functions
function [A] = inputCoefficientPerLine(n)
A = zeros(n, n);

for i = 1:n
    for j = 1:n
        A(i, j) = input(sprintf("Insert coefficient element %d for equation %d (a%d%d): ", j, i, i, j));
    end
end
end

function [b] = inputConstantPerLine(n)
b = zeros(n, 1);

for i = 1:n
    b(i) = input(sprintf("Insert constant for equation %d (b%d): ", i, i));
end
end

function [A, b] = inputSemua(n)
printMiddle("Input in form of matrix", 80)
printMiddle("Use ';' to separate rows", 80)
printMiddle("and ',' to separate columns", 80)

inputMatrix = input("", "s");

try
    dataPerBaris = string(split(inputMatrix, ";"));
    
    dataStr = strings(n, n + 1);
    dataNum = zeros(n, n + 1);
    
    for i = 1:n
        dataStr(i, :) = string(split(dataPerBaris(i), ","));
    
        for j = 1:(n + 1)
            dataNum(i, j) = str2double(dataStr(i, j));
        end
    end
    
    A = dataNum(:,1:n);
    b = dataNum(:,n+1);
catch
    error("Invalid input. Please assign according to inputted number of equations.")

    A = 0;
    b = 0;
end

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

%% Gauss Seidel Functions
function [stateResult, strResult, xResult] = gaussSeidel(A, b, n, iterLim, es)
% Check diagonally dominance of the matrix
isDiagonallyDominant = checkDiagonallyDominant(A, n);
% isDiagonallyDominant = true;

% Do some row exchange if not diagonally dominant
% Limited only to 
if ~isDiagonallyDominant
    % Check every diagonal element except last one
    for i = 1:(n-1)
        % Check below the diagonal element
        for j = 1:n
            % If the abs value of element below the diagonal element is
            % larger than the abs value of current diagonal element and
            % value of element on our diagonal row but on that diagonal
            % column is not zero, do row swapping
            if (abs(A(j, i)) > abs(A(i, i)) && A(i, j))
                A = swapEro(A, n, i, j);
                b = swapEro(b, 1, i, j);
            end
        end
    end
end

% Check again and give warning based on diagonally dominance
isDiagonallyDominant = checkDiagonallyDominant(A, n);
if ~isDiagonallyDominant
    warning("The matrix is not diagonally dominant, potentially broke the algorithm")
end

% Initialize array of solution, a helping array, and array of errors
x = ones(n, 1);
xStore = ones(n, 1);
xErrors = ones(n, 1);

% Initialize iteration variable
iter = 0;
done = false;

while(~done)
    % Do the iterative scheme
    for i = 1:n
        sumAx = 0;

        for j = 1:n
            if (i ~= j)
                sumAx = sumAx + A(i, j) * x(j);
            end
        end

        xStore(i) = (1 / A(i, i)) * (b(i) - sumAx);
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

%% Formatting Procedures
function printWelcome()
fprintf("  __  __       _        _         ____      _            _       _             \n")
fprintf(" |  \\/  | __ _| |_ _ __(_)_  __  / ___|__ _| | ___ _   _| | __ _| |_ ___  _ __ \n")
fprintf(" | |\\/| |/ _` | __| '__| \\ \\/ / | |   / _` | |/ __| | | | |/ _` | __/ _ \\| '__|\n")
fprintf(" | |  | | (_| | |_| |  | |>  <  | |__| (_| | | (__| |_| | | (_| | || (_) | |   \n")
fprintf(" |_|  |_|\\__,_|\\__|_|  |_/_/\\_\\  \\____\\__,_|_|\\___|\\__,_|_|\\__,_|\\__\\___/|_|   \n")

fprintf("\nby: Hafizh Renanto Akhmad\n")
fprintf("      13621060\n")
end

function printMiddle(string, total)
chars = convertStringsToChars(string);

n = length(chars);

beforePrint = floor((total - n) / 2);

for i = 1:beforePrint
    fprintf(" ");
end
for i = 1:n
    fprintf(chars(i));
end

fprintf("\n")
end

function printRightInline(string, total)
chars = convertStringsToChars(string);

n = length(chars);

beforePrint = total - n;

for i = 1:beforePrint
    fprintf(" ");
end
for i = 1:n
    fprintf(chars(i));
end
end

function printEquation(Ab, n)
stringAb = matrixNumToStr(abs(Ab), n, n + 1);

maxSizeCol = zeros(1, n + 1);
for j = 1:(n+1)
    for i = 1:n
        if (j > 1)
            currentStringLength = length(char(stringAb(i, j)));
        else
            currentStringLength = length(char(sprintf("%g", Ab(i, j))));
        end

        if (currentStringLength > maxSizeCol(1, j))
            maxSizeCol(1, j) = currentStringLength;
        end
    end
end

for i = 1:n
    for j = 1:(n+1)
        if (j < n + 1)
            if (Ab(i, j))
                if (j == 1)
                    if ((abs(Ab(i, j)) ~= 1) || Ab(i, j) == 0)
                        stringToPrint = sprintf("%g",Ab(i, j));
                    else
                        stringToPrint = " ";
                    end
                else
                    absSumOfPreviousColumn = 0;

                    for k = 1:(j - 1)
                        absSumOfPreviousColumn = absSumOfPreviousColumn + abs(Ab(i, k));
                    end

                    if (absSumOfPreviousColumn)
                        if (Ab(i, j) > 0)
                            fprintf(" + ")
                        else
                            fprintf(" - ")
                        end
                    else
                        if (Ab(i, j) > 0)
                            fprintf("   ")
                        else
                            fprintf(" - ")
                        end
                        
                    end

                    if ((abs(Ab(i, j)) ~= 1) || Ab(i, j) == 0)       
                        stringToPrint = stringAb(i, j);
                    else
                        stringToPrint = " ";
                    end
                end

                printRightInline(stringToPrint, maxSizeCol(1, j));
                fprintf(" X%d", j)
            else
                nSpaces = maxSizeCol(1, j) + 3;

                if (~(j == 1))
                    nSpaces = 3 + nSpaces;
                end

                for k = 1:nSpaces
                    fprintf(" ")
                end
            end
        else
            fprintf(" = %g\n", Ab(i, j))
        end
    end
end

end

function printChangeSection()
for i=1:80
    fprintf(">")
end
fprintf("\n")
end

function matrixString = matrixNumToStr(matrixNum, nrow, ncol)
matrixString = strings(nrow, ncol);
for i = 1:nrow
    for j = 1:ncol
        matrixString(i, j) = sprintf("%g", matrixNum(i, j));
    end
end
end

function printResult(x, n)
for i = 1:n
    fprintf("X%d = %g\n", i, x(i, 1));
end
end
