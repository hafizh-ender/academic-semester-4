%% Main Program
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
n = 15;

% Initialize stop conditions
stopIter = 1000000000;
stopError = 0.01;

% Intialize overrelaxation
lambda = 1.5;

% Get the temperature matrix
[liebmanT, liebmanE] = liebmann(n, upperBC, lowerBC, leftBC, rightBC, stopIter, stopError, lambda);
gaussSeidelT = gaussSeidelSPL(upperBC, lowerBC, leftBC, rightBC, stopIter, stopError, 1);

% Flip for display
flippedLiebmanT = zeros(n + 2, n + 2);
for i = 1:(n + 2)
    flippedLiebmanT(n + 3 - i, :) = liebmanT(i, :);
end

% Plot the temperature matrix
plotChoose = input(sprintf("Please choose your plot below\n1: Heatmap\n2: 3D Surface\n3: 3D Bar\n\nChoice: "));
if plotChoose == 1
    fig = figure();
    ax = axes(fig);
    X = 0:n+1;
    Y = 0:n+1;
    h = heatmap(X,Y,liebmanT);
    h.YDisplayData=flip(h.YDisplayData);
    h.Title = "Temperature Distribution Heatmap in Heated Pipe (°C)";
elseif plotChoose == 2
    fig = figure();
    ax = axes(fig);
    X = 0:n+1;
    Y = 0:n+1;
    h = surf(X,Y,liebmanT);
    zlabel('Temperature (°C)');
    title('Temperature Distribution in Heated Pipe');
else
    fig = figure();
    ax = axes(fig);
    X = 0:n+1;
    Y = 0:n+1;
    h = bar3(Y,liebmanT);
    set(ax,'XTickLabel',X);
    set(ax,'YTickLabel',Y);
    for k = 1:n+1
        zdata = h(k).ZData;
        h(k).CData = zdata;
        h(k).FaceColor = 'interp';
    end
    zlabel('Temperature (°C)');
    title('Temperature Distribution in Heated Pipe');
end




%% Liebman Method
function [resultT, e] = liebmann(n, upperBC, lowerBC, leftBC, rightBC, stopIter, stopError, lambda)
% Initialize T
T = zeros(n + 2, n + 2);

% Fill the boundary conditions
for i = 2:(n+1)
    T(1, i)     = lowerBC;
    T(n + 2, i) = upperBC;
    T(i, 1)     = leftBC;
    T(i, n + 2) = rightBC;
end

% % Insert corner as half BC (for pretty plot purpose)
% T(1, 1)         = (T(1, 2) + T(2, 1)) / 2;
% T(n + 2, n + 2) = (T(n + 2, n + 1) + T(n + 1, n + 2)) / 2;
% T(n + 2, 1)     = (T(n + 2, 2) + T(n + 1, 1)) / 2;
% T(1, n + 2)     = (T(1, n + 1) + T(2, n + 2)) / 2;

% Insert corner as null (for pretty plot purpose)
T(1, 1)         = NaN;
T(n + 2, n + 2) = NaN;
T(n + 2, 1)     = NaN;
T(1, n + 2)     = NaN;

% Initialize storeT (for storing previous value)
storeT = T;

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
            T(i, j) = lambda * T(i, j) + (1 - lambda) * storeT(i, j);
            
            e(i - 1, j - 1) = abs((T(i, j) - storeT(i, j)) / T(i, j));
            
            % Replace errorMaxElement if larger error is shown
            if errorMaxElement < e(i - 1, j - 1)
                errorMaxElement = e(i - 1, j - 1);
            end
        end
    end

    % Update storeT
    for i = 2:(n + 1)
        for j = 2:(n + 1)
            storeT(i, j) = T(i, j);
        end
    end

    % Increment iter
    iter = iter + 1;

    % Replace max error with max error in the iteration
    if maxError > errorMaxElement
        maxError = errorMaxElement;
    end
end

resultT = T;
end

%% Gauss Seidel for SPL
function resultT = gaussSeidelSPL(upperBC, lowerBC, leftBC, rightBC, stopIter, stopError, lambda)
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

[~, ~, tResult] = gaussSeidel(T, T0, 9, stopIter, stopError, lambda);

resultT = zeros(5, 5);

for i = 2:4
    for j = 2:4
       resultT(i, 1) = leftBC;
       resultT(i, 5) = rightBC;
       resultT(1, i) = lowerBC;
       resultT(5, i) = upperBC;
    end
end

for i = 2:4
    resultT(2, i) = tResult(i - 1, 1); 
    resultT(3, i) = tResult(i + 2, 1); 
    resultT(4, i) = tResult(i + 5, 1); 
end

resultT(1, 1) = NaN;
resultT(5, 5) = NaN;
resultT(5, 1) = NaN;
resultT(1, 5) = NaN;
end

function [stateResult, strResult, xResult] = gaussSeidel(A, b, n, iterLim, es, lambda)
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

        xStore(i) = (lambda / A(i, i)) * (b(i) - sumAx) + (1 - lambda) * x(i);
        % xStore(i) = (1 / A(i, i)) * (b(i) - sumAx);
        
        % Apply overrelaxation
        % xStore(i) = lambda * xStore(i) + (1 - lambda) * x(i);
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