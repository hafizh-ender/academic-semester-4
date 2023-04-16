% Gauss Eliminination Program
% by Hafizh Renanto Akhmad (13621060)

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

% Initialize augmented matrix
Ab = [A b];

% Forward elminiation
n = length(A);

i = 1;
unsolvable = false;

while ((i < n) && ~unsolvable)
    pembagi = Ab(i, i);

    if (pembagi == 0)
        % Exchange row if pivot = 0
        dapet_pivot = false;
        pivot = i + 1;

        while (pivot <= n)
            nilai_pivot = Ab(pivot, i);

            if (nilai_pivot ~= 0)
                dapet_pivot = true;
                break
            else
                pivot = pivot + 1;
            end
        end

        if dapet_pivot
            disp("Switch Row")
            Ab = switch_row(Ab, i, pivot);
        else
            unsolvable = true;
        end
    else
        pengali = zeros(n, 1);
    
        for k = (i+1):n
            pengali(k) = Ab(k, i);
        end
    
        for j = i:(n+1)
            for k = (i+1):n
                Ab(k, j) = Ab(k, j) - Ab(i, j) * pengali(k) / pembagi;
            end
        end
    
        i = i + 1;
    end
end

% Tampilkan hasil matrix augmented
disp("Result of Forward Elimination (Augmented Matrix):")
Ab

if unsolvable || Ab(n, n) == 0
    disp("Unsolvable")
else
    % Back Substitution
    x = Ab(:, n+1);
    
    x(n) = x(n) / Ab(n, n);
    
    for i = (n-1):-1:1
        pengurang = 0;
    
        for j = n:-1:(i+1)
            pengurang = pengurang + Ab(i, j) * x(j);
        end
        
        x(i) = (x(i) - pengurang) / Ab(i, i);
    end
    
    % Tampilkan x
    disp("Result of Back Substitution:")
    x
end

function M = switch_row(A, a, b)

% Store the initial row a
N = A(a,:);
M = A(:,:);

% Replace value in row a with value in row b
for i = 1:length(N)
    M(a, i) = M(b, i);
end

% Replace value in row b with store value in row a
for i = 1:length(N)
    M(b, i) = N(i);
end

end