% Gauss Seidel Program
% by Hafizh Renanto Akhmad (13621060)

% Variable for the problem
s = 0.8;
c = 0.6;

% Define the matrices
disp("Ax = b")
A = [
    c,  1,   0, 0, 0;
    0,  s,   0, 0, 1;
    0,  0, 2*s, 0, 0;
    0, -c,   c, 1, 0;
    0,  s,   s, 0, 0
]

b = [
    0;
    0;
    1;
    0;
    0
]

n = length(A);

% First, check if the matrix is diagonally dominant or not
is_diagonal = is_diagonal_dominant(A);
fprintf("Checking Diagonally Dominance: %g\n", is_diagonal)

A

b

% As we can see, the matrix is not diagonally dominant. We can see that the
% value of diagonal element in 5th row is zero. Let's exchange it with row
% 2. Also, we can see that at row 1, the diagonal element is smaller. So
% let's do an ERO by substracting with the value of row 2 from row 1.
A_changed = A(:,:);
A_changed = switch_row(A_changed, 2, 5);
A_changed = ero(A_changed, 1, 2, -1);

% Check again
is_diagonal_second = is_diagonal_dominant(A_changed);
fprintf("Checking Diagonally Dominance After Changing Manually: %g\n", is_diagonal_second)

% Now, it's diagonally dominant, thus we may safely do Gaus-Seidel
% iteration. Apply the same to b
A = A_changed(:,:);
b = switch_row(b, 2, 5);
b = ero(b, 1, 2, -1);

% Initialize P as initial guess and the error array
P = ones(n, 1);
P_simpan = ones(n, 1);
P_error = ones(n, 1);

% Declare iteration limit and stop error
iter = 0;
stop_iter = 100;
stop_error = 0.01;

finish = false;

while(~finish)
    % Do the iterative scheme
    for i = 1:n
        sum_to_AP = 0;

        for j = 1:n
            if (i ~= j)
                sum_to_AP = sum_to_AP + A(i, j) * P(j);
            end
        end

        P_simpan(i) = (1 / A(i, i)) * (b(i) - sum_to_AP);
    end

    % Increment iteration
    iter = iter + 1;

    % Store approximate error
    for i = 1:n
        P_error(i) = (P_simpan(i) - P(i)) / P_simpan(i);
    end

    % Update P
    P = P_simpan;

    if (iter >= stop_iter || max(abs(P_error)) < stop_error)
        finish = true;
    end
end

% Show the result
P
P_error

function result = is_diagonal_dominant(matrix)

n = length(matrix);

sum = zeros(n, 1);

for i = 1:n
    for j = 1:n
        if (i ~= j)
            sum(i) = sum(i) + abs(matrix(i, j));
        end
    end
end

result = true;

for k = 1:n
    if (abs(matrix(k, k)) < sum(k))
        result = false;
    end
end

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

function new_matrix = ero(matrix, row_target, row_used, multiplier)

new_matrix = matrix(:,:);

columns = size(matrix, 2);

for i = 1:columns
    new_matrix(row_target, i) = new_matrix(row_target, i) + new_matrix(row_used, i) * multiplier;
end

end