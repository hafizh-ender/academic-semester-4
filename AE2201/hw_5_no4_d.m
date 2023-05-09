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

% Get the lower-triangular matrix
L = choleski(B)

n = 8;

% Get the H
H = (L^-1) * A * transpose(L^-1)

% Solve with Jacobi
[eVals, Z] = jacobi(H)
X = transpose(L^-1) * Z;
for i = 1:n % Normalize eigenvectors.
    xMag = sqrt(dot(X(:,i),X(:,i)));
    X(:,i) = X(:,i)/xMag;
end

[eVals,X] = sortEigen(eVals,X); % Sort in ascending order.

eigenvalues = eVals(1:3)' % Extract 3 smallest
eigenvectors = X(:,1:3) % eigenvalues & vectors.

n = length(eigenvectors);

% Plot the 3 lowest eigenvalue
x = 0:9;
x_1 = zeros(1,10);
x_2 = zeros(1,10);
x_3 = zeros(1,10);

for i = 1:8
    x_1(i+1) = X(i, 1);
end

for i = 1:8
    x_2(i+1) = X(i, 2);
end

for i = 1:8
    x_3(i+1) = X(i, 3);
end

clf
u1 = plot(x, x_1, "Marker","o");
hold on
u2 = plot(x, x_2, "Marker","square");
u3 = plot(x, x_3, "Marker","diamond");
grid on
xlabel('x')
ylabel('u')
ylim([-0.4 0.6])
legend([["0.2431" "0.6897" "1.2913"]])
hold off

%%
function L = choleski(matrix)

% Initiate lower-triangular matrix
n = length(matrix);
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

function [eVals, eVecs] = jacobi(A,tol)

% Jacobi method for computing eigenvalues and
% eigenvectors of a symmetric matrix A.
% USAGE: [eVals,eVecs] = jacobi(A,tol)
% tol = error tolerance (default is 1.0e-9).
if nargin < 2; tol = 1.0e-9; end
n = size(A,1);
maxRot = 5*(n^2); % Limit number of rotations
P = eye(n); % Initialize rotation matrix
for i = 1:maxRot % Begin Jacobi rotations
    [Amax,k,L] = maxElem(A);
    if Amax < tol
        eVals = diag(A);
        eVecs = P;
        return
    end
    [A,P] = rotate(A,P,k,L);
end

end

function [Amax,k,L] = maxElem(A)

% Finds Amax = A(k,L) (largest off-diag. elem. of A).
n = size(A,1);
Amax = 0;
for i = 1:n-1
    for j = i+1:n
        if abs(A(i,j)) >= Amax
            Amax = abs(A(i,j));
            k = i; L = j;
        end
    end
end

end

function [A,P] = rotate(A,P,k,L)

% Zeros A(k,L) by a Jacobi rotation and updates
% transformation matrix P.
n = size(A,1);
diff = A(L,L) - A(k,k);
if abs(A(k,L)) < abs(diff)*1.0e-36
    t = A(k,L);
else
    phi = diff/(2*A(k,L));
    t = 1/(abs(phi) + sqrt(phi^2 + 1));

    if phi < 0
        t = -t;
    end
end
c = 1/sqrt(t^2 + 1); s = t*c;
tau = s/(1 + c);
temp = A(k,L); A(k,L) = 0;
A(k,k) = A(k,k) - t*temp;
A(L,L) = A(L,L) + t*temp;
for i = 1:k-1 % For i < k
    temp = A(i,k);
    A(i,k) = temp -s*(A(i,L) + tau*temp);
    A(i,L) = A(i,L) + s*(temp - tau*A(i,L));
end
for i = k+1:L-1 % For k < i < L
    temp = A(k,i);
    A(k,i) = temp - s*(A(i,L) + tau*A(k,i));
    A(i,L) = A(i,L) + s*(temp - tau*A(i,L));
end
for i = L+1:n % For i > L
    temp = A(k,i);
    A(k,i) = temp - s*(A(L,i) + tau*temp);
    A(L,i) = A(L,i) + s*(temp - tau*A(L,i));
end
for i = 1:n % Update transformation matrix
    temp = P(i,k);
    P(i,k) = temp - s*(P(i,L) + tau*P(i,k));
    P(i,L) = P(i,L) + s*(temp - tau*P(i,L));
end

end

function [eVals,eVecs] = sortEigen(eVals,eVecs)

% Sorts eigenvalues & eigenvectors into ascending order of eigenvalues.
% USAGE: [eVals,eVecs] = sortEigen(eVals,eVecs)
n = length(eVals);
for i = 1:n-1
    index = i; val = eVals(i);
    for j = i+1:n
        if eVals(j) < val
            index = j; val = eVals(j);
        end
    end
    if index ~= i
        eVals = switch_row(eVals,i,index);
        eVecs = switch_col(eVecs,i,index);
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

function M = switch_col(A, a, b)

% Store the initial col a
N = A(:,a);
M = A(:,:);

% Replace value in col a with value in col b
for i = 1:length(N)
    M(i, a) = M(i, b);
end

% Replace value in col b with store value in col a
for i = 1:length(N)
    M(i, b) = N(i);
end

end