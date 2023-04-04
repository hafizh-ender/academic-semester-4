% Baca matriks A
A = [ 10  -2  -1   2   3   1  -4   7;
       5  11   3  10  -3   3   3  -4;
       7  12   1   5   3 -12   2   3;
       8   7  -2   1   3   2   2   4;
       2 -15  -1   1   4  -1   8   3;
       4   2   9   1  12  -1   4   1;
      -1   4  -7  -1   1   1  -1  -3;
      -1   3   4   1   3  -4   7   6
    ]

% Baca matriks b
b = [  0;
      12;
      -5;
       3;
     -25;
     -26;
       9;
      -7;
    ]

% Ambil panjang/dimensi dari b
n = length(b)

% Buat matriks augmented
Ab = [A b];

% Loop over all column (for pivot)
for i = 1:n
    % Simpan bilanga di pivot
    pembagi = Ab(i, i);

    % Simpan bilangan pengali
    pengali = zeros(n, 1);

    for k = 1:n
        pengali(k) = Ab(k, i);
    end

    % Lakukan ERO pada matriks
    for j = i:(n+1)
        for k = 1:n
            if (k ~= i)
                Ab(k, j) = Ab(k, j) - Ab(i, j) * pengali(k) / pembagi;
            end
        end
    end

    for j = 1:(n+1)
        Ab(i, j) = Ab(i, j) / pembagi;
    end
end

Ab