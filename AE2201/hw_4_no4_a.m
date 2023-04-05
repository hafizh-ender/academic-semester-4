% Baca matrix
A = [ 10  -2  -1   2   3   1  -4   7;
       5  11   3  10  -3   3   3  -4;
       7  12   1   5   3 -12   2   3;
       8   7  -2   1   3   2   2   4;
       2 -15  -1   1   4  -1   8   3;
       4   2   9   1  12  -1   4   1;
      -1   4  -7  -1   1   1  -1  -3;
      -1   3   4   1   3  -4   7   6
    ]

b = [  0;
      12;
      -5;
       3;
     -25;
     -26;
       9;
      -7;
    ]

n = length(b)

% Augmentasi matrix
Ab = [A b];

% Mulai looping
for i = 1:n
    pembagi = Ab(i, i);
    pengali = zeros(n, 1);

    for k = 1:n
        pengali(k) = Ab(k, i);
    end

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

% Lihat hasil akhir
Ab