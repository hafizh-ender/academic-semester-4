% Gauss-Seidel
A = [ 10  -2  -1   2   3   1  -4   7;
       5  11   3  10  -3   3   3  -4;
       7  12   1   5   3 -12   2   3;
       8   7  -2   1   3   2   2   4;
       2 -15  -1   1   4  -1   8   3;
       4   2   9   1  12  -1   4   1;
      -1   4  -7  -1   1   1  -1  -3;
      -1   3   4   1   3  -4   7   6
    ];

b = [  0;
      12;
      -5;
       3;
     -25;
     -26;
       9;
      -7;
    ];

n = length(b);

% Inisialisasi array of x
x = ones(n, 1);
x_simpan = ones(n, 1);

% Inisialisasi jumlah iterasi dan batas iterasi
iter = 0;
stop_iter = 50;

finish = false;

while(~finish)
    for i = 1:n
        sum_to_mult = 0;

        for j = 1:n
            if (i ~= j)
                sum_to_mult = sum_to_mult + A(i, j) * x(j);
            end
        end

        x_simpan(i) = (1 / A(i, i)) * (b(i) - sum_to_mult);
    end

    iter = iter + 1;
    x = x_simpan;

    if (iter >= stop_iter)
        finish = true;
    end

    fprintf("Iterasi ke-%d:\n x = \n", iter);
    disp(x);
end