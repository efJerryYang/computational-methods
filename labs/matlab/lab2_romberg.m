function romberg(f, xlim, n, eps)
    a = xlim(1);
    b = xlim(2);
    h = b - a;
    T = zeros(n, n);
    T(1, 1) = 1 / 2 * h * (f(a) + f(b));
    for i = 1:n
        tmpsum = 0;
        jmax = 2^(i - 1);
        for j = 1:jmax
            tmpsum = tmpsum + f(a + (j - 1 / 2) * h);
        end
        T(i+1, 1) = 1 / 2 * T(i, 1) + 1 / 2 * h * tmpsum;

        for m = 1:i
            T(i+1, m+1) = (4^m * T(i+1, m) - T(i, m)) / (4^m - 1);
        end
        for m = 1:i
            fprintf("%12.9f\t", T(i, m));
        end
        fprintf("\n");
        if i > 1 && abs(T(i+1, i+1) - T(i, i)) < eps
            fprintf("Accuracy requirement satisfied.\n\n");
            break;
        end
        h = h / 2;
    end
end

% Problem 1.1: f(x) = x^2 * exp(x)
f = @(x) x^2 * exp(x);
eps = 1e-6;
xlim = [0, 1];
iter_num = 30;
fprintf("\nProblem 1.1: f(x) = x^2 * exp(x)\n");
romberg(f, xlim, iter_num, eps);

% Problem 1.2: f(x) = exp(x)sin(x)
f = @(x) exp(x) * sin(x);
eps = 1e-6;
xlim = [1, 3];
fprintf("Problem 1.2: f(x) = exp(x)sin(x)\n");
romberg(f, xlim, iter_num, eps);

% Problem 1.3: f(x) = 4 / (1 + x^2)
f = @(x) 4 / (1 + x^2);
eps = 1e-6;
xlim = [0, 1];
fprintf("Problem 1.3: f(x) = 4 / (1 + x^2)\n");
romberg(f, xlim, iter_num, eps);

% Problem 1.4: f(x) = 1 / (x + 1)
f = @(x) 1 / (x + 1);
eps = 1e-6;
xlim = [0, 1];
fprintf("Problem 1.4: f(x) = 1 / (x + 1)\n");
romberg(f, xlim, iter_num, eps);
