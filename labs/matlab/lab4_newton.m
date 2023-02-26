% multi-root newton method
function [n, x0] = newton(f, df, eps1, eps2, N, x0, lambda)
    n = 1;
    while n <= N
        F = f(x0);
        DF = df(x0);
        if abs(F) < eps1
            return
        end
        if abs(DF) < eps2
            fprintf('Reach a critical point!\n');
            return
        end
        x1 = x0 - lambda * F / DF;
        Tol = abs(x1 - x0);
        if Tol < eps1
            return
        end
        n = n + 1;
        x0 = x1;
    end
    fprintf('Fail to converge within %d iterations!\n', N);
end

% newton method
function [n, x0] = newton(f, df, eps1, eps2, N, x0)
    [n, x0] = newton(f, df, eps1, eps2, N, x0, 1.0);
end

function [f, df] = redefine_func(f, df)
    f = @(r, x) f(x);
    df = @(J, x) diag(df(x));
end
