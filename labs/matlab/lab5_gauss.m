function X = swaprows(X, i, j)
    X([i j], :) = X([j i], :);
end

function [val, idx] = pivoting(A, k, n)
    [val, idx] = max(A(k:n, k));
    idx = idx + k - 1;  % index must add previous length that omitted by slice operator
end

function [val, idx] = pivoting_with_scaling(A, b, k, n, implicit)
    s = max(A(:, k:n), [], 2);
    if any(s == 0)
        fprintf('Cannot solve a singular matrix!\n');
        return
    end
    if implicit
        [val, idx] = max(A(k:n, k) ./ s(k:n));
    else
        A(k:n, k:n) = A(k:n, k:n) ./ s(k:n, k:n);
        b(k:n) = b(k:n) ./ s(k:n);
        [val, idx] = max(A(k:n, k));
    end
    idx = idx + k - 1;  % index must add previous length that omitted by slice operator
end

function x = gauss(n, A, b)
    for k = 1:n-1
        % select pivot in columns
        [val, idx] = pivoting(A, k, n);
        if val == 0
            fprintf('Cannot solve a singular matrix!\n');
            return
        end
        % swap rows
        if idx ~= k
            A = swaprows(A, idx, k);
            b([idx k]) = b([k idx]);
        end
        % elimination
        for i = k+1:n
            m = A(i, k) / A(k, k);
            A(i, :) = A(i, :) - A(k, :) * m;
            b(i) = b(i) - b(k) * m;
        end
    end
    if A(n, n) == 0
        fprintf('Cannot solve a singular matrix!\n');
        return
    end
    x = zeros(n, 1);
    x(n) = b(n) / A(n, n);
    for k = n-1:-1:1  % the usage of reverse sequence
        x(k) = (b(k) - dot(A(k, k+1:n), x(k+1:n))) / A(k, k);  % something really annoying 
    end
end
