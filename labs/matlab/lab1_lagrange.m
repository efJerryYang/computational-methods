function [x, y] = lagrange(xs, fxs, x)
    num = size(xs, 1);
    y = 0;
    for i = 1:num
        li = 1;
        for j = 1:num
            if j == i
                continue
            end
            li = li .* (x - xs(j)) ./ (xs(i) - xs(j));
        end
        y = y + li * fxs(i);
    end
    x = x * ones(size(y));
end

function [x, y] = lagrange(xs, fxs, x_vec)
    num = size(xs, 1);
    y = zeros(size(x_vec));
    for i = 1:num
        li = ones(size(x_vec));
        for j = 1:num
            if j == i
                continue
            end
            li = li .* (x_vec - xs(j)) ./ (xs(i) - xs(j));
        end
        y = y + li .* fxs(i);
    end
    x = x_vec;
end
