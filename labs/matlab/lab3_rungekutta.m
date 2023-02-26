function [xs, ys] = rungekutta(f, xspan, y0, num)
    a = xspan(1);
    b = xspan(2);
    x0 = a;
    h = (b - a) / num;
    xs = zeros(num, 1);
    ys = zeros(num, 1);
    for n = 1:num
        K1 = h * f(x0, y0);
        K2 = h * f(x0 + h / 2, y0 + K1 / 2);
        K3 = h * f(x0 + h / 2, y0 + K2 / 2);
        K4 = h * f(x0 + h, y0 + K3);
        x1 = x0 + h;
        y1 = y0 + 1 / 6 * (K1 + 2 * K2 + 2 * K3 + K4);
        xs(n) = x0;
        ys(n) = y0;
        x0 = x1;
        y0 = y1;
    end
end

function [p, xs, ys] = show_plot(p, f, xspan, y0, iternum)
    [xs, ys] = rungekutta(f, xspan, y0, iternum);
    p = plot(xs, ys);
end

function [p, xs, ys] = show_plot(p, f, xs, show, text)
    ys = f(xs);
    if show
        p = plot(xs, ys);
        title(text);
    end
end

function show_result(f1, f2, f3, xspan, y0, iternums, show, dense, title, text)
    disp(title);
    for iternum = iternums
        fprintf("\nIternum: %d\n", iternum);
        p = 0;
        [p, xs, ys] = show_plot(p, f2, xspan, y0, iternum);
        [p, xt, yt] = show_plot(p, f3, xs, show, text);
        data = [xt yt ys];
        header = {'x', 'True y', 'Pred y'};
        disp(table(data(:, 1), data(:, 2), data(:, 3), 'VariableNames', header));
    end
end
