
using Markdown

display(md"""
## 实验题目3 四阶龙格-库塔(Runge-Kutta)方法
""")
using LaTeXStrings
using PrettyTables

# ╔═╡ 27d3fe58-815f-4550-8553-ccaa93073eff
function rungekutta(f::Function, xspan, y0, num)
    a, b = xspan
    x0 = a
    h = (b - a) / num
    xs, ys = zeros(num), zeros(num)
    for n = 1:num
        K1 = h * f(x0, y0)
        K2 = h * f(x0 + h / 2, y0 + K1 / 2)
        K3 = h * f(x0 + h / 2, y0 + K2 / 2)
        K4 = h * f(x0 + h, y0 + K3)
        x1 = x0 + h
        y1 = y0 + 1 / 6 * (K1 + 2K2 + 2K3 + K4)
        xs[n], ys[n] = x0, y0 = x1, y1
    end
    xs, ys
end

function show_plot(p, f::Function, xspan, y0::Float64, iternum::Integer)
    xs, ys = rungekutta(f, xspan, y0, iternum)
    p, xs, ys
end
function show_plot(p, f::Function, xs, show::Bool, text)
    p, xs, f.(xs)
end
function show_result(f1::Function, f2::Function, f3::Function, xspan, y0, iternums, show::Bool, dense::Bool, title, text)
    println("\n\n" * title)
    for iternum in iternums
        print("\nIternum: $iternum\n")
        p = 0
        p, xs, ys = show_plot(p, f2, xspan, y0, iternum)
        p, xt, yt = show_plot(p, f3, xs, show, text)
        data = [xt yt ys]
        header = (["x", "True y", "Pred y"])
        pretty_table(
            data;
            alignment=[:c, :c, :c],
            header=header,
            formatters=ft_printf("%14.8f"))
    end
end


# ╔═╡ ca710b63-74dc-4801-b353-f7ac826ffe5d
iternums = [5, 10, 20]

f1(y, p, x) = x + y    # lib RK4() solver
xspan = (0.0, 1.0)
y0 = -1.0
f2(x, y) = x + y       # my rungekutta() solver
f3(x) = -x - 1         # true result 
title = L"Problem\ 1.1: \frac{\mathrm{d} y}{\mathrm{d} x} = x + y"
text = L"y = -x - 1"
show_result(f1, f2, f3, xspan, y0, iternums, true, true, title, text) # show=true, dense=true

# ╔═╡ 58f8d978-2a36-4dfd-8c2d-590a508e8dce
iternums = [5, 10, 20]

f1(y, p, x) = -y^2
xspan = (0.0, 1.0)
y0 = 1.0
f2(x, y) = -y^2
f3(x) = 1 / (x + 1)
title = L"Problem\ 1.2: \frac{\mathrm{d} y}{\mathrm{d} x} = -y^2"
text = L"y = \frac{1}{x + 1}"
show_result(f1, f2, f3, xspan, y0, iternums, true, false, title, text) # show=true, dense=true

# ╔═╡ e548f440-96a4-449e-95e4-3adabf1006dc
iternums = [5, 10, 20]

f1(y, p, x) = 2 * y / x + x^2 * exp(x)
xspan = (1.0, 3.0)
y0 = 0.0
f2(x, y) = 2 * y / x + x^2 * exp(x)
f3(x) = x^2 * (exp(x) - exp(1))
title = L"Problem\ 2.1:\frac{\mathrm{d} y}{\mathrm{d} x}=\frac{2y}{x}+x^2 e^x"
text = L"y=x^2(e^x - e)"
show_result(f1, f2, f3, xspan, y0, iternums, true, false, title, text) # show=true, dense=true

# ╔═╡ b17d6cbc-3e64-4a0a-aa35-948158dd1c3d
iternums = [5, 10, 20]

f1(y, p, x) = (y^2 + y) / x
xspan = (1.0, 3.0)
y0 = -2.0
f2(x, y) = (y^2 + y) / x
f3(x) = 2x / (1 - 2x)
title = L"Problem\ 2.2:\frac{\mathrm{d} y}{\mathrm{d} x}=\frac{1}{x}(y^2+y)"
text = L"y=\frac{2x}{1-2x}"
show_result(f1, f2, f3, xspan, y0, iternums, true, false, title, text) # show=true, dense=true

# ╔═╡ 9b374966-a212-468e-949b-4950f04df864
iternums = [5, 10, 20]

f1(y, p, x) = -20(y - x^2) + 2x
xspan = (0.0, 1.0)
y0 = 1 / 3
f2(x, y) = -20(y - x^2) + 2x
f3(x) = x^2 + 1 / 3 * exp(-20x)
title = L"Problem\ 3.1: \frac{\mathrm{d} y}{\mathrm{d} x}=-20(y-x^2)+2x"
text = L"y=x^2+\frac{1}{3}e^{-20x}"
show_result(f1, f2, f3, xspan, y0, iternums, true, false, title, text) # show=true, dense=true

# ╔═╡ ba6e68e2-3fa7-419d-9442-04ac13e2ea7f
iternums = [5, 10, 20]

f1(y, p, x) = -20y + 20sin(x) + cos(x)
xspan = (0.0, 1.0)
y0 = 1.0
f2(x, y) = -20y + 20sin(x) + cos(x)
f3(x) = exp(-20x) + sin(x)
title = L"Problem\ 3.2: \frac{\mathrm{d} y}{\mathrm{d} x}=-20y+20\sin(x)+\cos(x)"
text =  L"y=e^{-20x}+\sin(x)"
show_result(f1, f2, f3, xspan, y0, iternums, true, false, title ,text) # show=true, dense=true

# ╔═╡ 3e5facb5-6e75-41e6-b8de-2f40c4bf2835
iternums = [5, 10, 20]

f1(y, p, x) = -20(y - exp(x)sin(x)) + exp(x) * (sin(x) + cos(x))
xspan = (0.0, 1.0)
y0 = 0.0
f2(x, y) = -20(y - exp(x)sin(x)) + exp(x) * (sin(x) + cos(x))
f3(x) = exp(x) * sin(x)
title = L"Problem\ 3.3: \frac{\mathrm{d} y}{\mathrm{d} x}=-20(y-e^x \sin(x))+e^x (\sin(x) + \cos(x))"
text = L"y=e^x \sin(x)"
show_result(f1, f2, f3, xspan, y0, iternums, true, false, title, text) # show=true, dense=true

# ╔═╡ 3fd2f570-3323-4260-bff5-048e26478c5d
md"""
### 思考题
"""

# ╔═╡ feecb444-040d-46c0-8e2a-ce290e8887c4
md"""

1. 对实验 1，数值解和解析解相同吗？为什么？试加以说明。
   
   对于问题1.1，数值解和解析解是相同的，因为本题的解是线性函数，能够通过所得数值解的两个点确定直线的方程，即等价于得到了解析解。
   
   本例中，待求解微分方程为$\frac{\mathrm{d} y}{\mathrm{d} x} = x + y$，解为$y = -x - 1$，而`rungekutta()`函数求解的任意两点（如`(0.2,-1.2)`, `(1.0,-2.0)`）所决定的直线方程即为$y = -x - 1$。
    
   而对于问题1.2，虽然数值解和解析解之间差异已经极小（绝对误差在1e-7~1e-5数量级，仅仅对比相同x所在的y取值，如下表所示），但对于非线性函数$y=\frac{1}{1+x}$，在未知函数解析式类型的情况下，是几乎不可能仅仅通过数值解所求得的点，来推断准确的函数解析式的，此时不能认为所求得的数值解就是解析解。
   ```
   ┌────────────────┬────────────────┬────────────────┬────────────────┬────────────────┐
   │     Test x     │     True y     │ 5-Iter Pred y  │ 10-Iter Pred y │ 20-Iter Pred y │
   ├────────────────┼────────────────┼────────────────┼────────────────┼────────────────┤
   │     0.20000000 │     0.83333333 │     0.83333904 │     0.83333373 │     0.83333336 │
   │     0.40000000 │     0.71428571 │     0.71429213 │     0.71428615 │     0.71428574 │
   │     0.60000000 │     0.62500000 │     0.62500589 │     0.62500040 │     0.62500003 │
   │     0.80000000 │     0.55555556 │     0.55556069 │     0.55555590 │     0.55555558 │
   │     1.00000000 │     0.50000000 │     0.50000441 │     0.50000030 │     0.50000002 │
   └────────────────┴────────────────┴────────────────┴────────────────┴────────────────┘
   ```

2. 对实验 2，N 越大越精确吗？试加以说明。
   
   虽然确实N越大越精确，但从本例实验的结果来看，因为当n=5的时候已经获得足够精确的数值解了，再增大n的值只是增加了计算量，却不能再明显提高结果的精度，此时我们不能一味的增大N，而要根据所需要达到的精度要求及时终止计算。
   
   本例中，$y=x^2(e^x - e)$，在迭代次数从5增加到20的时候，数值上的精度只增加了2位，继续增大n对于所求数值解精度改变很小，很难继续使用Runge-Kutta方法继续进行求解，并且这样的计算资源成本是不可忽略的。
   ```
   ┌────────────────┬────────────────┬────────────────┬────────────────┬────────────────┐
   │     Test x     │     True y     │ 5-Iter Pred y  │ 10-Iter Pred y │ 20-Iter Pred y │
   ├────────────────┼────────────────┼────────────────┼────────────────┼────────────────┤
   │     1.40000000 │     2.62035955 │     2.61394279 │     2.61974052 │     2.62031131 │
   │     1.80000000 │    10.79362466 │    10.77631317 │    10.79201760 │    10.79350178 │
   │     2.20000000 │    30.52458129 │    30.49165420 │    30.52159814 │    30.52435589 │
   │     2.60000000 │    72.63928396 │    72.58559861 │    72.63450354 │    72.63892578 │
   │     3.00000000 │   156.30529585 │   156.22519828 │   156.29825744 │   156.30477188 │
   └────────────────┴────────────────┴────────────────┴────────────────┴────────────────┘
   ```
   
3. 对实验 3，N 较小会出现什么现象？试加以说明
 
   当n较小的时候所得数值解和正确结果相差较大，结果失真，说明在一定条件下确实需要更大的n来更好的获得数值解。而具体这个n的大小如何选取则取决于待求解微分方程性质，这里应该涉及到更深入的课程或者研究。

   对本例而言，从下表以及所绘制的图像都很容易能看到，当n较小的时候会导致求得数值解偏差极大，甚至于几乎就完全是错误的（大约与正确结果相差1e3的量级），所以选择充分大的n，并设置结果收敛的措施，才能确保最终可以得到精度合适的数值解的同时不会造成太大的计算资源浪费。
   
   下表为了便于对齐，略去了多余的x数据，方程的解析解为$y=e^{-20x}+\sin(x)$，数值解如下所示：
   ```
   ┌────────────────┬────────────────┬────────────────┬────────────────┬────────────────┐
   │       x        │     True y     │ 5-Iter Pred y  │ 10-Iter Pred y │ 20-Iter Pred y │
   ├────────────────┼────────────────┼────────────────┼────────────────┼────────────────┤
   │     0.20000000 │     0.04610521 │     1.76000000 │     0.07925926 │     0.04667348 │
   │     0.40000000 │     0.16011182 │     8.81333333 │     0.16658436 │     0.16021366 │
   │     0.60000000 │     0.36000205 │    43.68000000 │     0.36295382 │     0.36008591 │
   │     0.80000000 │     0.64000004 │   217.29333333 │     0.64255042 │     0.64008338 │
   │     1.00000000 │     1.00000000 │  1084.32000000 │     1.00250560 │     1.00008333 │
   └────────────────┴────────────────┴────────────────┴────────────────┴────────────────┘
   ```
   以下为方程$\frac{\mathrm{d} y}{\mathrm{d} x}=-20(y-e^x \sin(x))+e^x (\sin(x) + \cos(x))$的部分数值解表格，为便于集中观察而总结如下，解析解为$y=e^x \sin(x)$，

   ```
   ┌────────────────┬────────────────┬────────────────┬────────────────┬────────────────┐
   │       x        │     True y     │ 5-Iter Pred y  │ 10-Iter Pred y │ 20-Iter Pred y │
   ├────────────────┼────────────────┼────────────────┼────────────────┼────────────────┤
   │     0.20000000 │     0.24265527 │     0.29864621 │     0.24511651 │     0.24274900 │
   │     0.40000000 │     0.58094390 │     0.92721987 │     0.58409696 │     0.58105449 │
   │     0.60000000 │     1.02884567 │     2.83547734 │     1.03241831 │     1.02896834 │
   │     0.80000000 │     1.59650534 │    10.71088533 │     1.60032101 │     1.59663402 │
   │     1.00000000 │     2.28735529 │    47.94144638 │     2.29115692 │     2.28748035 │
   └────────────────┴────────────────┴────────────────┴────────────────┴────────────────┘
   ```
   
"""

# ╔═╡ f80f51b3-6a23-477a-86bd-639693748b63
display(md"""
### 参考资料

1. julia ordinary differential equations tutorial https://diffeq.sciml.ai/stable/tutorials/ode_example/

2. intro to solving differential equations in julia https://www.youtube.com/watch?v=KPEqYtEd-zY

3. julia ode solver type: Runge-Kutta https://diffeq.sciml.ai/stable/solvers/ode_solve/#Explicit-Runge-Kutta-Methods

4. julia ode problem type https://diffeq.sciml.ai/stable/types/ode_types/#ode_prob

5. julia ode speed up perf https://diffeq.sciml.ai/stable/features/performance_overloads/#performance_overloads

6. julia ode common solver option https://diffeq.sciml.ai/stable/basics/common_solver_opts/#solver_options

7. 《计算方法实验指导》实验题目 3 四阶龙格—库塔(Runge—Kutta)方法
""")
