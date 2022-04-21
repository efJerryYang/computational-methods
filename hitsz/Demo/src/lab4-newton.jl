
using Markdown

# ╔═╡ 08ed10aa-cd4a-4732-b8b4-aebe9b0a1226
display(md"""
## 实验题目4 牛顿(Newton)迭代法
""")

# ╔═╡ 96cdaee8-8eef-4485-a2c2-da2d42a1a4e0
using Printf
# using Plots
using NLsolve
using SymPy
using Roots
using PrettyTables
using LaTeXStrings

# ╔═╡ 898e89eb-d53e-4b0d-846f-9210703a02df
md"""
#### newton

本部分为牛顿法的代码实现，使用重根的迭代法作为外层代码，牛顿法实现为重根修正的特殊形式，即``\lambda=1``
"""

# ╔═╡ 2a245c46-f12d-4950-9738-fc18b37f3cd3
# multi-root newton method
function newton(f::Function, df::Function, ϵ1, ϵ2, N, x0, λ)
    n = 1
    while n <= N
        F = f(x0)
        DF = df(x0)
        if abs(F) < ϵ1
            return n - 1, x0
        end
        if abs(DF) < ϵ2
            @printf("Reach a critical point!\n")
            return
        end
        x1 = x0 - λ * F / DF
        Tol = abs(x1 - x0)
        if Tol < ϵ1
            return n - 1, x1
        end
        n = n + 1
        x0 = x1
    end
    @printf("Fail to converge within %d iterations!\n", N)
end
# newton method
function newton(f::Function, df::Function, ϵ1, ϵ2, N, x0)
    newton(f, df, ϵ1, ϵ2, N, x0, 1.0)
end

# ╔═╡ 1b10f5fb-e597-4cd8-a2e8-03f974053d65
function get_func_diff(symf::Sym)
    @syms x
    symdf = diff(symf)
    f = lambdify(symf)
    df = lambdify(symdf)
    f, df
end
function get_func_diff(f::Function)
    @syms x
    symf = f(x)
    symdf = diff(symf)
    df = lambdify(symdf)
    f, df
end
function redefine_func(f::Function, df::Function)
    function f!(r, x)
        r .= f.(x)
    end
    function j!(J, x)
        (s1, s2) = size(J)
        J .= zeros(s1, s1)
        for i in 1:s1
            J[i, i] = df(x[i])
        end
    end
    f!, j!
end

function collect_data(f::Function, df::Function, f!::Function, j!::Function,
    x0, ϵ1, ϵ2, N, λ)

    t1 = @elapsed rs = fzeros(f, x0 - 0.5, x0 + 0.5)
    r1 = rs[1]
    i1 = NaN

    t2 = @elapsed sol = nlsolve(f!, j!, [x0]; method=:newton, ftol=ϵ1, iterations=N)
    r2 = sol.zero[1]
    # println(sol)
    i2 = sol.iterations[1]

    t3 = @elapsed i3, r3 = newton(f, df, ϵ1, ϵ2, N, x0)
    if λ > 1
        t4 = @elapsed i4, r4 = newton(f, df, ϵ1, ϵ2, N, x0, λ)
        xs = [r1, r2, r3, r4]
        ts = [t1, t2, t3, t4]
        is = [i1, i2, i3, i4]
        method_name = ["lib Roots", "lib NLsolve", "my newton", "my newton (λ=$λ)"]
        data = [method_name xs f.(xs) ts is]
        return data
    else
        xs = [r1, r2, r3]
        ts = [t1, t2, t3]
        is = [i1, i2, i3]
        method_name = ["lib Roots", "lib NLsolve", "my newton"]
        data = [method_name xs f.(xs) ts is]
        return data
    end
end

# ╔═╡ 3c3d940c-2869-467c-8230-68336edf8f3b
function show_result(f::Function, x0, ϵ1, ϵ2, N, title)
    f, df = get_func_diff(f)
    f!, j! = redefine_func(f, df)

    header = (["Method", "x", "f(x)", "Time Cost (s)", "Iter"])
    data = collect_data(f, df, f!, j!, x0, ϵ1, ϵ2, N, 1)
    pretty_table(
        data;
        alignment=[:c, :c, :c, :c, :c],
        header=header,
        header_crayon=crayon"bold"
    )

end
function show_result(f::Function, x0, λ, ϵ1, ϵ2, N, title)
    f, df = get_func_diff(f)
    f!, j! = redefine_func(f, df)

    header = (["Method", "x", "f(x)", "Time Cost (s)", "Iter"])
    data = collect_data(f, df, f!, j!, x0, ϵ1, ϵ2, N, λ)
    pretty_table(
        data;
        alignment=[:c, :c, :c, :c, :c],
        header=header,
        header_crayon=crayon"bold"
    )
end


# ╔═╡ 740751b5-cd64-430e-a8d6-5acfea82056b
f(x) = cos(x) - x
ϵ1 = 1e-6
ϵ2 = 1e-4
N = 10
x0 = pi / 4
title = L"Problem\ 1.1:\ f(x)=\cos x -x=0,~~~~x_0=\frac{\pi}{4}\approx 0.78539816"
println()
println(title)
show_result(f, x0, ϵ1, ϵ2, N, title)

# ╔═╡ d0e762e5-19c9-4a5c-84f0-076cd746d9b7
f(x) = exp(-x) - sin(x)
ϵ1 = 1e-6
ϵ2 = 1e-4
N = 10
x0 = 0.6
title = L"Problem\ 1.2:\ f(x)=e^{-x}-\sin x=0,~~~~x_0=0.6"
println()
println(title)
show_result(f, x0, ϵ1, ϵ2, N, title)


# ╔═╡ 51539216-4021-401b-9b29-5b207df57f07
f(x) = x - exp(-x)
ϵ1 = 1e-6
ϵ2 = 1e-4
N = 10
x0 = 0.5
title = L"Problem\ 2.1:\ f(x)=x-e^{-x}=0,~~~~x_0=0.5"
println()
println(title)
show_result(f, x0, ϵ1, ϵ2, N, title)


# ╔═╡ a9d91d31-4460-480f-9f2c-387c73c95e53
f(x) = x^2 - 2x * exp(-x) + exp(-2x)
ϵ1 = 1e-6
ϵ2 = 1e-4
N = 20
x0 = 0.5
λ = 2
title = L"Problem\ 2.2:\ f(x)=x^2 -2x e^{-x} + e^{-2x}=0,~~~~ x_0=0.5"
println()
println(title)
show_result(f, x0, λ, ϵ1, ϵ2, N, title)

# ╔═╡ bf8ffa46-acbe-4464-9cff-e39f661c52fc
md"""
### 思考题
"""

# ╔═╡ 0dbdadd9-e908-476c-a0d2-4cb028a98563
md"""
1. 对问题 1 确定初值的原则是什么？实际计算中应如何解决？

   选择一个有根区间，本例中容易得到
   
   ``f(x)=\cos(x)-x``，有``f(0)=1>0,f(\frac{\pi}{2})=-\frac{\pi}{2}<0``，取区间中点即``x=\frac{\pi}{4}``为初值
   
   ``f(x)=e^{-x} - \sin(x)``，有``f(0)=1>0,f(1.2)\approx-0.631<0``，同样取区间中点即``x=0.6``为初值
   
   实际计算中，根据其他算法求出多个精度较粗的有根区间，然后使用牛顿法逼近获得较为精确的数值解。
   
   通常，对于我们而言，可能在给定区间直接作出函数图像是最简单最可行的方式<sup>8,9</sup>。

   不过也存在一些数值的算法，可以计算给定区间函数根的情况，如`Julia`的`Roots.jl`库中的`fzeros()`函数，是本实验中尝试的所有方法中效率最高的，精度最高和计算耗时最短。

   除此以外，在我所查找的资料中提及，对于更一般的情形，试图通过程序自动化来计算函数的根的话，情况会变得更加复杂，涉及到多个领域的研究<sup>10,11</sup>。
   
2. 对问题 2 如何解释在计算中出现的现象？试加以说明

   本例中，
   
   ``(1)f_1(x) = x - e^{-x}=0``

   ```
   ┌─────────────┬──────────┬────────────┬───────────────┬──────┐
   │   Method    │    x     │    f(x)    │ Time Cost (s) │ Iter │
   ├─────────────┼──────────┼────────────┼───────────────┼──────┤
   │  lib Roots  │ 0.567143 │    0.0     │   0.0038645   │ NaN  │
   │ lib NLsolve │ 0.567143 │ -1.9648e-7 │   0.0644682   │ 2.0  │
   │  my newton  │ 0.567143 │ -1.9648e-7 │   0.0200168   │ 2.0  │
   └─────────────┴──────────┴────────────┴───────────────┴──────┘
   ```
   ``(2)f_2(x) = x^2 - 2xe^{-x} + e^{-2x}=(x - e^{-x})^{2}={f_1}^2(x)=0``
      
   ```
   ┌─────────────────┬──────────┬─────────────┬───────────────┬──────┐
   │     Method      │    x     │    f(x)     │ Time Cost (s) │ Iter │
   ├─────────────────┼──────────┼─────────────┼───────────────┼──────┤
   │    lib Roots    │ 0.567143 │     0.0     │   0.0050046   │ NaN  │
   │   lib NLsolve   │ 0.566606 │ 7.09902e-7  │   0.0675075   │ 7.0  │
   │    my newton    │ 0.566606 │ 7.09902e-7  │   0.0201382   │ 7.0  │
   │ my newton (λ=2) │ 0.567143 │ 3.86358e-14 │   0.0156408   │ 2.0  │
   └─────────────────┴──────────┴─────────────┴───────────────┴──────┘
   ```
   
   显然，方程(2)在方程(1)的根位置有重根，可以看到直接应用牛顿迭代法计算轮数为7轮，耗时通常情况下稍微增加，但由于当前实例计算简单、精度要求低，耗时变化不明显，不过能看到在给定精度要求情况下所得精度低于无重根牛顿迭代法。

   由理论课知识可知，当存在重根时牛顿迭代法的收敛速度为线性收敛。在后续使用修正的牛顿法可以使收敛速度重新达到平方收敛，耗时几乎与无重根时一致，迭代次数和精度也相同。

   除此以外，在本实验中，直接手写的牛顿法`my newton`运行效率意外的高于库`NLsolve.jl`指定`nlsolve()`参数`method=:newton`后的求解效率，耗时更短，这是让人十分意外的。当然，事实上在同一数量级时间的差异并不大，而`NLsolve.jl`库通用于求解非线性系统，仍然是指定了使用牛顿法来求解实际问题时的优选。而对于一般求解数值根，在给定区间时，使用求解方式得到优化的`Roots.jl`库中的`fzeros()`函数应当为效率最高的选择。

3. 略
"""

# ╔═╡ 0efe3358-381d-4ae6-bb3d-26c61a5359d1
display(md"""
### 参考资料

1. julia covert sym to func https://stackoverflow.com/questions/27357579/julia-how-do-i-convert-a-symbolic-expression-to-a-function
2. julia SymPy Tutorial https://docs.juliahub.com/SymPy/KzewI/1.0.28/Tutorial/basic_operations/#lambdify-1
3. julia SymPy Tutorial https://docs.juliahub.com/SymPy/KzewI/1.0.28/Tutorial/solvers/
4. julia time elapse https://discourse.julialang.org/t/difference-between-tic-toc-time-or-elapse-in-julia/1177/2
5. julia NLsolve github https://github.com/JuliaNLSolvers/NLsolve.jl
6. 《数值分析原理》吴勃英 18-19,27-29
7. 《计算方法实验指导》实验题目 4 牛顿(Newton)迭代法
8. how to find initial guess https://computingskillset.com/solving-equations/how-to-find-the-initial-guess-in-newtons-method/ 
9. how to choose starting point https://math.stackexchange.com/questions/743373/how-to-choose-the-starting-point-in-newtons-method
10. wikipedia newton fractal https://en.wikipedia.org/wiki/Newton_fractal
11. How to find all roots of complex polynomials by Newton’s method https://link.springer.com/article/10.1007/s002220100149
""")
