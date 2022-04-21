using Markdown

display(md"""
## 实验题目1 拉格朗日(Lagrange)插值
""")


using Printf
using LinearAlgebra
using LaTeXStrings
using PrettyTables
function lagrange(xs, fxs, x::Number)
    num = size(xs, 1)
    y, i = 0.0, 1
    while i <= num
        li = 1.0
        for j in 1:num
            if j == i
                continue
            end
            li *= (x - xs[j]) / (xs[i] - xs[j])
        end
        y += li * fxs[i]
        i += 1
    end
    x, y
end
function lagrange(xs, fxs, x::Vector)
    num = size(xs, 1)
    y, i = zeros(size(x, 1)), 1
    while i <= num
        li = fill(1.0,size(x,1))
        for j in 1:num
            if j == i
                continue
            end
            li = li .* (x .- xs[j]) / (xs[i] - xs[j])
        end
        y = y + li .* fxs[i]
        i += 1
    end
    x, y
end

function show_result(f::Function, split_nums::Vector, test_x::Vector, xlim::Vector, ylim::Vector, prefix, text)
    for n in split_nums
        # initialization
        x_min, x_max = xlim
        # x_range = x_min-0.2:0.02:x_max+0.2
        xs = x_min:(x_max-x_min)/n:x_max
        ys = f.(xs)

        # plot(x_range, f.(x_range), label="f(x)")  # plot f(x)
        # plot!(legend=:outertopright, title=prefix * " $n-Order Interpolation")

        # series_x = Vector(x_range)
        # _, series_y = lagrange(xs, ys, series_x)  # compute the interpolation function points
        # plot!(series_x, series_y, color=:violet, label="p(x)")  # add p(x) function curve

        # plot!(ylim=ylim, yflip=false)  # add ylim
        # add sample for lagrange interpolation
        # plot!(xs, ys, seriestype=:scatter, markersize=3, msw=1, color=:deepskyblue, label="sample")  

        test_y = f.(test_x)
        # add test x & y, plot true points
        # p = plot!(test_x, test_y, seriestype=:scatter, markersize=3, msw=1, color=:blue, label="true")  
        _, pred_y = lagrange(xs, ys, test_x)
        println()
        println(prefix * " $n-Order Interpolation:")
        
        data = [test_x test_y pred_y]
        header = (["Test x", "Test y", "Pred y"])
        pretty_table(
            data;
            alignment=[:c, :c, :c],
            header=header,
            formatters=ft_printf("%11.6f"))
    end
end

function show_result(f::Function, split_nums::Nothing, split_xs::Vector, test_x, xlim, ylim, prefix, text, comment)
    xs = split_xs
    ys = f.(xs)
    test_y = f.(test_x)
    _, pred_y = lagrange(xs, ys, test_x)

    println()
    data = [test_x test_y pred_y]
    header = (["Test x", "Test y", "Pred y"])
    pretty_table(
        data;
        alignment=[:c,:c,:c],
        header=header,
        formatters=ft_printf("%11.6f"))
end


# ╔═╡ 0c406504-4d55-4716-8e00-85b978ff0e6a
f(x) = 1 / (1 + x^2)
split_nums = [5, 10, 20]
test_x = [0.75, 1.75, 2.75, 3.75, 4.75]
xlim = [-5, 5]
ylim = [-1, 2]
println()
println("f(x) = 1 / (1 + x^2)")
prefix = "Problem 1.1 "
text = L"f(x)=\frac{1}{1+x^2}"
show_result(f, split_nums, test_x, xlim, ylim,prefix,text)

# ╔═╡ f5d1cdcd-6b2c-48ed-8a53-5083fd777c9b
f(x) = exp(x)
split_nums = [5, 10, 20]
test_x = [-0.95, -0.05, 0.05, 0.95]
xlim = [-1, 1]
# ylim = [-1, 10]  # the good-looking ylim is defined manually
ylim = []
println()
println("f(x) = exp(x)")
prefix = "Problem 1.2 "
text = L"f(x)=e^x"
show_result(f, split_nums, test_x, xlim, ylim, prefix, text)

# ╔═╡ c7bf369b-29e3-4c57-93b8-c70d54d997a2
f(x) = 1 / (1 + x^2)
split_nums = [5, 10, 20]
test_x = [-0.95, -0.05, 0.05, 0.95]
xlim = [-1, 1]
# ylim = [-1, 2]
ylim = []
println()
println("f(x) = 1 / (1 + x^2)")
prefix = "Problem 2.1 "
text = L"f(x) = \frac{1}{1+x^2}"
show_result(f, split_nums, test_x, xlim, ylim, prefix, text)

# ╔═╡ 11470f6a-5dc5-4314-9628-439d2eac0f6b
f(x) = exp(x)
split_nums = [5, 10, 20]
test_x = [0.75, 1.75, 2.75, 3.75, 4.75]
xlim = [-5, 5]
# ylim = [-1, 10]  # the good-looking ylim is defined manually
ylim = []
println()
println("f(x) = exp(x)")
prefix = "Problem 2.2 "
text = L"f(x) = e^x"
show_result(f, split_nums, test_x, xlim, ylim, prefix, text)

# ╔═╡ 8fa9e9fb-577f-4a26-8639-b34f3eaeaa05
f(x) = sqrt(x)
split_xs = [1, 4, 9]
test_x = [5, 50, 115, 185]
xlim = [0, 200]
# ylim = [-1, 2]
ylim = []
println()
print("Problem 4.1  f(x) = sqrt(x)")
prefix = "Problem 4.1 "
text = L"f(x) = \sqrt{x}"
comment = "This is a test commen"
comment = L"last~3~points~are~extrapolation~results"
show_result(f, nothing, split_xs, test_x, xlim, ylim, prefix, text, comment)

f(x) = sqrt(x)
split_xs = [36, 49, 64]
test_x = [5, 50, 115, 185]
xlim = [0, 200]
# ylim = [-1, 2]
ylim = []
println()
print("Problem 4.2  f(x) = sqrt(x)")
prefix = "Problem 4.2 "
text = L"f(x) = \sqrt{x}"
comment = L"only~the~second~one~is~interpolation~result"
show_result(f, nothing, split_xs, test_x, xlim, ylim, prefix, text, comment)

f(x) = sqrt(x)
split_xs = [100, 121, 144]
test_x = [5, 50, 115, 185]
xlim = [0, 250]
# ylim = [-1, 2]
ylim = []
println()
print("Problem 4.3  f(x) = sqrt(x)")
prefix = "Problem 4.3 "
text = L"f(x) = \sqrt{x}"
comment = L"only~the~third~one~is~interpolation~result"
show_result(f, nothing, split_xs, test_x, xlim, ylim, prefix, text,comment)

f(x) = sqrt(x)
split_xs = [169, 196, 225]
test_x = [5, 50, 115, 185]
xlim = [0, 250]
# ylim = [-1, 2]
ylim = []
println()
print("Problem 4.4  f(x) = sqrt(x)")
prefix = "Problem 4.4 "
text = L"f(x) = \sqrt{x}"
comment = L"first~3~points~are~extrapolation~results"
show_result(f, nothing, split_xs, test_x, xlim, ylim, prefix, text, comment)

# ╔═╡ 50494629-ec79-4842-8441-a2f42e98d4f4
md"""
### 思考题
"""

# ╔═╡ 05069d8d-651d-4cc3-bc51-594ed8222de9
md"""
1. 对实验 1 存在的问题，应如何解决？

   当插值多项式次数过高的时候会出现Runge现象，插值多项式在距离已知点位置较远处会剧烈震荡，越靠近端点，逼近的效果越差，这表明了节点的密集不一定能保证在两节点间插值函数逼近程度的上升。

   对函数``f(x)=\frac{1}{1+x^2}``进行5阶、20阶插值的相对误差的变化，如下图所示：

   很明显看到，误差在距离已知点之外的急剧增长，相对误差高达10^4^量级，结果已经严重失真。目前而言，这一问题的解决方案主要有两类：

   一个是从插值函数的二阶导数剧烈变化出发，修改插值条件对插值函数的二阶导数进行限制，如使用Hermite型插值；

   另一种是将长区间划分为若干个小区间，在每一个小区间上分别做低次插值来避免Runge现象，逼近效果要比在整个区间上用高阶光滑差值效果更好，即使用分段插值和样条插值。

2. 对实验 2 存在的问题的回答，试加以说明

   首先不一定，从精度上考虑虽然有一定的合理性，但插值节点过于密集时，一方面计算量增大却没提高对于精度计算的收益，另一方面区间缩短、节点增加并不能保证两节点间能很好的逼近函数，反而有可能出现Runge现象。但合理的对区间长度进行选择，同时采用低次插值来避免Runge现象，能够得到较好的拟合效果。

   不过，实例中对于函数``f(x)=\frac{1}{1+x^2}``，较短区间的插值效果比长区间插值更好，而且短区间插值甚至缓解了Runge现象对于精度的影响，虽然相对而言边界点仍有较大的误差，但其相对误差的绝对大小也明显低于低阶插值的结果，如下图所示；而函数``f(x)=e^x``无论是长区间还是短区间插值，都能得到相对较好的拟合效果，但短区间插值相对误差更低，此处不再赘述。

3. 略

4. 如何理解插值问题中的内插和外推？

   通常我们认为对于连续函数内插的可靠程度高于外推，因为对于未知的连续函数而言，我们无法预知任何在已知点信息之外的有关函数的信息，无法简单通过多项式插值来对函数趋势进行判断。
   
   外推等价于根据已知点预测完全未知点的函数值，但我们所得的插值多项式不含有任何有关待拟合函数的已知点外的信息，我们没有理由认为实际函数的变化必须符合多项式函数的变化，根据多项式函数的特性进行外推是不合理的，故这里我们可以简单的认为即时出现外推的结果相比于内插更好的情况，也很大程度上是函数本身特性导致的巧合，当然，或许也可以根据函数的性质来决定函数是否适合进行多项式差值的外推。
   
   认为内插比外推可靠的原因是，在内插的区间我们认为连续函数在相邻点间的变化不会过于剧烈，从而可以简单的认为内插更加可靠，而外推时我们没有任何而先验知识可以用于对完全未知的区间函数值进行推断。
   
   从实验结果来看，第一个实例体现的是外推的严重错误（如下图所示），接下来的几个实例中外推所得误差和内插结果相差不明显，但在即便如此也只是所选区间拟合的巧合，因为第一个实例中的极端例子就说明了，不加考量的直接外推会造成的灾难性后果。
"""

# ╔═╡ e7fc0eac-142c-469f-9829-226b701b50ac
display(md"""
### 参考资料

1. julia plot xlim and ylim https://stackoverflow.com/questions/53230969/how-to-scale-a-plot-in-julia-using-plots-jl
2. julia fill array with specific value https://www.geeksforgeeks.org/fill-an-array-with-specific-values-in-julia-array-fill-method/
3. julia prettytabels https://ronisbr.github.io/PrettyTables.jl/stable/
4. 《数值分析原理》吴勃英 105-106,123
5. markdown多张图片并排显示 https://www.cnblogs.com/jaycethanks/p/12201959.html
6. markdown图片大小设定 https://www.cnblogs.com/jaycethanks/p/12202169.html
7. 《计算方法实验指导》实验题目 1 拉格朗日(Lagrange)插值
""")
