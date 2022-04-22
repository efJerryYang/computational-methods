using Markdown

display(md"""
## 实验题目2 龙贝格(Romberg)积分法
""")
using Printf
function romberg(f::Function, xlim, n, ϵ)
    a, b = xlim
    h = b - a
    T = zeros(n, n)
    T[1, 1] = 1 / 2 * h * (f(a) + f(b))
    for i = 1:n
        tmpsum = 0
        jmax = 2^(i - 1)
        for j = 1:jmax
            tmpsum += f(a + (j - 1 / 2) * h)
        end
        T[i+1, 1] = 1 / 2 * T[i, 1] + 1 / 2 * h * tmpsum

        for m = 1:i
            T[i+1, m+1] = (4^m * T[i+1, m] - T[i, m]) / (4^m - 1)
        end
        for m = 1:i
            @printf("%12.9f\t", T[i, m])
        end
        @printf("\n")
        if i > 1 && abs(T[i+1, i+1] - T[i, i]) < ϵ
            @printf("Accuracy requirement satisfied.\n\n")
            break
        end
        h /= 2
    end
end

# ╔═╡ 5090d873-7ff5-4b0e-8476-b6d788c62462
iter_num = 30

f(x) = x^2 * exp(x)
ϵ = 1e-6
xlim = 0, 1
println()
println("Problem 1.1: f(x) = x^2 * exp(x)")
romberg(f, xlim, iter_num, ϵ)

f(x) = exp(x)sin(x)
ϵ = 1e-6
xlim = 1, 3
println("Problem 1.2: f(x) = exp(x)sin(x)")
romberg(f, xlim, iter_num, ϵ)


f(x) = 4 / (1 + x^2)
ϵ = 1e-6
xlim = 0, 1
println("Problem 1.3: f(x) = 4 / (1 + x^2)")
romberg(f, xlim, iter_num, ϵ)

f(x) = 1 / (x + 1)
ϵ = 1e-6
xlim = 0, 1
println("Problem 1.3: f(x) = 1 / (x + 1)")
romberg(f, xlim, iter_num, ϵ)

# ╔═╡ d187517e-5f21-4b4c-ae81-d519777b145d
md"""
### 思考题
"""

# ╔═╡ 18463200-a64a-400a-a046-b215fd059e7a
md"""
1. 略

2. 在实验 1 中二分次数和精度的关系如何？

   二分次数越多所求的精度越高，通常预设较大的二分次数来确保计算结果有足够的精度，同时也设定早停需要满足的精度要求，避免达到所需精度之后继续计算导致增加的运算量

3. 略

4. 略
"""

# ╔═╡ 40e75753-9cdc-44cc-a8e3-b8953d56fb37
display(md"""
### 参考资料

1. julia 数值积分 https://blog.csdn.net/m0_37816922/article/details/103475445
2. Romberg Integration-Numerical Analysis http://homepages.math.uic.edu/~jan/mcs471/romberg.pdf
3. 《数值分析》吴勃英 196-199

""")
