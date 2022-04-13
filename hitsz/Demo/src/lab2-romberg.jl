### A Pluto.jl notebook ###
# v0.10.0

using Markdown

# ╔═╡ 1198114b-4ee5-479d-a539-4770886bbaa4
display(md"""
## 实验题目2 龙贝格(Romberg)积分法
""")

# ╔═╡ f3390402-8050-4b60-b93e-31f42897a94d
md"""
### 实验简介

本实验为龙贝格积分法，需要完成龙贝格积分代码的编写，并对实验题目进行求解。

本次实验过程中，主要是对龙贝格积分法的算法流程进行调试，从数学原理直接构造算法流程用于计算，充分体会了`Julia`编程语言的流畅性，感受到了编写代码时同`Python`一样自如，却能拥有和`C`相比的循环结构。

实验的目的为使用龙贝格积分法计算定积分，并输出T数表。

该实验报告主要分7个部分，大纲罗列如下：

- 实验简介：即本部分分所有内容
- 数学原理：即龙贝格积分法的数学公式，用于改写为算法流程
- 代码实现：使用`Julia`编程语言，根据数学原理，编写实验代码
- 测试代码：教材上的例题用于对程序结果进行初步的检验
- 实验题目：实验指导书中所要求的题目
  - 问题1：直接使用龙贝格积分法计算定积分的值
- 思考题：本部分为实验指导书中所要求完成的思考题解答
- 参考资料：本部分为本次实验过程中查阅的参考资料
"""

# ╔═╡ 423c9cdd-ea9f-49ce-9f84-252c63dfa582
md"""
### 数学原理

教材中给出的计算公式如下
```math
\begin{cases}
	T_{0,0}&=\frac{b-a}{2}\left[ f\left( a \right) +f\left( b \right) \right] ,\\
	T_{0,i}&=\frac{1}{2}T_{0,i-1}+\frac{1}{2}\frac{b-a}{2^{i-1}}\sum_{j=1}^{2^{i-1}}{f\left[ a+\left( j-\frac{1}{2} \right) \cdot \frac{b-a}{2^{i-1}} \right]}, i=1,2,3,...,\\
	T_{m,k}&=\frac{4^mT_{m-1,k+1}-T_{m-1,k}}{4^m-1},m=1,2,...,i; k=i-m.\\
\end{cases}
```
因`Julia`语言数组类下标的起点为1，同时实验指导书所给T数表为下三角形，故将原公式改写如下
```math
\begin{aligned}
\begin{cases}
	T_{1,1}&=\frac{b-a}{2}\left[ f\left( a \right) +f\left( b \right) \right] ,\\
	T_{i+1,1}&=\frac{1}{2}T_{i,1}+\frac{1}{2}\frac{b-a}{2^{i-1}}\sum_{j=1}^{2^{i-1}}{f\left[ a+\left( j-\frac{1}{2} \right) \cdot \frac{b-a}{2^{i-1}} \right]}, i=1,2,3,...,n\\
	T_{i+1,m+1}&=\frac{4^mT_{i+1,m}-T_{i,m}}{4^m-1},m=1,2,...,i.\\
\end{cases}
\end{aligned}
```
随后可对照公式完成代码的编写
"""

# ╔═╡ 654ffab7-7e9d-4066-97f5-edc69b1b619f
md"""
### 代码实现

使用`Julia`编程语言，根据上述数学原理，编写`romberg`积分法实验代码。

以下部分为`romberg()`函数定义：
"""

# ╔═╡ 83f83660-08c9-4e28-aedb-81c283097079
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

# ╔═╡ 506efb71-374f-42f8-87e4-2b1f5ef50c06
md"""
### 测试代码

本部分使用教材上的例题用于对程序结果进行初步的检验，计算结果和教材给出数表类似，可以认为测试通过。
"""

# ╔═╡ 7c05fca8-a31b-4bcb-96bb-bd4d8bf41aa1
# f(x) = x^2 * exp(x)
# f(x) = 1 / x
# ϵ = 1e-6
# xlim = 1, 3

# romberg(f, xlim, 10, ϵ)

# ╔═╡ 0da756fa-1da5-4cdd-989e-ed517cf9306c
md"""
### 实验题目
"""

# ╔═╡ cc4e4ed4-5acc-4137-8831-a7f05003414b
md"""
#### 问题 1
"""

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

# ╔═╡ Cell order:
# ╟─1198114b-4ee5-479d-a539-4770886bbaa4
# ╟─f3390402-8050-4b60-b93e-31f42897a94d
# ╟─423c9cdd-ea9f-49ce-9f84-252c63dfa582
# ╟─654ffab7-7e9d-4066-97f5-edc69b1b619f
# ╠═83f83660-08c9-4e28-aedb-81c283097079
# ╟─506efb71-374f-42f8-87e4-2b1f5ef50c06
# ╠═7c05fca8-a31b-4bcb-96bb-bd4d8bf41aa1
# ╟─0da756fa-1da5-4cdd-989e-ed517cf9306c
# ╟─cc4e4ed4-5acc-4137-8831-a7f05003414b
# ╠═5090d873-7ff5-4b0e-8476-b6d788c62462
# ╟─d187517e-5f21-4b4c-ae81-d519777b145d
# ╟─18463200-a64a-400a-a046-b215fd059e7a
# ╟─40e75753-9cdc-44cc-a8e3-b8953d56fb37
