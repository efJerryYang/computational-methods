using Markdown

# ╔═╡ 1d7eb8a3-b4b0-4ce8-bcac-24d3c261146d
display(md"""
## 实验题目5 高斯(Gauss)列主元消去法
""")

using Printf
using LinearAlgebra
using PrettyTables


function swaprows!(X::AbstractMatrix, i::Integer, j::Integer)
    @inbounds for k = 1:size(X, 2)
        X[i, k], X[j, k] = X[j, k], X[i, k]
    end
end

function pivoting!(A::Matrix{Float64}, k::Integer, n::Integer)
    val, idx = findmax(A[k:n, k])
    idx += k - 1  # index must add previous length that omitted by slice operator
    return val, idx
end
function pivoting!(A::Matrix{Float64}, b::Vector{Float64}, k::Integer, n::Integer, implicit::Bool)
    s = [maximum(A[i, k:n]) for i in k:n]
    if 0 in s
        println("Cannot solve a singular matrix!")
        return
    end
    if implicit
        val, idx = findmax(A[k:n, k] ./ s[1:n-k+1])
    else
        A[k:n, k:n] = A[k:n, k:n] ./ s
        b[k:n] = b[k:n] ./ s
        val, idx = findmax(A[k:n, k])
    end
    idx += k - 1  # index must add previous length that omitted by slice operator
    return val, idx
end

# Gauss列主元消去法
function gauss(n, A::Matrix{Float64}, b::Vector{Float64})
    for k = 1:n-1
        # select pivot in columns
        val, idx = pivoting!(A, k, n)
        if val == 0
            println("Cannot solve a singular matrix!")
            return
        end
        # swap rows
        if idx != k
            swaprows!(A, idx, k)
            b[idx], b[k] = b[k], b[idx]
        end
        # elimination
        for i = k+1:n
            m = A[i, k] / A[k, k]
            A[i, :] -= A[k, :] * m
            b[i] -= b[k] * m
        end
    end
    if A[n, n] == 0
        println("Cannot solve a singular matrix!")
        return
    end
    x = similar(b, Float64)
    x[n] = b[n] / A[n, n]
    for k = n-1:-1:1  # the usage of reverse sequence
        x[k] = (b[k] - dot(A[k, k+1:n], x[k+1:n])) / A[k, k]  # something really annoying 
    end
    x
end


# # ╔═╡ b0fe18e0-a512-4345-8181-d4577c2d0efb
# md"""
# #### explicit/implicit
# """

# # ╔═╡ 1f11b02a-1e31-4369-be69-b97015fcdcaa
# # Gauss列主元消去法
# function gauss(n, A::Matrix{Float64}, b::Vector{Float64}, implicit::Bool)
#     for k = 1:n-1
#         # select pivot in columns
#         val, idx = pivoting!(A, b, k, n, implicit)
#         if val == 0
#             println("Cannot solve a singular matrix!")
#             return
#         end
#         # swap rows
#         if idx != k
#             swaprows!(A, idx, k)
#             b[idx], b[k] = b[k], b[idx]
#         end
#         # elimination
#         for i = k+1:n
#             m = A[i, k] / A[k, k]
#             A[i, :] -= A[k, :] * m
#             b[i] -= b[k] * m
#         end
#     end
#     if A[n, n] == 0
#         println("Cannot solve a singular matrix!")
#         return
#     end
#     x = similar(b, Float64)
#     x[n] = b[n] / A[n, n]
#     for k = n-1:-1:1  # the usage of reverse sequence
#         x[k] = (b[k] - dot(A[k, k+1:n], x[k+1:n])) / A[k, k]  # something really annoying 
#     end
#     x
# end

# ╔═╡ 63034364-b678-4831-8ec1-0b1f8a193d43
md"""
#### 执行代码

本部分是实验代码进行运行时封装的部分，将函数的调用细节隐藏在show_result()函数内部，便于直接从外部使用特定参数对函数进行调用。
"""

# ╔═╡ 258b5589-2aae-485e-b884-0f0562db5627
function show_result(A, b)
    n = size(A, 1)
    A = Float64.(A)
    b = Float64.(b)
    
    println("input matrix: [A | b]")
    data=[A repeat([|], inner=(n, 1)) b]
    pretty_table(
        data;
        header_crayon=crayon"bold",
        tf = tf_matrix,
        noheader=true,
        formatters=ft_printf("%11.8f"))
    
    println("lib solver result:")
    pretty_table(
        @time A \ b;
        header_crayon=crayon"bold",
        tf = tf_matrix,
        noheader=true,
        formatters=ft_printf("%11.8f"))
    
    println("my gauss solver result:")
    pretty_table(
        @time gauss(n, A, b);
        header_crayon=crayon"bold",
        tf = tf_matrix,
        noheader=true,
        formatters=ft_printf("%11.8f"))
end

# ╔═╡ d6b9467c-8b94-4a44-9fe1-0a6483c151c0
md"""
#### 问题 1
"""

# ╔═╡ cd8e7161-b62e-4ee9-9fd5-a52d9192716d
md"""
##### 1.1
"""
println()
print("Problem 1.1\t")
# ╔═╡ 6aef44f6-190c-468a-ba6c-6015bb7e1978
A = [0.4096 0.1234 0.3678 0.2943
     0.2246 0.3872 0.4015 0.1129
     0.3645 0.1920 0.3781 0.0643
     0.1784 0.4002 0.2786 0.3927]
b = [1.1951; 1.1262; 0.9989; 1.2499]
show_result(A, b)

# ╔═╡ a03fb506-1004-458f-bccb-bf856dfd8442
md"""
##### 1.2
"""
println()
print("Problem 1.2\t")
# ╔═╡ 20649c76-4c82-462f-8d1d-ce6e9bdd18f0
A = [136.01  90.860       0       0
     90.860  98.810 -67.590       0
          0 -67.590  132.01  46.260
          0       0  46.260  177.17]
b = [226.87; 122.08; 110.68; 223.43]
show_result(A, b)

# ╔═╡ db024297-2c38-4a53-93a3-2e04e9eebf85
md"""
##### 1.3
"""
println()
print("Problem 1.3\t")
# ╔═╡ 313befa9-50c2-4a40-8722-ff919258293d
A = [  1 1/2 1/3 1/4
     1/2 1/3 1/4 1/5
     1/3 1/4 1/5 1/6
     1/4 1/5 1/6 1/7]
b = [25 / 12; 77 / 60; 57 / 60; 319 / 420]
show_result(A, b)

# ╔═╡ a1771b50-4842-4608-9ae9-17beb86f9664
md"""
##### 1.4
"""
println()
print("Problem 1.4\t")
# ╔═╡ bc47fa10-f5f8-4c11-9b78-fd994bfd549e
A = [10  7  8  7
      7  5  6  5
      8  6 10  9
      7  5  9 10]
b = [32; 23; 33; 31]
show_result(A, b)

# ╔═╡ 6a89bbfa-608e-45d0-8ba6-28fbb7790638
md"""
#### 问题 2

"""

# ╔═╡ 07090ad7-24a6-400a-9a26-5da0ffa35f12
md"""
##### 2.1
"""
println()
print("Problem 2.1\t")
# ╔═╡ fc6d4d4a-05bf-4670-8157-a4398dc7ff4f
A = [ 197   305  -206  -804
     46.8  71.3 -47.4  52.0
     88.6  76.4 -10.8  802
     1.45  5.90  6.13  36.5]
b = [136; 11.7; 25.1;  6.60]
show_result(A, b)

# ╔═╡ 32a45274-3924-465e-8960-c07f51d9f82b
md"""
##### 2.2
"""
println()
print("Problem 2.2\t")
# ╔═╡ 22c1a3fc-4da4-4283-85c9-25d4ea3a06e6
A = [0.5398  0.7161 -0.5554 -0.2982
     0.5257  0.6924  0.3565 -0.6255
     0.6465 -0.8187 -0.1872  0.1291
     0.5814  0.9400 -0.7779 -0.4042]
b = [0.2058; -0.0503; 0.1070; 0.1859]
show_result(A, b)

# ╔═╡ 4672ca8c-44b4-4688-9c29-1b8003ccfacc
md"""
##### 2.3
"""
println()
print("Problem 2.3\t")
# ╔═╡ 3f4f2f5e-9418-4180-b147-01efb9e3d0f5
A = [10  1  2
      1 10  2
      1  1  5]
b = [13; 13; 7]
show_result(A, b)

# ╔═╡ 6e217fa4-6d31-4e07-8fa9-474c3175f5bd
md"""
##### 2.4
"""
println()
print("Problem 2.4\t")
# ╔═╡ e42ed867-374c-4c14-93fd-c2d743d7fdd2
A = [4 -2 -4
    -2 17 10
    -4 10  9]
b = [-2; 25; 15]
show_result(A, b)

# ╔═╡ 9670d76b-abef-41be-991e-771351ce06ce
md"""
### 总结

本实验为高斯列主元消元法的实现，以及运用所写代码完成问题的求解，同时熟悉了`Julia`语言的一些内置函数，以提高代码运行效率的细节用法，例如：
 - 使用`findmax()`同时获取向量最大元素及其下标，
 - 使用`@inbounds`宏减少不必要的边界检查，以节约时间
 - 使用`similar()`返回内容任意形状相同的矩阵或者向量
 - 使用`norm()`和`opnorm()`分别计算向量和矩阵的范数
 - 使用`@elapsed`宏获取对应代码运行的时间数值

本次实验代码随能正确完成任务，且在矩阵阶数较低时有高于库函数的运行效率，但较为零散和频繁的内存分配是的函数在处理高阶矩阵时的时间效率远远低于库函数解法。主要的优化方向为分辨出代码中不必要的内存分配部分，但由于当前`Julia`几乎就是默认使用引用传递，不存在`python`语言中的内存消耗问题，优化可能还需要参考编译出的汇编代码，可在函数前加`@code_native`宏来查看。

这是本实验课的最后一个实验，虽然并非最后完成的实验，但在实验报告的整理过程中，通过之前实验完成过程中学习的`PrettyTables.jl`库重写了`show_result()`函数以获得更好的结果呈现方式。

通过这半学期的练习，对于`Julia`语言各个领域库的使用有了基本的了解，方便了日后的深入使用。
"""

# ╔═╡ 8abcee4c-7840-4ee1-b287-583921330037
display(md"""
### 参考资料

1. julia swapcols fast https://stackoverflow.com/questions/58667332/is-there-a-way-to-swap-columns-in-o1-in-julia
2. julia \_swapcol fast https://discourse.julialang.org/t/swap-cols-rows-of-a-matrix/47904/9
3. julia pivoting https://stackoverflow.com/questions/45396685/what-does-an-exclamation-mark-mean-after-the-name-of-a-function
4. julia pivoting https://people.richland.edu/james/lecture/m116/matrices/pivot.html
5. julia similar https://stackoverflow.com/questions/62142717/julia-quick-way-to-initialise-an-empty-array-thats-the-same-size-as-another
6. moving average pseudocode https://stackoverflow.com/questions/28820904/how-to-efficiently-compute-average-on-the-fly-moving-average
7. julia repeat method https://www.geeksforgeeks.org/creating-array-with-repeated-elements-in-julia-repeat-method/
8. julia repeat usage http://www.jlhub.com/julia/manual/en/function/repeat
9. moving average https://stackoverflow.com/questions/12636613/how-to-calculate-moving-average-without-keeping-the-count-and-data-total
10. 《数值分析原理》吴勃英 46-48
11. 《计算方法实验指导》实验题目 5 高斯(Gauss)列主元消去法
""")
