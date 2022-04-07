https://stackoverflow.com/questions/53230969/how-to-scale-a-plot-in-julia-using-plots-jl
```
(xlim=(xMin,xMax), ylim=(yMin, yMax), yflip = false)
```
```
for j in filter(x -> x != i, 1:num)  # much slower than continue
```
作图添加函数表达式的位置
```
xmin, xmax = xlims(p)
x = xmax + (xmax - xmin) * 0.04
y = mean(ylims(p))
ymax = ylims(p)[2]
annotate!(x, y, text, :black)
```


DataFrames使用：https://dataframes.juliadata.org/stable/man/getting_started/

Dataframe展示，代码写法
https://discourse.julialang.org/t/display-more-decimals-in-dataframe/22545
```
julia> Base.show(io::IO, t::Float64) = @printf io  "%1.9f" t

julia> x = DataFrame(a=[1.142300004,1.142300051])
2×1 DataFrame
│ Row │ a           │
│     │ Float64     │
├─────┼─────────────┤
│ 1   │ 1.142300004 │
│ 2   │ 1.142300051 │
```

!不需要换元，也不需要改xspan，当x最左边的值确定的时候，就是初值的位置

```julia
# from: https://stackoverflow.com/questions/58667332/is-there-a-way-to-swap-columns-in-o1-in-julia
function swapcols!(X::AbstractMatrix, i::Integer, j::Integer)
    @inbounds for k = 1:size(X, 1)
        X[k, i], X[k, j] = X[k, j], X[k, i]
    end
end
# from: https://discourse.julialang.org/t/swap-cols-rows-of-a-matrix/47904/9
function _swapcol!(x, i, j)
    for k in axes(x, 1)  # <- give dimension as input to axes function
        x[k, i], x[k, j] = x[k, j], x[k, i]
    end
end

```

```julia
# https://stackoverflow.com/questions/45396685/what-does-an-exclamation-mark-mean-after-the-name-of-a-function
# https://people.richland.edu/james/lecture/m116/matrices/pivot.html
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
```

```julia
# https://stackoverflow.com/questions/62142717/julia-quick-way-to-initialise-an-empty-array-thats-the-same-size-as-another
    x = similar(b, Float64)
    
```

```julia
# x = [@elapsed(rand(i,i)\rand(i)) for i = 10:15]
# for i = 10:15
#     A, b = rand(i,i), rand(i)
#     @time A\b
# end
# # compute moving average
# # https://stackoverflow.com/questions/28820904/how-to-efficiently-compute-average-on-the-fly-moving-average
# # n=1;
# # curAvg = 0;
# # loop{
# #   curAvg = curAvg + (newNum - curAvg)/n;
# #   n++;
# # }
# n = 2000
# curAvg = zeros(n)
# # curAvg[i] = 0
# for i = 2:n
#     curAvg[i] = curAvg[i-1] + (@elapsed(rand(i,i)\rand(i)) - curAvg[i-1])/(i-1)
# end
# plot(curAvg)

# function ma(n)
#     curAvg = zeros(n)
#     # curAvg[i] = 0
#     for i = 2:n
#         curAvg[i] = curAvg[i-1] + (@elapsed(rand(i,i)\rand(i)) - curAvg[i-1])/(i-1)
#     end
#     plot(curAvg)
# end

```

```julia
function ema!(x)
    curAvg = x[1]
    n = size(x,1)
    for i=2:n
        x[i] = curAvg + (x[i] - curAvg)/(i-1)
        curAvg = x[i]
    end
    x
end

function sma!(x)
    k = 20
    n = size(x,1)
    for i=1:k-1
        x[i] = mean(x[1:i])
    end
    for i=k:n
        x[i] = mean(x[i-k+1:i])
    end
    x
end

```
