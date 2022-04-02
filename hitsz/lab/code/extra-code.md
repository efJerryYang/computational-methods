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
