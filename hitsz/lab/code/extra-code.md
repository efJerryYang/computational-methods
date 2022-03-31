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