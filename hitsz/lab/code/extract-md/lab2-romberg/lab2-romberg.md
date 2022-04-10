## 实验题目2 龙贝格(Romberg)积分法

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



### 数学原理

教材中给出的计算公式如下
$$
\begin{cases}
	T_{0,0}&=\frac{b-a}{2}\left[ f\left( a \right) +f\left( b \right) \right] ,\\
	T_{0,i}&=\frac{1}{2}T_{0,i-1}+\frac{1}{2}\frac{b-a}{2^{i-1}}\sum_{j=1}^{2^{i-1}}{f\left[ a+\left( j-\frac{1}{2} \right) \cdot \frac{b-a}{2^{i-1}} \right]}, i=1,2,3,...,\\
	T_{m,k}&=\frac{4^mT_{m-1,k+1}-T_{m-1,k}}{4^m-1},m=1,2,...,i; k=i-m.\\
\end{cases}
$$
因`Julia`语言数组类下标的起点为1，同时实验指导书所给T数表为下三角形，故将原公式改写如下
$$
\begin{aligned}
\begin{cases}
	T_{1,1}&=\frac{b-a}{2}\left[ f\left( a \right) +f\left( b \right) \right] ,\\
	T_{i+1,1}&=\frac{1}{2}T_{i,1}+\frac{1}{2}\frac{b-a}{2^{i-1}}\sum_{j=1}^{2^{i-1}}{f\left[ a+\left( j-\frac{1}{2} \right) \cdot \frac{b-a}{2^{i-1}} \right]}, i=1,2,3,...,n\\
	T_{i+1,m+1}&=\frac{4^mT_{i+1,m}-T_{i,m}}{4^m-1},m=1,2,...,i.\\
\end{cases}
\end{aligned}
$$
随后可对照公式完成代码的编写



### 代码实现

使用`Julia`编程语言，根据上述数学原理，编写`romberg`积分法实验代码。

以下部分为`romberg()`函数定义：


```julia
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
```



### 测试代码

本部分使用教材上的例题用于对程序结果进行初步的检验，计算结果和教材给出数表类似，可以认为测试通过。


```julia
f(x) = x^2 * exp(x)
f(x) = 1 / x
ϵ = 1e-6
xlim = 1, 3

romberg(f, xlim, 10, ϵ)
```

     1.333333333	
     1.166666667	 1.111111111	
     1.116666667	 1.100000000	 1.099259259	
     1.103210678	 1.098725349	 1.098640372	 1.098630548	
     1.099767702	 1.098620043	 1.098613022	 1.098612588	 1.098612518	
    Accuracy requirement satisfied.

### 实验题目

#### 问题 1


```julia
iter_num = 30

f(x) = x^2 * exp(x)
ϵ = 1e-6
xlim = 0, 1
println("f(x) = x^2 * exp(x)")
romberg(f, xlim, iter_num, ϵ)

f(x) = exp(x)sin(x)
ϵ = 1e-6
xlim = 1, 3
println("f(x) = exp(x)sin(x)")
romberg(f, xlim, iter_num, ϵ)


f(x) = 4 / (1 + x^2)
ϵ = 1e-6
xlim = 0, 1
println("f(x) = 4 / (1 + x^2)")
romberg(f, xlim, iter_num, ϵ)

f(x) = 1 / (x + 1)
ϵ = 1e-6
xlim = 0, 1
println("f(x) = 1 / (x + 1)")
romberg(f, xlim, iter_num, ϵ)
```

    f(x) = x^2 * exp(x)
     1.359140914	
     0.885660616	 0.727833850	
     0.760596332	 0.718908238	 0.718313197	
     0.728890177	 0.718321459	 0.718282340	 0.718281850	
    Accuracy requirement satisfied.
    
    f(x) = exp(x)sin(x)
     5.121826420	
     9.279762907	10.665741736	
    10.520554284	10.934151409	10.952045388	
    10.842043468	10.949206529	10.950210203	10.950181074	
    10.923093890	10.950110697	10.950170975	10.950170352	10.950170310	
    Accuracy requirement satisfied.
    
    f(x) = 4 / (1 + x^2)
     3.000000000	
     3.100000000	 3.133333333	
     3.131176471	 3.141568627	 3.142117647	
     3.138988494	 3.141592502	 3.141594094	 3.141585784	
     3.140941612	 3.141592651	 3.141592661	 3.141592638	 3.141592665	
    Accuracy requirement satisfied.
    
    f(x) = 1 / (x + 1)
     0.750000000	
     0.708333333	 0.694444444	
     0.697023810	 0.693253968	 0.693174603	
     0.694121850	 0.693154531	 0.693147901	 0.693147478	
    Accuracy requirement satisfied.


​    

### 思考题

1. 略

2. 在实验 1 中二分次数和精度的关系如何？

   二分次数越多所求的精度越高，通常预设较大的二分次数来确保计算结果有足够的精度，同时也设定早停需要满足的精度要求，避免达到所需精度之后继续计算导致增加的运算量

3. 略

4. 略

### 参考资料

1. julia 数值积分 https://blog.csdn.net/m0_37816922/article/details/103475445
2. Romberg Integration-Numerical Analysis http://homepages.math.uic.edu/~jan/mcs471/romberg.pdf
3. 《数值分析》吴勃英 196-199

