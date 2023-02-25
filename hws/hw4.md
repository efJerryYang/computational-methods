9. 
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

```julia
f(x) = 2 / sqrt(pi) * exp(-x)
ϵ = 1e-6
xlim = 0, 1
romberg(f, xlim, 20, ϵ)
f(x) = exp(-x^2)
ϵ = 1e-6
xlim = 0, 0.8
romberg(f, xlim, 20, ϵ)
```

```
Problem 9.1
 0.771743332	
 0.728069946	 0.713512151	
 0.716982762	 0.713287034	 0.713272026	
Accuracy requirement satisfied.
Problem 9.2
 0.610916970	
 0.646316000	 0.658115677	
 0.654851153	 0.657696204	 0.657668239	
 0.656966396	 0.657671477	 0.657669829	 0.657669854	
Accuracy requirement satisfied.
```

- 9.1，取 $I=0.713272026$

- 9.2，取$I=0.657669854$
