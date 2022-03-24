using Printf
function newton(f::Function, df::Function, x0, ϵ1, ϵ2, N)
    n = 1
    while n<=N
        F= f(x0)
        DF = df(x0)
        if abs(F) < ϵ1
            @printf("%4d%18.12f\n",n-1,x0)
            return
        end
        if abs(DF) < ϵ2
            @printf("Reach a critical point!\n")
            return
        end
        x1 = x0 - F/DF
        Tol = abs(x1-x0)
        if Tol < ϵ1
            @printf("%4d%18.12f\n",n-1,x1)
            return
        end
        n = n +1
        x0 = x1
    end
    @printf("Fail to converge within %d iterations!\n",N)
end