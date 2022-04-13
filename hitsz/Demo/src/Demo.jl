module Demo
using LaTeXStrings
using Markdown
using Printf
using LinearAlgebra
using LaTeXStrings
using PrettyTables
using NLsolve
using Roots
using SymPy
using PyCall
function julia_main()::Cint
    len = length(ARGS)
    if len == 0
        println("An option is needed!")
        println("Type --help or -h to see available options.")
    else
        if ARGS[1] == "--help" || ARGS[1] == "-h"
            print("""Demo.exe [switches]
             -h, --help\t\t\tPrint this help infomation
             -1, --lab1, --lagrange\t\tDisplay lab1-lagrange result
             -2, --lab2, --romberg\t\tDisplay lab2-romberg result
             -3, --lab3, --rungekutta\tDisplay lab3-rungekutta result
             -4, --lab4, --newton\t\tDisplay lab4-newton result
             -5, --lab5, --gauss\t\tDisplay lab5-gauss result
             -a, --all\t\t\tDisplay all the results
            """)
        elseif ARGS[1] == "--lab1" || ARGS[1] == "-1" || ARGS[1] == "--lagrange"
            println()
            include("lab1-lagrange.jl")
        elseif ARGS[1] == "--lab2" || ARGS[1] == "-2" || ARGS[1] == "--romberg"
            println()
            include("lab2-romberg.jl")
        elseif ARGS[1] == "--lab3" || ARGS[1] == "-3" || ARGS[1] == "--runge-kutta" || ARGS[1] == "--rungekutta"
            println()
            include("lab3-rungekutta.jl")
        elseif ARGS[1] == "--lab4" || ARGS[1] == "-4" || ARGS[1] == "--newton"
            println()
            include("lab4-newton.jl")
        elseif ARGS[1] == "--lab5" || ARGS[1] == "-5" || ARGS[1] == "--gauss"
            println()
            include("lab5-gauss.jl")
        elseif ARGS[1] == "--all" || ARGS[1] == "-a"
            println()
            include("lab1-lagrange.jl")
            println()
            include("lab2-romberg.jl")
            println()
            include("lab3-rungekutta.jl")
            println()
            include("lab4-newton.jl")
            println()
            include("lab5-gauss.jl")
            # invalid option
        else
            println("Invalid option!")
        end
    end
    return 0
end
if abspath(PROGRAM_FILE) == @__FILE__
    julia_main()
end
end # module

