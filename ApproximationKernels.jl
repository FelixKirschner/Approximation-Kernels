#=
#########  ##   ##  ######  ######  ###### 
#########  #### ##  ##  ##    ##    ##      
#########  ## ####  ##  ##    ##    ####    
#########  ##   ##  ######    ##    ######  

maybe improve function names
=#

cd(@__DIR__)
include("libs.jl")
include("SupportFunctions.jl")
include("SymmetryAdaptedBasis.jl")

foo(-2)
##
function sigmaSymRedEff(n::Integer, r::Integer, fixr::Integer, silent::Bool, solver = "Mosek", bar = false)

    @time baseProducts, dictBaseproducts, fbas, rMons, dictMons = returnProductsDict(n, r, fixr)

    @info("Setting up model...")

    @time m, vrbl, obj, gcoeffMat = setUpModel(n, r, fixr, bar, silent, solver, dictBaseproducts, baseProducts, fbas, rMons, dictMons)
    
    @info("Optimize called...")
    
    @time optimize!(m)

    display((r, termination_status(m)))

    if bar
        result = value.(obj)
    else
        result = n - (n / 2) * value.(obj)
    end

    listVars = []

    for i = 1:length(vrbl)

        push!(listVars, value.(vrbl[i]))

    end

    @info("Objective value: $result")

    return result, listVars, gcoeffMat

end


#=
Let us have a look at an example.

First fix n and r to the desired values
=#


n = 2
r = 4
fixr = r
#bar = false
silent = false
@time list = sigmaSymRedEff(n, r, r, silent, "Mosek");

#=

=#

writeKernel(list, n, fixr, false)

compList = []

compList2 = []
for fixr = 5:9##
    rList = []
    for r = 19:20
        push!(rList, (r,sigmaSymRedEff(n, r, fixr, true, "CSDP")[1]))
    end
    push!(compList2, rList)
end

show(compList)

function gJ(n, r)
    return ((r + 1 - n + 1) * cos((pi * n) / (r + 2)) + sin((pi * n) / (r + 2)) * cot(pi / (r + 2))) / (r + 2)
end

1-gJ(1,12)

using Plots
for el in compList
    x = []
    y = []
    for i = 1:length(el)
        push!(x, el[i][1])
        push!(y, el[i][2])
    end
    if x[1] >= 3
        if x[1] == 3
            q = x[1]
            plot(x, y, xaxis = ("r", (0,20), 0:1:20), yaxis = "Resolution", label = "r' = $q")
        else
            q=x[1]
            plot!(x, y, xaxis = ("r", (0,20), 0:1:20), yaxis = "Resolution", label = "r' = $q")
        end
    end
end

savefig("test_resolution_comp_4")

show(compList)

for el in compList
    for el2 in el
        print("(",el2[1],",",round(el2[2]; digits = 4),")")
    end
    println("\n \n \n")
end

x = [1,2]
y = [2,3]
plot(x,y)

el = compList[5]

x = []
tmp = []
for i = 1:length(el)
    push!(tmp, el[i][2])
    push!(x, el[i][1]+i-1)
end
plot(x,tmp)

list

writeKernel(list, 2, 10, false)
vals = []
for i = 1:length(compList)
    push!(vals, compList[i][end][2])
end

plot([i for i = 1:10],vals)

B = initializeBasis(2, 50)
B[2]
binomial(3,2)

rap(n,r) = n*(1-cos(n*pi/(n+r)))
rap(6,6)