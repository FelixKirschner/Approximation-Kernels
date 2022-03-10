#=
#########  ##   ##  ######  ######  ######
#########  #### ##  ##  ##    ##    ##
#########  ## ####  ##  ##    ##    ####
#########  ##   ##  ######    ##    ######

This code 

=#
cd(@__DIR__)
include("SupportFunctions.jl")
include("SymmetryAdaptedBasis.jl")
using JuMP
using MosekTools
using CSDP
using Hypatia

##
function sigmaSymRedEff(n::Integer, r::Integer, bar::Bool, silent::Bool, solver = "Mosek")


    @time baseProducts, dictBaseproducts, fbas, rMons, dictMons = returnProductsDict(n, r)

    @info("Setting up model...")
    @time m, vrbl, obj, gcoeffMat = setUpModel(n, r, bar, silent, solver, dictBaseproducts, baseProducts, fbas, rMons, dictMons)

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
    return value.(obj), listVars, gcoeffMat


end


n = 2
r = 15
bar = true
@time list = sigmaSymRedEff(2,13, true, false, "CSDP");
writeKernel(list, n, r, bar)
