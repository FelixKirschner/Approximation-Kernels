#=
#########  ##   ##  ######  ######  ###### 
#########  #### ##  ##  ##    ##    ##      
#########  ## ####  ##  ##    ##    ####    
#########  ##   ##  ######    ##    ######  

No notes
=#

cd(@__DIR__)
include("libs.jl")
include("SupportFunctions.jl")
include("SymmetryAdaptedBasis.jl")

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


n = 2
r = 4
fixr = r
#bar = false
silent = false
@time list = sigmaSymRedEff(n, r, r, silent, "Mosek");
