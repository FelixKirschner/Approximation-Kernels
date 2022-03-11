#=
#########  ##   ##  ######  ######  ######
#########  #### ##  ##  ##    ##    ##
#########  ## ####  ##  ##    ##    ####
#########  ##   ##  ######    ##    ######


=#
cd(@__DIR__)
include("libs.jl")
include("SupportFunctions.jl")
include("SymmetryAdaptedBasis.jl")


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

bar = true
for r = 11:12
    list = sigmaSymRedEff(4, r, bar, false, "CSDP");
    writeKernel(list, 4, r, bar)
end

initializeBasis(2,2)