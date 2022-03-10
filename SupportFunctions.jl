using Combinatorics

function fillmonomials(n, r)
    monlist = []
  for p in Combinatorics.combinations(1:(n+r), r)
      sort!(p)
      c = zeros(Int64, 1, n + 1)
      pos = 1
      lastPos = 0
      for i in p
          pos = pos + (i - lastPos - 1)
          c[pos] += 1
          lastPos = i
      end
      push!(monlist, c[1:n])
  end
  return monlist
end


function flips(gamma)
    n = length(gamma)
    binarr = [[in(a, c) ? 1 : -1 for a = 1:n] for c in combinations(1:n)]
    push!(binarr, [-1 for i = 1:n])
    binarr2 = []
    for el in binarr
        if !in((-1) * el, binarr2)
            push!(binarr2, el)
        end
    end
    
    list = []
    for el in binarr2
        push!(list, [gamma[i] * el[i] for i = 1:n])
    end
    return unique(list)
end

function symmetrizeBaseProduct(baseElement1::Vector{Any}, baseElement2::Vector{Any})
    
    n = length(baseElement1[1]) - 1

    perms = collect(permutations([i for i = 1:n]))

    facn = factorial(n)

    pairs = Vector{Int64}[]

    coeffs = []

    for el1 in baseElement1
        for el2 in baseElement2
            tmp = el1[1:n] - el2[1:n]
            tmp2 = el1[n+1] * el2[n+1]

            for p in perms
                tmp3 = [tmp[p[i]] for i = 1:n]
                push!(coeffs, tmp2)
                push!(pairs, tmp3)
            end
        end
    end

    auxlistPairs = unique(pairs)
    auxlistCoeffs = []
    for el in auxlistPairs
        push!(auxlistCoeffs, sum(coeffs[i] for i in findall(x -> x == el, pairs)))
    end
    retListMon = Vector{Int64}[]
    retListCoeff = []
    for i = 1:length(auxlistCoeffs)
        if auxlistCoeffs[i] != 0.0
            push!(retListMon, auxlistPairs[i])
            push!(retListCoeff, auxlistCoeffs[i] / facn)
        end
    end
    return [retListMon, retListCoeff]
end


function checkOrth(base)
    #test function, are the elements of the base pairwise orthogonal ? true : false 
    orth = true
    for i = 1:length(base)
        for j = i+1:length(base)
            for l = 1:length(base[i])
                for m = 1:length(base[j])
                    list = symmetrizeBaseProduct(base[i][l], base[j][m])
                    if !isempty(list[1])
                        orth = false
                        display((i, j, k, l))
                    end
                end
            end
        end
    end
    return orth
end

function initializeBasis(n, r)
    λ = [AbstractAlgebra.Partition(el) for el in collect(partitions(n))]
    
    βf, μf = init(n, r)
    
    fbas = [getPartOfBasis(λ[i], μf, βf, r) for i = 1:length(λ)]

    for base in fbas
        unique!(base)
    end
    retbas = []
    for base in fbas
        if length(base) != 0
            push!(retbas, base)
        end
    end
    return retbas
end


function orbitSize(mon)
    perms = collect(permutations([i for i = 1:length(mon)]))
    list = []
    for p in perms
        push!(list, [mon[p[i]] for i = 1:length(mon)])
    end
    return length(unique!(list))
end

function writeKernel(list, n::Integer, r::Integer, bar::Bool)
    auxMon = fillmonomials(n, r)
    for el in auxMon
        sort!(el)
    end
    
    unique!(auxMon)
    
    gList = []
    
    for i = 1:length(auxMon)-1
        push!(gList, 0.5 * sum(dot(vec(list[2][j]), vec(list[3][i][j])) for j = 1:length(list[2])))
    end
    
    push!(gList, 1.0)
    
    dict = Dict([auxMon[i] => gList[i] for i = 1:length(auxMon)])
    
    newMon = fillmonomials(n, r)
    
    gListFin = []
    for el in newMon
        push!(gListFin, dict[sort!(el)])
    end
    reverse!(gListFin)
    cd(@__DIR__)
    if bar 
        str1 = "sigmaBarkernel"
    else
        str1 = "sigmakernel"
    end
    str1 *= string(r)
    io = open(str1, "w")
    for i = 1:length(gListFin)
        str = string(gListFin[i])
        write(io, str * "\n")
    end
    close(io)
end

function setUpModel(n::Integer, r::Integer, bar::Bool, silent::Bool, solver::String, dictBaseproducts, baseProducts, fbas, rMons, dictMons)

    vrbl = Array{VariableRef}[]

    geOne = Vector{SparseMatrixCSC{Float64,Int64}}[]

    constrMat = Vector{SparseMatrixCSC{Float64,Int64}}[]
    constrMatTwo = Vector{SparseMatrixCSC{Float64,Int64}}[]
    gcoeffMat = Vector{SparseMatrixCSC{Float64,Int64}}[]



    if solver == "Mosek"
        m = Model(Mosek.Optimizer)
    elseif solver == "CSDP"
        m = Model(CSDP.Optimizer)
    elseif solver == "Hypatia"
        m = Model(Hypatia.Optimizer)
    end

    if silent
        MOI.set(m, MOI.Silent(), true)
    end
    




    @variable(m, t >= 0)
    eOne = [0 for i = 1:n]
    eOne[n] = 1

    if bar
        geTwo = []
        canUnit2 = []
        push!(canUnit2, 2 * eOne)
    end

    @info("Initializing variables...")

    for i = 1:length(fbas)
        if length(fbas[i]) == 0
            continue
        elseif length(fbas[i]) == 1
            NN = @variable(m, [1:1, 1:1])
            push!(vrbl, NN)
            @constraint(m, vrbl[i][1, 1] .>= 0)
        else
            NN = @variable(m, [1:length(fbas[i]), 1:length(fbas[i])], base_name = "y$i", PSD)
            push!(vrbl, NN)
        end
    end

    auxMon = fillmonomials(n, r)

    for el in auxMon
        sort!(el)
    end

    unique!(auxMon)


    @info("Creating matrices for γ in IN_r^n...")

    for gamma in auxMon
        if gamma == [0 for i = 1:n]
            tmplist = SparseMatrixCSC{Float64,Int64}[]
            for i1 = 1:length(fbas)
                push!(
                    tmplist,
                    spzeros(Float64, length(fbas[i1]), length(fbas[i1])),
                )
            end

            for el in dictMons[gamma]
                tmp = dictBaseproducts[(el[1], el[2], el[3])]
                loc = findall(x -> x == [0 for i = 1:n], tmp[1])
                tmplist[el[1]][el[2], el[3]] = sum(tmp[2][loci] for loci in loc)
                tmplist[el[1]][el[3], el[2]] = sum(tmp[2][loci] for loci in loc)
            end


            @constraint(
                m,
                sum(
                    dot(((vrbl[i])), (tmplist[i])) for i = 1:length(vrbl)
                ) == 1
            )
        else

            list = flips(gamma)

            sublist = SparseMatrixCSC{Float64,Int64}[]

            for i1 = 1:length(fbas)
                push!(
                    sublist,
                    spzeros(Float64, length(fbas[i1]), length(fbas[i1])),
                )
            end

            for el in dictMons[gamma] #indicesLookUp[locGamma[1]]
                tmp = dictBaseproducts[(el[1], el[2], el[3])]
                loc = findall(x -> x == gamma, tmp[1])
                sublist[el[1]][el[2], el[3]] = sum(tmp[2][loci] for loci in loc)
                sublist[el[1]][el[3], el[2]] = sum(tmp[2][loci] for loci in loc)
            end

            if gamma == eOne
                push!(geOne, sublist)
                #display(geOne)
            end

            if bar && in(gamma, canUnit2)
                push!(geTwo, sublist)
            end
            push!(gcoeffMat, sublist)

            #now we look at all the flips and make sure we can factor out the g_alpha
            for i = 1:length(list)
                if list[i] == gamma
                    continue
                else
                    addSubList = SparseMatrixCSC{Float64,Int64}[]
                    for cnt = 1:length(fbas)
                        push!(
                            addSubList,
                            spzeros(
                                Float64,
                                length(fbas[cnt]),
                                length(fbas[cnt]),
                            ),
                        )
                    end


                    for el in dictMons[list[i]]
                        symmProd = dictBaseproducts[(el[1], el[2], el[3])]
                        location = findall(
                            x ->
                                x == list[i] ||
                                    x == (-1) * list[i],
                            symmProd[1],
                        )
                        addSubList[el[1]][el[2], el[3]] = sum(
                            (symmProd[2][location[i]]) for
                            i = 1:length(location)
                        )
                        addSubList[el[1]][el[3], el[2]] = sum(
                            (symmProd[2][location[i]]) for
                            i = 1:length(location)
                        )
                    end

                    push!(
                        constrMat,
                        [sublist[i5] - addSubList[i5] for i5 = 1:length(fbas)],
                    )


                end

            end
        end
    end


    @info("Creating matrices for γ in Γ(n,r)...")


    for gamma in rMons

        sublist = SparseMatrixCSC[]
        for i = 1:length(fbas)
            push!(sublist, spzeros(Float64, length(fbas[i]), length(fbas[i])))
        end


        for el in dictMons[gamma]
            tmp = dictBaseproducts[(el[1], el[2], el[3])]
            loc = findall(x -> x == gamma || x == (-1)*gamma, tmp[1])
            sublist[el[1]][el[2], el[3]] = sum(tmp[2][loci] for loci in loc)
            sublist[el[1]][el[3], el[2]] = sum(tmp[2][loci] for loci in loc)
        end

        push!(constrMatTwo, sublist)

    end

    unique!(constrMat)
    unique!(constrMatTwo)

    @info("Adding constraints...")

    for el in constrMat
        boo = true
        for mat in el
            if mat != zeros(Float64, size(mat, 1), size(mat, 2))
                boo = false
            end
        end
        if !boo
            expr = AffExpr()
            for i = 1:length(vrbl)
                @views @inbounds add_to_expression!(expr, vec(vrbl[i])' * vec(el[i]))
            end
            @constraint(
                m,
                expr == 0
            )
        end
    end
    for el in constrMatTwo

        expr = AffExpr()
        for i = 1:length(vrbl)
            @views @inbounds add_to_expression!(expr, vec(vrbl[i])' * vec(el[i]))
        end

        @constraint(
            m,
            expr == 0
        )
    end
    if bar
        if r >= 2
            @constraint(
                m,
                t >=
                1 -
                sum(dot((vrbl[i]), (geOne[1][i])) for i = 1:length(fbas)) +
                1 / 2 *
                sum(dot((vrbl[i]), (geTwo[1][i])) for i = 1:length(fbas))
            )
            @constraint(
                m,
                t >=
                -(
                    1 -
                    sum(
                        dot((vrbl[i]), (geOne[1][i])) for i = 1:length(fbas)
                    ) +
                    1 / 2 *
                    sum(dot((vrbl[i]), (geTwo[1][i])) for i = 1:length(fbas))
                )
            )
        else
            @constraint(
                m,
                t >=
                1 -
                sum(dot(((vrbl[i])), ((geOne[1][i]))) for i = 1:length(fbas))
            )
            @constraint(
                m,
                t >=
                -(
                    1 -
                    sum(
                        dot(((vrbl[i])), ((geOne[1][i]))) for i = 1:length(fbas)
                    )
                )
            )
        end
        obj =
            n * (
                1 / 2 * t + 1 -
                1 / 2 *
                sum(dot(((vrbl[i])), (geOne[1][i])) for i = 1:length(fbas))
            )
        @objective(m, Min, obj)
    else
        obj = sum(dot(((vrbl[i])), ((geOne[1][i]))) for i = 1:length(fbas))
        @objective(m, Max, obj)
    end

    numconstr = (num_constraints(m, AffExpr, MOI.EqualTo{Float64}))
    numconstr += (num_constraints(m, AffExpr, MOI.GreaterThan{Float64}))
    numconstr += (num_constraints(m, VariableRef, MOI.GreaterThan{Float64}))
    blcksizes = Int64[]
    for i = 1:length(vrbl)
        push!(blcksizes, size(vrbl[i], 1))
    end
    @info("Blocksizes: $blcksizes")
    @info("Number of constraints: $numconstr")
    return m, vrbl, obj, gcoeffMat

end


function returnProductsDict(n::Integer, r::Integer)

    dictMons = Dict()
    rMons = []
    baseProducts = []
    dictBaseproducts = Dict()

    @info("Generating symmetry adapted basis")
    fbas = initializeBasis(n, r)

    @info("Generating symmetrized dictBaseproducts")

    for i1 = 1:length(fbas)
        for i2 = 1:length(fbas[i1])
            for i3 = 1:length(fbas[i1])
                symmProd = symmetrizeBaseProduct(fbas[i1][i2], fbas[i1][i3])
                push!(baseProducts, symmProd)
                dictBaseproducts[(i1, i2, i3)] = symmProd
                for el in symmProd[1]
                    tmpSum = sum(abs(el[i]) for i = 1:n)
                    if !haskey(dictMons, el)
                        dictMons[el] = [[i1, i2, i3]]
                    else
                        tmp = dictMons[el]
                        dictMons[el] = push!(tmp, [i1, i2, i3])
                    end
                    if tmpSum > r
                        push!(rMons, el)
                    end
                end
            end
        end
    end
    return baseProducts, dictBaseproducts, fbas, unique(rMons), dictMons
end

