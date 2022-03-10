
function init(n, k)
    β = []
    μ = []

    β = fillmonomials(n, k)

    for el in β
        sort!(el)
    end
    unique!(β)

    for b in β
        btmp = unique(b)
        mu = []
        for el in btmp
            cnt = count(x -> x == el, b)
            push!(mu, cnt)
        end
        sort!(mu, rev = true)
        push!(μ, AbstractAlgebra.Partition([mu[i] for i = 1:length(mu)])) #
    end

    return β, μ
end


# Split permutation into transpositions of form (j,j+1)
function permDecomp(p::Perm)
    res = []
    arr = copy(p.d)
    n = length(arr)
    for i = 1:n
        for j = n:-1:(i+1)
            if arr[j-1] > arr[j]
                push!(res, perm(union(1:(j-2), [j, j - 1], (j+1):n)))
                tmp = arr[j-1]
                arr[j-1] = arr[j]
                arr[j] = tmp
            end
        end
    end
    return res
end

function remove!(a, item)
    deleteat!(a, findall(x->x==item, a))
end


#set entry p in tableau Y to value c
function setTabEntry(Y, p, c)
    linPos = 0
    for i = 1:(p[1]-1)
        linPos += Y.part[i]
    end
    linPos += p[2]
    Y.fill[linPos] = c
end

#fill tableau with entries p
function permTab(p, Y)
    return YoungTableau(Y.part, [p[i] for i in Y.fill])
end

#returns true if tableau Y is semi standard
function isSemiStandard(Y)
    s = size(Y)
    #check rows
    for row = 1:s[1]
        last = Y[row, 1]
        for col = 2:(Generic.rowlength(Y, row, 1)+1)
            if last > Y[row, col]
                return false
            end
            last = Y[row, col]
        end
    end
    # Check columns
    for col = 1:s[2]
        last = Y[1, col]
        for row = 2:(Generic.collength(Y, 1, col)+1)
            if last >= Y[row, col]
                return false
            end
            last = Y[row, col]
        end
    end
    return true
end

#generates all semi standard tableaux of shape lam and content mu
function generateGenSemStd(lam, mu)
    if lam.n != mu.n
        return nothing
    end

    Y = YoungTableau(lam.part)
    mu = mu.part
    content = [i for i = 1:length(mu) for j = 1:mu[i]]
    posComb = unique(permutations(content) |> collect)
    cnt = 0
    semiStdTab = []

    for el in posComb
        tmpY = permTab(el, Y)
        if isSemiStandard(tmpY)
            cnt += 1
            push!(semiStdTab, tmpY)
        end
    end

    return [cnt, semiStdTab]
end

#get b, important in order to construct the specht polynomials
function getb(beta, mu , r)
    if length(beta) != sum(mu)
        return nothing
    end

    b = [0 for it = 1:length(mu)]
    vec = [0 for it = 1:r+1]
    for j = 1:r+1
        vec[j] = count(x -> x == j - 1, beta)
    end
    for i = 1:length(mu)
        indx = findfirst(x -> x == mu[i], vec)
        b[i] = indx - 1
        vec[indx] = 0
    end
    return b
end

#get the Column Stabilizer Set CStabt
function getCStabt(t)
    s = size(t)
    colStabt = []
    colStabtfin = []
    for i = 1:s[2]
        col = t[:, i]
        z = findall(x -> x == 0, col)
        for k = length(z):-1:1
            deleteat!(col, z[k])
        end
        colPer = unique(collect(permutations(col)))
        if i == 1 && s[2] != 1
            colStabt = colPer
            continue
        elseif i == 1 && s[2] == 1
            colStabt = colPer
            return colStabt
        else
            for el in colStabt
                for el2 in colPer
                    push!(colStabtfin, [el; el2])
                end
            end
        end
        if i != s[2]
            colStabt = colStabtfin
            colStabtfin = []
        elseif i == s[2]
            return colStabtfin
        end
    end
end

#returns the sign of a permutation where the entries are [1:n]
function getsignDan(p)
    res = []
    arr = copy(p)
    n = length(arr)
    for i = 1:n
        for j = n:-1:(i+1)
            if arr[j-1] > arr[j]
                push!(res, perm(union(1:(j-2), [j, j - 1], (j+1):n)))
                tmp = arr[j-1]
                arr[j-1] = arr[j]
                arr[j] = tmp
            end
        end
    end
    return iseven(length(res)) ? 1 : -1
end

#returns the sign of an element of CStabt
function getsign(p, tab)
    s = size(tab)
    cols = [0 for i = 1:s[2]]
    cols2 = []
    for i = 1:s[2]
        colcont = []
        cnt = 1
        for j = 2:s[1]
            if tab[j,i] != 0 && j != s[1]
                cnt += 1
                push!(colcont, tab[j,i])
            elseif tab[j,i] != 0 && j == s[1]
                cols[i] = cnt + 1
                push!(colcont, tab[j,i])
            elseif tab[j,i] == 0
                cols[i] = cnt
                push!(cols2, colcont)
                break
            end
        end
    end
    sign = 1
    ColStab = []
    for i = 1:s[2]
        if i == 1
            ColStab = p[1:cols[1]]
        else
            summe = sum(cols[k] for k=1:i-1)
            ColStab = p[summe+1:summe + cols[i]]
        end
        vect = [0 for k = 1:length(ColStab)]
        for j = 1:length(ColStab)
            vect[j] = findfirst( x -> x == tab[j,i], ColStab)
        end
        sign = sign*getsignDan(vect)
    end
    return sign
end

#returns {T} i.e. the row equivalence class of T
function getEquClaT(Tab)
    equClass = []
    list = unique(getCStabt(conj(Tab)))
    for el in list
        tmpTab = YoungTableau(Tab.part, el)
        push!(equClass, tmpTab)
    end
    return equClass
end

#returns sign(sigma) * sigma (X^(t,T)) i think
function getsubBaseEl(tab, STab, perm, bExp)
    n = length(perm)
    tmp = YoungTableau(conj(tab).part, perm)
    tmp = conj(tmp)
    s = size(tmp)
    x = [0 for i = 1:n+1]
    x[n+1] = getsign(perm, tab)
    for i = 1:s[1]
        for j = 1:s[2]
            if tab[i,j] != 0
                ind = tmp[i, j]
                expo = bExp[STab[i, j]]
                x[ind] = expo
            end
        end
    end
    return x
end

#for given lambda, mu, beta this return the subbasis belonging to those
function getPartOfBasis(lambda, mu, beta,r)
    t = YoungTableau(lambda)
    retBase = []
    for i = 1:length(mu)
        if generateGenSemStd(lambda, mu[i])[1] == 0
            continue
        end
        b = getb(beta[i], mu[i].part,r)
        listssTabs = generateGenSemStd(lambda, mu[i])[2]
        for T in listssTabs
            Slist = getEquClaT(T)
            sigmaCSt = getCStabt(t)
            bas = []
            for sigma in sigmaCSt
                for S in Slist
                    x = getsubBaseEl(t, S, sigma, getb(beta[i], mu[i].part,r))
                    push!(bas, x)
                end
            end
            push!(retBase, bas)
        end
    end
    return retBase
end