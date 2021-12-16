function make_g(D)
    u = ones(D)
    u[1] = -1.0
    g = Diagonal(u) + zeros(D,D)
end

#g = make_g(4)
#latexify(g)
#typeof(g)

function vector2symmat!(v,sm)
    M = length(sm[1,:]) # N = M*(M+1)/2 
    N = length(v)
    if M*(M+1)÷2 != N
        error("Dimension missmatch, M = $M, N= $N")
    end
    L = 1
    for i in 1:M
        for j in i:M
            sm[i,j] = v[L]
            sm[j,i] = sm[i,j]
            L=L+1
        end
    end
end




function symmat2vector(sm)
    M = length(sm[1,:])
    N = M*(M+1)÷2
    v = zeros(N)
    l = 1
    for i in 1:M
        for j in i:M
            v[l] = sm[i,j]
            l = l+1
        end
    end
    return v
end


function vector2symmat(v)
    N = length(v)
    MM = (-1 + sqrt(1+8N))/2 
    if MM % 1 ≈ 0
        M = Int64(MM)
    else
        error("the vector dimension does not corresponds to the number of free entries of a symmetric matrix")
    end
    #println("M = $M")
    #sm = Array{Real,2}(M,M)
    sm = zeros(M,M)
    
    L = 1
    for i in 1:M
        for j in i:M
            sm[i,j] = v[L]
            sm[j,i] = sm[i,j]
            L=L+1
        end
    end
    return sm
end


function vector_unpack(ζ)
    N = length(Symbolics.scalarize(ζ)) # N = D + D*(D+1)/2
    DD = (-3 + sqrt(9+8N))/2 
    if DD % 1 ≈ 0
        D = Int64(DD)
    else
        error("In function vector_unpack: the vector dimension does not corresponds to the number of free entries of a symmetric matrix N = $N, D = $DD")
    end
    return Symbolics.scalarize(ζ)[1:D], vector2symmat(Symbolics.scalarize(ζ)[D+1:end])[:,:]
end

#vector_unpack(rand(14))
#vector_unpack(rand(5))

function upper2symm!(s)
        M = length(s[:,1])
        for i ∈ 1:M
                for j ∈ i:M
                        s[j,i] = s[i,j]
                end
        end
end

function upper2symm(s)
        ss = copy(s)
        M = length(s[:,1])
        for i ∈ 1:M
                for j ∈ i:M
                        ss[j,i] = s[i,j]
                end
        end
        return ss[:,:]
end

function symm_A(AA)
    AS = zeros(D,D,D)
    for k ∈ 1:D
        for j ∈ 1:D
            for i ∈ 1:D
                AS[i,j,k] = (AA[i,j,k] + AA[j,i,k])/2
            end
        end
    end
    return AS[:,:,:]
end

function tr_A(A,g)
    D = length(A[1,1,:])
    TrAS = zeros(D)
    for k ∈ 1:D
        TrAS[k] = tr(g*A[:,:,k])
    end
    return TrAS[:]
end

function STF_A(A,g)
    D = length(A[1,1,:])
    A_STF = zeros(D,D,D)
    for k ∈ 1:D
        T = tr_A(A,g)[k]
        A_STF[:,:,k] = A[:,:,k] - T*g/D 
    end
    return A_STF[:,:,:]
end

function full_symmetry_check(A)
    D = length(A[1,1,:])
    E = 0
    for i ∈ 1:D
        for j ∈ 1:D
            for k ∈ 1:D
                E = E = abs(A[i,j,k] - A[j,i,k]) + abs(A[i,j,k] - A[k,j,i])
            end
        end
    end
    return E
end


# Functions

function Φ_vt(ζ_v, ζ_s, p)
    D = length(ζ_v)
    g = make_g(D)
    #println(g)
    χ₀ = p[1]
    χ₁ = p[2]
    χ₂ = p[3]
    #ζ_s = my_symmetrization!(Symbolics.scalarize(ζ_t))
    μ = Symbolics.scalarize(ζ_v)'*g*Symbolics.scalarize(ζ_v)
    ν = (g*Symbolics.scalarize(ζ_v))'*Symbolics.scalarize(ζ_s)*Symbolics.scalarize(g*ζ_v)
    l = Symbolics.scalarize(ζ_s)*Symbolics.scalarize(g*ζ_v)
    𝚶 = l'* g * l
    τ₂= tr(Symbolics.scalarize(ζ_s)*g'*Symbolics.scalarize(ζ_s)*g)
    D1 = D//2+1
    D2 = D//2+2
    return χ₀*μ^(1-D//2) + χ₁*ν*μ^(-1-D//2) + χ₂*(τ₂ - 4*D1*𝚶*μ^(-1) + 2*D1*D2*ν^2*μ^(-2))*μ^(-1-D//2) 
    # the function seems to be correct: tr(g*T)=0 and also A is fully symmetric after symmetrization 
    # on first two indices
    #return 𝚶
end

function Φ_v(ζ,p)
    ζ_v, ζ_t = vector_unpack(ζ)
    #ζ_t = vector2symmat(ζ[D+1:end])
    #return Φ_vt(ζ[1:D], ζ_t, p)
    return Φ_vt(ζ_v, ζ_t, p)
end

Φ_vt(ζ_v,ζ_t, [1.,1.,1.])#, ζ_t, [1.;1;1])
Φ_v(ζ, [1.,1.,1.])
#ζ_v

function Φ(ζ,p)
    #χ = zeros(3)
    χs, gs = p
    χ = Symbolics.scalarize(χs)
    g = Symbolics.scalarize(gs)
    #D = length(g[:,1])
    ζ_v, ζ_t = vector_unpack(ζ)
    μ = Symbolics.scalarize(ζ_v' * gs * ζ_v)
    #μ = Symbolics.scalarize(ζ_v)'*Symbolics.scalarize(gs)*Symbolics.scalarize(ζ_v)
    #ν = (g*Symbolics.scalarize(ζ_v))'*Symbolics.scalarize(ζ_s)*Symbolics.scalarize(g*ζ_v)
    #l = Symbolics.scalarize(ζ_s)*Symbolics.scalarize(g*ζ_v)
    #𝚶 = l'* g * l
    #τ₂= tr(Symbolics.scalarize(ζ_s)*g'*Symbolics.scalarize(ζ_s)*g)
    D1 = D//2+1
    D2 = D//2+2
    return χ[1]*μ^(1-D//2)# + χ₁*ν*μ^(-1-D//2) + χ₂*(τ₂ - 4*D1*𝚶*μ^(-1) + 2*D1*D2*ν^2*μ^(-2))*μ^(-1-D//2) 
end
    

function get_dim(N)
    DD = (-3 + sqrt(9+8N))/2 
    if DD % 1 ≈ 0
        return D = Int64(DD)
    else
        error("In function vector_unpack: the vector dimension does not corresponds to the number of free entries of a symmetric matrix N = $N, D = $DD")
    end
end

    
function mu(v)
    D = get_dim(length(v))
    g = make_g(D)
    return v[1:D]' * g * v[1:D]
end

function lambda(v)
    D = get_dim(length(v))
    g = make_g(D)
    return v[1:D]' * g * v[D+1:end]
end

#"""
#trasnforms indices of a MxM symmetric matrix into a vector index
#"""
l_ind(i::Int64,j::Int64,M::Int64)::Int64
    if i <= j 
    return Int64((i-1)*M - (i-1)*i÷2 + j)
    else
    return Int64((j-1)*M - (j-1)*j÷2 + i)
end




