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
    #return Symbolics.scalarize(ζ)[1:D], vector2symmat(Symbolics.scalarize(ζ)[D+1:end])[:,:]
    return ζ[1:D], vector2symmat(ζ[D+1:end])[:,:]
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

function tr_A_12(A,g)
    D = length(A[1,1,:])
    TrAS = zeros(D)
    for k ∈ 1:D
        TrAS[k] = tr(g*A[:,:,k])
    end
    return TrAS[:]
end

function tr_A_13(A,g)
    D = length(A[1,1,:])
    TrAS = zeros(D)
    for k ∈ 1:D
        for i in 1:D
            for j in 1:D
                TrAS[k] += g[i,j]*A[i,k,j]
            end
        end
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
                E = E + abs(A[i,j,k] - A[j,i,k]) + abs(A[i,j,k] - A[k,j,i])
            end
        end
    end
    return E
end

function get_A(FaA)
    D = length(FaA[:,1])
    A = zeros(D,D,D)
    for i in 1:D
        A[:,:,i] = vector2symmat(FaA[i,D+1:end])
    end
    return A[:,:,:]
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
    return (χ₀*μ^2 + χ₁*ν + χ₂*(τ₂ - 4*D1*𝚶*μ^(-1) + 2*D1*D2*ν^2*μ^(-2)))*μ^(-D1) 
    # the function seems to be correct: tr(g*T)=0 and also A is fully symmetric after symmetrization 
    # on first two indices
    #return 𝚶
end

function Φ_v(ζ,p)
    #ζ_v, ζ_t = vector_unpack(ζ)
    get_dim(ζ)
    ζ_t = vector2symmat(ζ[D+1:end])
    return Φ_vt(ζ[1:D], ζ_t, p)
    #return Φ_vt(ζ_v, ζ_t, p)
end

#Φ_vt(ζ_v,ζ_t, [1.,1.,1.])#, ζ_t, [1.;1;1])
#Φ_v(ζ, [1.,1.,1.])
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
    

function Φ_new(ζ,p)
    #χ = zeros(3)
    χs = p
    χ = Symbolics.scalarize(χs)
    #g = Symbolics.scalarize(gs)
    #D = length(g[:,1])
    μ = mu(Symbolics.scalarize(ζ))
    ν = nu(Symbolics.scalarize(ζ))
    𝚶 = omicron(Symbolics.scalarize(ζ))
    τ₂= tau2(Symbolics.scalarize(ζ))
    D1 = D//2+1
    D2 = D//2+2
    if D > 2
        return (χ[1]*μ^2 + χ[2]*ν + χ[3]*(τ₂ - 4*D1*𝚶*μ^(-1) + 2*D1*D2*ν^2*μ^(-2)))*μ^(-D1) 
    elseif D==2 
        return χ[1]*μ*log(μ) + (χ[2]*ν + χ[3]*(τ₂ - 4*D1*𝚶*μ^(-1) + 2*D1*D2*ν^2*μ^(-2)))*μ^(-D1)
    end
end

function get_dim(N)
    DD = (-3 + sqrt(9+8N))/2 
    if DD % 1 ≈ 0
        return D = Int64(DD)
    else
        error("In function vector_unpack: the vector dimension does not corresponds to the number of free entries of a symmetric matrix N = $N, D = $DD")
    end
end

function dd(i,j)
    return (i==j) ? 1 : 0.5
end

#"""
#trasnforms indices of a MxM symmetric matrix into a vector index
#"""
function l_ind(i::Int64,j::Int64,M::Int64)::Int64
    if i <= j 
    return Int64((i-1)*M - (i-1)*i÷2 + j)
    else
    return Int64((j-1)*M - (j-1)*j÷2 + i)
    end
end

function make_vector_TF!(v)
    D = get_dim(length(v))
    T = -v[D + l_ind(1,1,D)]
    for k in 2:D
        T = T + v[D + l_ind(k,k,D)]
    end
    v[D + l_ind(1,1,D)] = v[D + l_ind(1,1,D)] + T/D
    for k in 2:D
        v[D + l_ind(k,k,D)] = v[D + l_ind(k,k,D)] - T/D
    end
    return v[:]
end

    
function mu(v)
    D = get_dim(length(v))
    g = make_g(D)
    sum = 0
    for i in 1:D 
        for j in 1:D 
                #v_l g^{lk}v_k
                sum = sum + v[i]*g[i,j]*v[j]
        end
    end
    return sum
    #return v[1:D]' * g * v[1:D]
end

function omicron(v)
    D = get_dim(length(v))
    g = make_g(D)
    sum = 0
    for i in 1:D 
        for j in 1:D 
            for k in 1:D
                for l in 1:D
                    for n in 1:D
                        for m in 1:D
                            #v_l g^{lk} z_{kj} g^{ji} z_{in} g^{nm} v_m
                            sum = sum + v[l]*g[l,k]*v[D+l_ind(k,j,D)]*g[j,i]*v[D+l_ind(i,n,D)]*g[n,m]*v[m]
                        end
                    end
                end
            end
        end
    end
    return sum
end

function nu(v)
    D = get_dim(length(v))
    g = make_g(D)
    sum = 0
    for i in 1:D 
        for j in 1:D 
            for k in 1:D
                for l in 1:D
                    #v_{i}g^{ij}z{jk}g^{kl}v_{l}    
                    sum = sum + v[i]*g[i,j]*v[D+l_ind(j,k,D)]*g[k,l]*v[l]
                end
            end
        end
    end
    return sum
end
    

function tau2(v)
    D = get_dim(length(v))
    g = make_g(D)
    sum = 0
    for i in 1:D 
        for j in 1:D 
            for k in 1:D
                for l in 1:D
                    #g^{lj} * z_{lk} * g^{ki} * z{ij}
                    sum = sum + g[l,j]*v[D+l_ind(l,k,D)]*g[k,i]*v[D+l_ind(i,j,D)]
                end
            end
        end
    end
    return sum
end

function χaA(O, vars_a, vars_A; simplify=true)
    vars_a = map(Symbolics.value, vars_a)
    vars_A = map(Symbolics.value, vars_A)
    first_derivs = map(Symbolics.value, vec(Symbolics.jacobian([Symbolics.values(O)], vars_a, simplify=simplify)))
    n = length(vars_a) #D
    m = length(vars_A) #L
    H = Array{Num, 2}(undef,(n, m))
    fill!(H, 0)
    for j=1:n
        for i=1:m
            H[j, i] = Symbolics.expand_derivatives(Symbolics.Differential(vars_A[i])(first_derivs[j]))
        end
        #Take out the trace
        # compute the g-trace 
        T = - H[j,n+l_ind(1,1,n)]
        for k=2:n
            T = T + H[j,n+l_ind(k,k,n)]
        end
        #substract from vector 
        H[j,n+l_ind(1,1,n)] = H[j,n+l_ind(1,1,n)] + T/n
        for k=2:n
            H[j,n+l_ind(k,k,n)] = H[j,n+l_ind(k,k,n)] - T/n
        end
        for k=1:n
            for l=1:k
                H[j,n+l_ind(k,l,n)] = H[j,n+l_ind(k,l,n)]*dd(k,l)
            end
        end
    end
    H
end

function χ0AB(O, vars_A; simplify=true)
    #vars_0 = map(Symbolics.value, vars_0)
    vars_A = map(Symbolics.value, vars_A)
    first_derivs = map(Symbolics.value, vec(Symbolics.jacobian([Symbolics.values(O)], vars_A, simplify=simplify)))
    m = length(vars_A) #L
    n = get_dim(m)
    H = Array{Num, 2}(undef,(m, m))
    H0 = Array{Num, 2}(undef,(m, m))
    fill!(H, 0)
    #take out the trace 
    T = - first_derivs[n+l_ind(1,1,n)]
    for k=2:n
        T = T + first_derivs[n+l_ind(k,k,n)]
    end
    first_derivs[n+l_ind(1,1,n)] = first_derivs[n+l_ind(1,1,n)] + T/n
    for k=2:n
        first_derivs[n+l_ind(k,k,n)] = first_derivs[n+l_ind(k,k,n)] - T/n
    end
    #take second derivative
    for j=1:m
        for i=1:m
            H[j, i] = Symbolics.expand_derivatives(Symbolics.Differential(vars_A[i])(first_derivs[j]))
        end
        #take out the trace of the second index
        T = - H[j,n+l_ind(1,1,n)]
        for k=2:n
            T = T + H[j,n+l_ind(k,k,n)]
        end
        H[j,n+l_ind(1,1,n)] = H[j,n+l_ind(1,1,n)] + T/n
        for k=2:n
            H[j,n+l_ind(k,k,n)] = H[j,n+l_ind(k,k,n)] - T/n
        end
    end
    

    for i=1:m
        for j=1:i
            H0[j, i] = H0[i,j] = Symbolics.expand_derivatives(Symbolics.Differential(vars_A[1])(H[j,i]))
        end
    end
    H0
end