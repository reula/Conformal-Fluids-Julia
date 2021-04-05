import Pkg; 
Pkg.activate("Conf_Fluids")
Pkg.add("Symbolics")
using Symbolics
using StaticArrays

function F(flu,con,χ)
    #the use of a dictionary can go away in any case.
    f = Dict(:μ => 1, :v => 2, :x1 => 3, :x2 => 4,:x3 => 5) #para acordarse de que campo asignamos a cada componente del vector
    c = Dict(:e => 1, :s => 2, :c1 => 3, :c2 => 4,:c3 => 5)

    μ = flu[f[:μ]]  # esto es -μ
    #μ = view(flu,1)
    T = (μ)^(-1//2) # use μ positive, so I changed μ -> -μ
    v = flu[f[:v]]
    χ₀ = χ[1]
    χ₁ = χ[2]
    χ₂ = χ[3]
    x1 = flu[f[:x1]]
    x2 = flu[f[:x2]]
    x3 = flu[f[:x3]]
    γ = (1 - v^2)^(-1//2)
    τ = -2χ₁ * x3 * T / (γ*μ^3) - 24χ₂*(1//2*(1-v^2)x3^2 + 14//3 * x2^2 + 7/5*(1-v^2)x1*x3)/μ^5
    ρ = -6χ₀ / μ^2 - 6χ₁*x1/(γ * μ^4 * T) - 42χ₂*(6//5 *x1^2 + 10γ^2*x2^2 + 3//2*(v^2-1)^2*x3^2)/(μ^5 * γ^2)
    Q = -10χ₀ * x2 * T / μ^3 - 168χ₂ * x2 * (x1 - (v^2 - 1)x3)/(γ * μ^5)

    A = (-12 * χ₂/μ^4)*[3(2γ^2 - 1)  3v*(6γ^2 - 1)    3v^2 ;
                v*(6γ^2 - 1)  (6γ^2*(1 + 2v^2) - 1)  v*(v^2 + 2) ;
                (6γ^2 * v^2 + 1)  v*(6γ^2*(2 + v^2) - 1)  (2v^2 + 1)]

    e = 4//3*ρ*(γ^2 - 1//4) + 2v*γ*Q + v^2*τ
    s = 4//3*ρ*γ^2*v + (v^2+1//1)γ*Q + v*τ;

    Y1 = -3χ₁*γ/T/μ^3 *(2γ^2 - 1//1)
    Y2 = -3χ₁*γ/T/μ^3 *v*(6γ^2 - 1//1)
    Y3 = -3χ₁*γ/T/μ^3 *(6γ^2*v^2 + 1//1)
    
    return  [- con[c[:e]] + 4//3* ρ*(γ^2 - 1//4) + 2//1* Q*γ*v + τ*v^2;
             - con[c[:s]] + 4//3* ρ*γ^2*v + Q*γ*(v^2 + 1//1) + τ*v;
             - con[c[:c1]] + A[1,1]*x1 + A[1,2]*x2 + A[1,3]*x3 + Y1;
             - con[c[:c2]] + A[2,1]*x1 + A[2,2]*x2 + A[2,3]*x3 + Y2;
             - con[c[:c3]] + A[3,1]*x1 + A[3,2]*x2 + A[3,3]*x3 + Y3]
end

N=5

@variables f[1:N], c[1:N], p[1:3]

JS = Symbolics.jacobian(F(f,c,p),f);
J_exp = Symbolics.build_function(JS, f, c, p);
Jac = eval(J_exp[1]);

function NR_step!(F, Jac, u0, y, p)
    u = u0 - Jac(u0,y, p) \ F(u0,y,p)
end

"""
    c_to_f!(flu, con, p)

Compute the fluid variables `flu` from the conservative `con` ones.
It inverts the function `F(flu,con,χ) = con(flu,χ) - con` using 
Newto-Raphson method.

As parameters have: 
    χ = Function parameters, 
    tol = error tolerance for the inversion, 
    iter_max = maximum iteration number for the NR method, 
    N = number of variables (5 here), 
    M = dimension of vector of variables

# Examples
```julia-repl
julia> 
tol = 10^(-16)
iter_max = 40
χ = [-1.0; -1.0; -5.0]
M=1
N=5
p = (χ, tol, iter_max, N, M)
flu = [-10.0; 0.2; 1.1; 1.5; 1.2]
con = [0.1366314976448222; 0.07009306769467444; 0.06115332989597844; 0.07178418128379448; 0.04927907295689818]
c_to_f!(flu + 0.1* rand(N), con, p) - con
5×1 Array{Float64,2}:
  2.1953994178147695e-11
 -1.393440918207034e-12
 -4.0323300254385686e-13
 -1.3405498933138915e-11
  1.2895684520231043e-11
```
"""
function c_to_f!(flu, con, p)
    χ, tol, iter_max, N, M = p
    for j ∈ 1:M
        flu[1,j] = -flu[1,j]
        iter = 1
        while F(flu[:,j],con[:,j], χ)'*F(flu[:,j],con[:,j], χ) > tol && iter < iter_max
            flu[:,j] = NR_step!(F, Jac, flu[:,j], con[:,j], χ)
            iter = iter + 1
            if iter == iter_max 
            println("iter_max reached j = $j")
            end
        end
    #println(iter)
        flu[1,j] = -flu[1,j]
    end
    return flu[:,:]
end



"""
    f_to_c!(flu, con, p)

Compute the conservative `con` variables for the fluid ones `flu`.
It uses the function `F(flu,con,χ) = con(flu,χ) - con` with `con = 0`.

As parameters have: 
    χ = Function parameters
    N = number of variables (5 here), 
    M = dimension of vector of variables

# Examples
```julia-repl
julia> 
χ = [-1.0; -1.0; -5.0]
M=1
N=5
p = (χ, N, M)
flu = [-10.0; 0.2; 1.1; 1.5; 1.2]
con = [0.1366314976448222; 0.07009306769467444; 0.06115332989597844; 0.07178418128379448; 0.04927907295689818]
f_to_c!(flu,con+rand(N),p) - con
5×1 Array{Float64,2}:
  0.0
 -1.3877787807814457e-17
  0.0
  0.0
  0.0
```
"""
function f_to_c!(flu, con, p)
    χ, N, M = p
    y = zeros(N)
    for j ∈ 1:M
        flu[1,j] = -flu[1,j]
        con[:,j] = F(flu[:,j],y, χ)
    end
    return con[:,:]
end

