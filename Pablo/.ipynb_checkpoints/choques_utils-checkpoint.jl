@inline function mysign_zero(a)
    return (1.0.*(a .> 0.0) + (-1.0).* (a .< 0.0))
end

@inline function minmod(a, b)
    return 0.5*(mysign_zero(a)+mysign_zero(b))*min(abs(a), abs(b))
end

@inline function minmod(a, b, c)
    sgnbc = (mysign_zero(b)+mysign_zero(c)) #this is 2 if both are positive, -2 if both are negative, 0 otherwise
    sgnac = (mysign_zero(a)+mysign_zero(c)) #this is 2 if both are positive, -2 if both are negative, 0 otherwise
    
    return 0.25*sgnbc*abs(sgnac)*min(abs(a), abs(b),abs(c))
end

@inline function minmod(a,b,c,d)
    return 0.125*(mysign_zero(a)+mysign_zero(b))*(abs((mysign_zero(a)+mysign_zero(c))*(mysign_zero(a)+mysign_zero(d))))*min(abs(a),abs(b),abs(c), abs(d))
end

function MM3(a::AbstractFloat, b::AbstractFloat, c::AbstractFloat, weight::AbstractFloat) #(2*D0,Dp,Dm)
    
    weight = weight * 2.
  
    if (abs(a) <= (weight*abs(b))) 
        return (abs(a) <= (weight*abs(c))) ? abs(a)*.5 : abs(c) 
    else 
        return (abs(b) <= abs(c)) ? abs(b) : abs(c)
    end 
end
#MM3(4.,-1.,4.,2.5)
function MM3N(a::AbstractFloat, b::AbstractFloat, c::AbstractFloat) #(2*D0,Dp,Dm)
    if (abs(a) <= abs(b))
        return (abs(a) <= abs(c)) ? abs(a) : abs(c) 
    else 
        return (abs(b) <= abs(c)) ? abs(b) : abs(c)
    end 
end
#MM3(4.,-1.,4.,2.5)


function DMM(a::AbstractFloat, b::AbstractFloat)
    return 0.5 * (mysign_zero(a) + mysign_zero(b)) * minimum([abs(a),abs(b)])
end
#DMM(2.,2.)

function DM4(a::AbstractFloat, b::AbstractFloat, c::AbstractFloat, d::AbstractFloat)
    return 0.125 * (mysign_zero(a) + mysign_zero(b)) * abs((mysign_zero(a) + mysign_zero(c)) * (mysign_zero(a) + mysign_zero(d))) * minimum([abs(a),abs(b),abs(c),abs(d)])
end
#DM4(1.,2.,3.,0.)

function MP5(F0::AbstractFloat, F1::AbstractFloat, F2::AbstractFloat, F3::AbstractFloat, F4::AbstractFloat)::AbstractFloat
    b1 = 0.0166666666667
    b2 = 1.3333333333333
    α = 4.
    ϵ = 1.e-10
    #remember that we are off 3 places, i-2 = 0 
    vor = b1 * (2. * F0 - 13. * F1 + 47. * F2 + 27. * F3 - 3. * F4)
                vmp = F2 + DMM(F3 - F2, α * (F2 - F1));
    if ((vor - F2) * (vor - vmp) <= ϵ) 
        return vor;
    else 
        djm1 = F0 - 2. * F1 + F2;
        dj = F1 - 2. * F2 + F3;
        djp1 = F2 - 2. * F3 + F4;
        dm4jph = DM4(4. * dj - djp1, 4 * djp1 - dj, dj, djp1);
        dm4jmh = DM4(4. * dj - djm1, 4 * djm1 - dj, dj, djm1);
        vul = F2 + α * (F2 - F1);
        vav = 0.5 * (F2 + F3);
        vmd = vav - 0.5 * dm4jph;
        vlc = F2 + 0.5 * (F2 - F1) + b2 * dm4jmh;
        vmin = maximum([minimum([F2,F3,vmd]), minimum([F2, vul, vlc])]);
        vmax = minimum([maximum([F2,F3,vmd]), maximum([F2, vul, vlc])]);
        return vor + DMM(vmin - vor, vmax - vor);
    end
end

function mp5(du::Array{Float64,1}, u::Array{Float64,1}, par, j, t) # j is the grid position
   
    h_1, U, N, par_flux, par_source, Fx, Speed_max, Source = par
    
    v = reshape(u,(N,U))
    dv = reshape(du,(N,U))
    
    F_RM = zeros(U)
    F_RP = zeros(U)
    F_LM = zeros(U)
    F_LP = zeros(U)
    H_p = zeros(U)
    H_m = zeros(U)
    
    u_l = v[j,:]
    u_l_p = v[mod1((j+1), N),:]
    u_l_m = v[mod1((j + N -1), N),:]
    u_l_p2 = v[mod1((j+2), N),:]
    u_l_m2 = v[mod1((j + N -2), N),:]
    u_l_p3 = v[mod1((j+3), N),:]
    u_l_m3 = v[mod1((j + N -3), N),:]
    
    S_MAX = maximum([Speed_max(u_l_p3, par_flux), Speed_max(u_l_m3, par_flux), 
            Speed_max(u_l_p2, par_flux), Speed_max(u_l_m2, par_flux), Speed_max(u_l_p, par_flux), 
            Speed_max(u_l_m, par_flux), Speed_max(u_l, par_flux)])
    
    F_Pp3 = 0.5 * (Fx(u_l_p3, par_flux) + S_MAX * u_l_p3)
    F_Mp3 = 0.5 * (Fx(u_l_p3, par_flux) - S_MAX * u_l_p3)
    F_Pp2 = 0.5 * (Fx(u_l_p2, par_flux) + S_MAX * u_l_p2)
    F_Mp2 = 0.5 * (Fx(u_l_p2, par_flux) - S_MAX * u_l_p2)
    F_Pp  = 0.5 * (Fx(u_l_p,  par_flux) + S_MAX * u_l_p)
    F_Mp  = 0.5 * (Fx(u_l_p,  par_flux) - S_MAX * u_l_p)
    F_P   = 0.5 * (Fx(u_l,    par_flux) + S_MAX * u_l)
    F_M   = 0.5 * (Fx(u_l,    par_flux) - S_MAX * u_l)
    F_Pm  = 0.5 * (Fx(u_l_m,  par_flux) + S_MAX * u_l_m)
    F_Mm  = 0.5 * (Fx(u_l_m,  par_flux) - S_MAX * u_l_m)
    F_Pm2 = 0.5 * (Fx(u_l_m2, par_flux) + S_MAX * u_l_m2)
    F_Mm2 = 0.5 * (Fx(u_l_m2, par_flux) - S_MAX * u_l_m2)
    F_Pm3 = 0.5 * (Fx(u_l_m3, par_flux) + S_MAX * u_l_m3)
    F_Mm3 = 0.5 * (Fx(u_l_m3, par_flux) - S_MAX * u_l_m3)
    
    for i in 1:U
        F_RM[i] = MP5(F_Mp2[i], F_Mp[i],  F_M[i],  F_Mm[i], F_Mm2[i])
        F_LM[i] = MP5(F_Pm3[i], F_Pm2[i], F_Pm[i], F_P[i],  F_Pp[i])
        F_LP[i] = MP5(F_Pm2[i], F_Pm[i],  F_P[i],  F_Pp[i], F_Pp2[i])
        F_RP[i] = MP5(F_Mp3[i], F_Mp2[i], F_Mp[i], F_M[i],  F_Mm[i])
        
        H_p[i] = F_LP[i] + F_RP[i]
        H_m[i] = F_LM[i] + F_RM[i]

        dv[j,i] = - h_1 * (H_p[i]-H_m[i]) + Source(u,t,par_source)[i]
    end
    return du[j,:]
end
#mp5(u,du,par,2)


function kt(du::Array{Float64,1}, u::Array{Float64,1}, par, j) # j is the grid position
    par_flux, h_1, U, N, Fx, Speed_max = par 
    theta = 1.5
    v = reshape(u,(N,U))
    dv = reshape(du,(N,U))
    
    v_x = zeros(U)
    v_xp = zeros(U)
    v_xm = zeros(U)
    
    u_l = v[j,:]
    u_l_p = v[mod1((j+1), N),:]
    u_l_m = v[mod1((j + N -1), N),:]
    u_l_p2 = v[mod1((j+2), N),:]
    u_l_m2 = v[mod1((j + N -2), N),:]
    
    Dp = u_l_p - u_l
    Dpp = u_l_p2 - u_l_p
    Dm = u_l - u_l_m
    Dmm = u_l_m - u_l_m2
    
    for i in 1:U
        v_x[i]  = 0.5*h_1*(mysign_zero(Dp[i])+mysign_zero(Dm[i]))*MM3N(0.5(Dp[i]+Dm[i]),Dp[i]*theta,Dm[i]*theta);
        v_xp[i] = 0.5*h_1*(mysign_zero(Dpp[i])+mysign_zero(Dp[i]))*MM3N(0.5(Dpp[i]+Dp[i]),Dpp[i]*theta,Dp[i]*theta);
        v_xm[i] = 0.5*h_1*(mysign_zero(Dm[i])+mysign_zero(Dmm[i]))*MM3N(0.5(Dm[i]+Dmm[i]),Dm[i]*theta,Dmm[i]*theta);
    end
    
    u_pp = u_l_p - 0.5 * v_xp / h_1
    u_pm = u_l   + 0.5 * v_x / h_1
    u_mp = u_l   - 0.5 * v_x / h_1
    u_mm = u_l_m + 0.5 * v_xm / h_1
    
    a_p = maximum([Speed_max(u_pp, par_flux), Speed_max(u_pm, par_flux)])
    a_m = maximum([Speed_max(u_mm, par_flux), Speed_max(u_mp, par_flux)])
        
    H_p = 0.5 * (Fx(u_pp, par_flux) + Fx(u_pm, par_flux)) - 0.5 * a_p * (u_pp - u_pm)
    H_m = 0.5 * (Fx(u_mm, par_flux) + Fx(u_mp, par_flux)) - 0.5 * a_m * (u_mp - u_mm)

    return du[j,:]  = - h_1 * (H_p[:]-H_m[:]) + Source(u,t,par_source)[i]
end
#mp5(u,du,par,2)


#MP5 Reconstruction
function MP5reconstruction(Vjmm, Vjm, Vj, Vjp, Vjpp)
    B1 = 0.0166666666666666667  #1/60
    B2 = 1.3333333333333333333  #4/3
    eps = 1e-10
    ALPHA = 4.0
    #=Vjmm = V[1]
    Vjm = V[2]
    Vj = V[3]
    Vjp = V[4]
    Vjpp = V[5]=#
    Vor = B1*(2.0*Vjmm - 13.0*Vjm + 47.0*Vj + 27*Vjp - 3.0*Vjpp) #=This is the original interpolation.
                                                                   All that follows is the application of 
                                                                   limiters to treat shocks=#
    Vmp = Vj + minmod(Vjp-Vj, ALPHA*(Vj-Vjm))  #mp = monotonicity preserving. It's the median between v_j, v_(j+1)
                                              #and an upper limit v^UL = v_j+ALPHA(v_j-v_(j-1))
    if ((Vor-Vj)*(Vor-Vmp)) < eps             #this condition is equivalent to asking vl in [vj, v^{MP}]
        Vl = Vor #vl = v^{L}_{j+1/2}
    else
        djm1 = Vjmm - 2.0*Vjm + Vj
        dj = Vjm - 2*Vj + Vjp
        djp1 = Vj - 2.0*Vjp + Vjpp
        dm4jph = minmod(4*dj - djp1, 4*djp1-dj, dj, djp1)  #ph = plus half (+1/2)
        dm4jmh = minmod(4*dj - djm1, 4*djm1-dj, dj, djm1)  #mh = minus half (-1/2)
        #d^{M4}_{j+1/2} = \minmod(4d_{j}-d_{j+1},4d_{j+1}-d_{j}, d_{j}, d_{j+1})
        Vul = Vj + ALPHA*(Vj - Vjm)   #upper limit
        Vav = 0.5*(Vj + Vjp)          #average
        Vmd = Vav - 0.5*dm4jph        #Vmedian
        Vlc = Vj + 0.5*(Vj-Vjm) + B2*dm4jmh
        Vmin = max(min(Vj, Vjp, Vmd), min(Vj, Vul, Vlc));
        Vmax = min(max(Vj, Vjp, Vmd), max(Vj, Vul, Vlc));
        Vl = Vor + minmod(Vmin-Vor, Vmax-Vor) #this places Vor between Vmin and Vmax
    end
    return Vl
end 


function MP5reconstruction!(Vl, Vjmm, Vjm, Vj, Vjp, Vjpp, N_Fields)
    B1 = 0.0166666666666666667  #1/60
    B2 = 1.3333333333333333333  #4/3
    eps = 1e-10
    ALPHA = 4.0
    #=Vjmm = V[1]
    Vjm = V[2]
    Vj = V[3]
    Vjp = V[4]
    Vjpp = V[5]=#
    for i in 1:N_Fields
        Vor = B1*(2.0*Vjmm[i] - 13.0*Vjm[i] + 47.0*Vj[i] + 27*Vjp[i] - 3.0*Vjpp[i]) #=This is the original interpolation.
                                                                       All that follows is the application of 
                                                                       limiters to treat shocks=#
        Vmp = Vj[i] + minmod(Vjp[i]-Vj[i], ALPHA*(Vj[i]-Vjm[i]))  #mp = monotonicity preserving. It's the median between v_j, v_(j+1)
                                                  #and an upper limit v^UL = v_j+ALPHA(v_j-v_(j-1))
        if ((Vor-Vj[i])*(Vor-Vmp)) < eps             #this condition is equivalent to asking vl in [vj, v^{MP}]
            Vl[i] = Vor #vl = v^{L}_{j+1/2}
        else
            djm1 = Vjmm[i] - 2.0*Vjm[i] + Vj[i]
            dj = Vjm[i] - 2*Vj[i] + Vjp[i]
            djp1 = Vj[i] - 2.0*Vjp[i] + Vjpp[i]
            dm4jph = minmod(4*dj - djp1, 4*djp1-dj, dj, djp1)  #ph = plus half (+1/2)
            dm4jmh = minmod(4*dj - djm1, 4*djm1-dj, dj, djm1)  #mh = minus half (-1/2)
            #d^{M4}_{j+1/2} = \minmod(4d_{j}-d_{j+1},4d_{j+1}-d_{j}, d_{j}, d_{j+1})
            Vul = Vj[i] + ALPHA*(Vj[i] - Vjm[i])   #upper limit
            Vav = 0.5*(Vj[i] + Vjp[i])          #average
            Vmd = Vav - 0.5*dm4jph        #Vmedian
            Vlc = Vj[i] + 0.5*(Vj[i]-Vjm[i]) + B2*dm4jmh
            Vmin = max(min(Vj[i], Vjp[i], Vmd), min(Vj[i], Vul, Vlc));
            Vmax = min(max(Vj[i], Vjp[i], Vmd), max(Vj[i], Vul, Vlc));
            Vl[i] = Vor + minmod(Vmin-Vor, Vmax-Vor) #this places Vor between Vmin and Vmax
        end
    end
end 



function mp5!(dv, v, par, t) # j is the grid position
    #asumimos u unidimensional por ahora
    h, N_Fields, N, par_flux, par_source, Fx!, Speed_max, Source!, par_mem = par
    #h_1, U, N, χ, par_flux, Flux, Speed_max, Source, par_mem = par_mp5 
    F_Mm3, F_Mm2, F_Mm1, F_M, F_Mp1, F_Mp2, F_Mp3, F_Pm3, F_Pm2, F_Pm1, F_P, F_Pp1, F_Pp2, F_Pp3, F_LP,   F_LM, F_RP, F_RM, H_m, H_p, sourcevec = par_mem
    
    fields = reshape(v,(N,N_Fields))
    dfields = reshape(dv,(N,N_Fields))
    #nota: f minuscula o u se usa para hablar de campos, F mayúscula para hablar de Flujos.
    
    for idx in 1:N
        #first we defined shifted indices
        idxm3 = mod(((idx-3) - 1),N) + 1
        idxm2 = mod(((idx-2) - 1),N) + 1
        idxm1 = mod(((idx-1) - 1),N) + 1
        idxp1 = mod(((idx+1) - 1),N) + 1
        idxp2 = mod(((idx+2) - 1),N) + 1
        idxp3 = mod(((idx+3) - 1),N) + 1
        
    
        um3 = @view fields[idxm3,:]
        um2 = @view fields[idxm2,:]
        um1 = @view fields[idxm1,:]
        u   = @view fields[idx,:]
        up1 = @view fields[idxp1,:]
        up2 = @view fields[idxp2,:]
        up3 = @view fields[idxp3,:]
        
        S_MAX = max(Speed_max(up3, par_flux), Speed_max(um3, par_flux), 
            Speed_max(up2, par_flux), Speed_max(um2, par_flux), Speed_max(up1, par_flux), 
            Speed_max(um1, par_flux), Speed_max(u, par_flux)) #maximum speed
        
        Fx!(F_Pm3, um3, par_flux)
        Fx!(F_Pm2, um2, par_flux)
        Fx!(F_Pm1, um1, par_flux)
        Fx!(F_P, u, par_flux)
        Fx!(F_Pp1, up1, par_flux)
        Fx!(F_Pp2, up2, par_flux)
        Fx!(F_Pp3, up3, par_flux)
        
        
        @. F_Mm3 = 0.5 * (F_Pm3 - S_MAX * um3)
        @. F_Mm2 = 0.5 * (F_Pm2 - S_MAX * um2)
        @. F_Mm1 = 0.5 * (F_Pm1 - S_MAX * um1)
        @. F_M   = 0.5 * (F_P   - S_MAX * u)
        @. F_Mp1 = 0.5 * (F_Pp1 - S_MAX * up1)
        @. F_Mp2 = 0.5 * (F_Pp2 - S_MAX * up2)
        @. F_Mp3 = 0.5 * (F_Pp3 - S_MAX * up3)
        @. F_Pm3 = 0.5 * (F_Pm3 + S_MAX * um3)
        @. F_Pm2 = 0.5 * (F_Pm2 + S_MAX * um2)
        @. F_Pm1 = 0.5 * (F_Pm1 + S_MAX * um1)
        @. F_P   = 0.5 * (F_P   + S_MAX * u)
        @. F_Pp1 = 0.5 * (F_Pp1 + S_MAX * up1)
        @. F_Pp2 = 0.5 * (F_Pp2 + S_MAX * up2)
        @. F_Pp3 = 0.5 * (F_Pp3 + S_MAX * up3)
    
        MP5reconstruction!(F_RM, F_Mp2, F_Mp1,  F_M,  F_Mm1, F_Mm2, N_Fields)
        MP5reconstruction!(F_LM, F_Pm3, F_Pm2, F_Pm1, F_P,  F_Pp1, N_Fields)
        MP5reconstruction!(F_LP, F_Pm2, F_Pm1,  F_P,  F_Pp1, F_Pp2, N_Fields)
        MP5reconstruction!(F_RP, F_Mp3, F_Mp2, F_Mp1, F_M,  F_Mm1, N_Fields)
        
        @. H_p = F_LP + F_RP
        @. H_m = F_LM + F_RM
        
        Source!(sourcevec, u, t, par_source)
        @. dfields[idx, :] = -h*(H_p - H_m) + sourcevec#+ Source(u,t,par_source)
        
    end
    
end