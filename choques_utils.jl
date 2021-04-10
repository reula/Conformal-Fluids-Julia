function mysign_zero(a::AbstractFloat)
    if (a > 0.) 
        return 1.0
	elseif (a < 0.) 
        return -1.0
	else 
        return 0.0
    end
end
#mysign_zero(-0.)

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

function mp5(du::Array{Float64,1}, u::Array{Float64,1}, par, j) # j is the grid position
   
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

        dv[j,i] = - h_1 * (H_p[i]-H_m[i])
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

    return du[j,:]  = - h_1 * (H_p[:]-H_m[:])
end
#mp5(u,du,par,2)
