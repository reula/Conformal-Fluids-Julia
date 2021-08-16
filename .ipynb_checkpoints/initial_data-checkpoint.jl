function create_initial_data(name::String, u_i, par)
    χ, U, M, dx, Euler, F, Jac = par
    
#    u_i=zeros(M*U)
    con_0 = view(reshape(u_i,(M,U)),:,1:U÷2)
    flu_0 = view(reshape(u_i,(M,U)),:,U÷2+1:U)

    if name == "speed_pulse"
    

        x0 = 0.0; x1 = 0.2 #x0 = 0.4; x1 = 0.6
        v0 = 0.0; δv = 0.5

        for i in 1:M
            x = dx*(i-1)
            if x > x0 && x < x1
                flu_0[i,2] = v0 + δv *(x - x0)^4 * (x - x1)^4 / (x1-x0)^8 * 250                   #Sz
                flu_0[i,1] = -6.0 # * sin(pi*(x - x0)/(x1-x0))^4 * sin(2*pi*(x - x0)/(x1-x0))     #By
            else
                flu_0[i,2] = v0
                flu_0[i,1] = -6. 
            end
    
        end
    f_to_c!(u_i, (χ, U, M, F)); # populate the conservative variables from the fluid ones

    elseif name == "small_pulse_to_the_right"

        e0 = 6.; δe = 0.1
        x0 = 0.; x1 = 0.2
        λ = 1. /sqrt(3.)
        for i in 1:M
            x = dx*(i-1)
            if x > x0 && x < x1
                con_0[i,1] = e0 + δe *(x - x0)^4 * (x - x1)^4 / (x1-x0)^8 * 250                   #Sz
                con_0[i,2] = λ*(con_0[i,1] - e0)
            else
                con_0[i,1] = e0
                con_0[i,2] = 0.
            end
                flu_0[i,1] = -1.
        end
        par_inv = (χ, 1e-8, 10000, U, M, F, Jac)
        return u_i[:] = c_to_f!(u_i,par_inv);
        
    elseif (name == "big_pulse_to_the_right") || (name == "big_pulse_negative_I")
    
        e0 = 6.; δe = 4.0
        x0 = 0.; x1 = 0.2
        λ = 1. /sqrt(3.)
        for i in 1:M
            x = dx*(i-1)
            if x > x0 && x < x1
                con_0[i,1] = e0 + δe *(x - x0)^4 * (x - x1)^4 / (x1-x0)^8 * 250                   #Sz
                con_0[i,2] = λ*(con_0[i,1] - e0)
            else
                con_0[i,1] = e0
                con_0[i,2] = 0.
            end
                flu_0[i,1] = -1. #-sqrt(con_0[i,1])/6
        end
    
        if !Euler
            println("Enter intermediate step")
            χ_int = [χ[1]; 0.0; χ[3]] #primero nos aproximamos a la solución con χ₁=0
            par_inv = (χ_int, 1e-8, 10, U, M, F, Jac)
            u_i = c_to_f!(u_i,par_inv);
            println("After intermediate step")
            #luego buscamos la buena con la semilla anterior
            par_inv = (χ, 1e-8, 10, U, M, F, Jac)
            u_i = c_to_f!(u_i,par_inv);
            println("After whole step")
        #check!
            if maximum(abs, f_to_c!(u_i,(χ,U,M,F)) - u_i) > 0.001
            error("Not converging")
            else
                return u_i
            end
            
        else
            return u_i = c_to_f_direct!(u_i,(U,M))
        end #!Euler
        
    elseif name == "constant_fields" 
    
        e0 = 6.; δe = 0.0
        s0 = 0.; δs = 0.0
        x0 = 0.; x1 = 0.2
        c1 = 0.1
        c2 = 0.2
        c3 = 0.3
        λ = 1. /sqrt(3.)
        for i in 1:M
            x = dx*(i-1)
            if x > x0 && x < x1
                con_0[i,1] = e0 + δe *(x - x0)^4 * (x - x1)^4 / (x1-x0)^8 * 250            #Sz
                con_0[i,2] = λ*(con_0[i,1] - e0) + s0
            else
                con_0[i,1] = e0
                con_0[i,2] = s0
            end
                flu_0[i,1] = -sqrt(con_0[i,1])/6
                con_0[i,3] = c1
                con_0[i,4] = c2
                con_0[i,5] = c3
        end
        par_inv = (χ, 1e-8, 100, U, M, F, Jac)
        return u_i = c_to_f!(u_i,par_inv);

                
            elseif name == "square"

        mu0 = -6.; δmu = 1.0
        x0 = 0.4; x1 = 0.6
        x10 = 0.0
        x20 = 0.0
        x30 = 0.0
        #λ = 1. /sqrt(3.)
        for i in 1:M
            x = dx*(i-1)
            if x > x0 && x < x1
                flu_0[i,1] = mu0 + δmu
                flu_0[i,3] = x10 *(x - x0)^4 * (x - x1)^4 / (x1-x0)^8 * 250                   #Sz
                flu_0[i,4] = x20 *(x - x0)^4 * (x - x1)^4 / (x1-x0)^8 * 250 
                flu_0[i,5] = x30 *(x - x0)^4 * (x - x1)^4 / (x1-x0)^8 * 250 
            else
            flu_0[i,1] = mu0
            end
        end
    
        return u_i = f_to_c!(u_i, (χ, U, M, F));
        
    elseif name == "only_diss"

        e0 = 6.; δe = 0.0
        x0 = 0.; x1 = 0.2
        x10 = 0.1
        x20 = 0.1
        x30 = 0.1
        #λ = 1. /sqrt(3.)
        for i in 1:M
            x = dx*(i-1)
            if x > x0 && x < x1
                con_0[i,1] = e0
                con_0[i,3] = x10 *(x - x0)^4 * (x - x1)^4 / (x1-x0)^8 * 250                   #Sz
                con_0[i,4] = x20 *(x - x0)^4 * (x - x1)^4 / (x1-x0)^8 * 250 
                con_0[i,5] = x30 *(x - x0)^4 * (x - x1)^4 / (x1-x0)^8 * 250 
            end
            con_0[i,1] = e0
            flu_0[i,1] = -1.
        end
        return u_i = c_to_f!(u_i,par_inv);
       
    else
        println("unknown data_name")
    end
end