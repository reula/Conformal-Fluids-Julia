function create_initial_data(name::String, u_i, par)
    χ, U, M, Euler = par
    
#    u_i=zeros(M*U)
    con_0 = view(reshape(u_i,(M,U)),:,1:U÷2)
    flu_0 = view(reshape(u_i,(M,U)),:,U÷2+1:U)

    if name == "speed_pulse"
    

        x0 = 0.0; x1 = 0.2 #x0 = 0.4; x1 = 0.6
        v0 = 0.0; δv = 0.1 

        for i in 1:M
            x[i] = dx*(i-1)
            if x[i] > x0 && x[i] < x1
                flu_0[i,2] = v0 + δv *(x[i] - x0)^4 * (x[i] - x1)^4 / (x1-x0)^8 * 250                   #Sz
                flu_0[i,1] = -6.0 # * sin(pi*(x[i] - x0)/(x1-x0))^4 * sin(2*pi*(x[i] - x0)/(x1-x0))     #By
            else
                flu_0[i,2] = v0
            end
    
            flu_0[i,1] = -1.   #
        #flu_0[i,2] = flu_0[i,2] + 0.1  # v
    
        end
    f_to_c!(u_i, (χ, U, M)); # populate the conservative variables from the fluid ones

    elseif name == "small_pulse_to_the_right"

        e0 = 6.; δe = 0.1
        x0 = 0.; x1 = 0.2
        λ = 1. /sqrt(3.)
        for i in 1:M
            x[i] = dx*(i-1)
            if x[i] > x0 && x[i] < x1
                con_0[i,1] = e0 + δe *(x[i] - x0)^4 * (x[i] - x1)^4 / (x1-x0)^8 * 250                   #Sz
                con_0[i,2] = λ*(con_0[i,1] - e0)
            else
                con_0[i,1] = e0
                con_0[i,2] = 0.
            end
                flu_0[i,1] = -1.
        end
        par_inv = (χ, tol, 10000, U, M)
        return u_i[:] = c_to_f!(u_i,par_inv);
        
    elseif (name == "big_pulse_to_the_right") || (name == "big_pulse_negative_I")
    
        e0 = 6.; δe = 4.0
        x0 = 0.; x1 = 0.2
        λ = 1. /sqrt(3.)
        for i in 1:M
            x[i] = dx*(i-1)
            if x[i] > x0 && x[i] < x1
                con_0[i,1] = e0 + δe *(x[i] - x0)^4 * (x[i] - x1)^4 / (x1-x0)^8 * 250                   #Sz
                con_0[i,2] = λ*(con_0[i,1] - e0)
            else
                con_0[i,1] = e0
                con_0[i,2] = 0.
            end
                flu_0[i,1] = -sqrt(con_0[i,1])/6
        end
    
        if !Euler
            χ_int = [- 1.0; -0.; - 10.0] #primero nos aproximamos a la solución con χ₁=0
            par_inv = (χ, tol, 1000, U, M)
            u_i = c_to_f!(u_i,par_inv);
            χ = [- 1.0; -0.5; - 10.0] #luego buscamos la buena con la semilla anterior
            par_inv = (χ, tol, 1000, U, M)
            u_i = c_to_f!(u_i,par_inv);
        #check!
            if max(abs, f_to_c!(u_i,(χ,U,M)) - u_i) > 0.001
            error("Not converging")
            else
                return u_i
            end
            
        else
            return u_i = c_to_f!(u_i,par_inv);
        end #!Euler

    elseif name == "only_diss"

        e0 = 6.; δe = 0.0
        x0 = 0.; x1 = 0.2
        x10 = 0.1
        x20 = 0.1
        x30 = 0.1
        #λ = 1. /sqrt(3.)
        for i in 1:M
            x[i] = dx*(i-1)
            if x[i] > x0 && x[i] < x1
                con_0[i,1] = e0
                con_0[i,3] = x10 *(x[i] - x0)^4 * (x[i] - x1)^4 / (x1-x0)^8 * 250                   #Sz
                con_0[i,4] = x20 *(x[i] - x0)^4 * (x[i] - x1)^4 / (x1-x0)^8 * 250 
                con_0[i,5] = x30 *(x[i] - x0)^4 * (x[i] - x1)^4 / (x1-x0)^8 * 250 
            end
            con_0[i,1] = e0
            flu_0[i,1] = -1.
        end
        return u_i = c_to_f!(u_i,par_inv);
       
    else
        println("unknown data_name")
    end
end