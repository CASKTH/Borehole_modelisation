

@variables t x
f_c_d2(x, dx)::Int64 = ceil(Int, (x/dx))
# f_c_d2(x) = Base.ceil(x, Int64)
# f_c_d2(2, 10)
@register_symbolic f_c_d2(x)

f_c_d(x, i) = ceil(x) == i
@register_symbolic f_c_d(x,i)

# func(t) = 1/(1+exp(-(t-20000)/2000))*20

#     @register_symbolic func(t)

@connector function Pin(; name)
    sts = @variables T(t)=1.0 q(t)=1.0 [connect = Flow]
    ODESystem(Equation[], t, sts, []; name = name)
end


@connector function Pin_fluid(; name)
    sts = @variables T(t)=1.0 i(t,x)=1.0 [connect = Flow]
    ODESystem(Equation[], t, sts, []; name = name)
end
@component function Pin_port(; name)
    @named p = Pin()

    sts = @variables T(t)=1.0 i(t,x)=1.0  
    eqs = [T ~ p.T
           q ~ p.q]
    compose(ODESystem(eqs, t, sts, []; name = name), p)
end
@component function Ground(; name)
    @named g = Pin()
    eqs = [g.T ~ 0]
    compose(ODESystem(eqs, t, [], []; name = name), g)
end

@component function OnePort(; name)
    @named p = Pin()
    @named n = Pin()
    sts = @variables T(t)=1.0 q(t)=1.0  
    eqs = [T ~ p.T - n.T
           0 ~ p.q + n.q
           q ~ p.q]
    compose(ODESystem(eqs, t, sts, []; name = name), p, n)
end

# @component function OnePort_Var(var; name)
    
    
#     @named n = Pin(),
#     @named p = Pin(),
#     Pin_array = [(p_n = Pin(name= Symbol(string("p_",index))) for index=1:var)...
#     ]
    
#     # All pin are supposed to be independante because they have their voltage
#     sts = @variables T(t)=1.0 q(t)=1.0 
#     eqs = [T ~ p.T - n.T 
#            0 ~ p.q + n.q
#            q ~ p.q]
#     compose(ODESystem(eqs, t, sts, []; name = name), p, n )
# end

@component function Distributor(constant, depth; name)
    @named p = Pin()
    @named n = Pin()
    ps = @parameters q(t) constant=constant depth=depth
    sts = @variables T(t)=1.0 q(t)=1.0  
    eqs = [
            # T ~ i/(2.0)
            -q ~ p.q + n.q
            (q*depth)/constant ~ p.T - n.T
           
           ]
    compose(ODESystem(eqs, t, sts, ps; name = name), p, n)
end

# @component function Pump(Delta_Func; name)
#     @named oneport = OnePort()
#     @unpack T, q = oneport

#     eqs = [
#         T ~ Delta_Func(t)
        
#     ]
#     extend(ODESystem(eqs, t, [],[]; name = name), oneport)
# end

@component function Resistor(; name, R = 1.0)
    @named oneport = OnePort()
    @unpack T, q = oneport
    ps = @parameters R = R
    eqs = [
        T ~ q * R,
    ]
    extend(ODESystem(eqs, t, [], ps; name = name), oneport)
end

@component function Joint(; name)
    @named p = Pin()
    @named n = Pin()

    sts = @variables T(t)=1.0 q(t)=1.0  
    eqs = [
            # T ~ i/(2.0)
            q ~ p.q - n.q
            T ~ p.T - n.T
           
           ]
    @named sys_1 = compose(ODESystem(eqs, t, sts, []; name=name), [p, n])
    
    eqs = [
        T ~ 0,
        q ~ 0
    ]
    extend(ODESystem(eqs, t; name = name), sys_1)
end

# Add voltage for default value 
@component function Capacitor(; name, C = 1.0, temp = 1.0)
    @named oneport = OnePort()
    @unpack T, q = oneport
    ps = @parameters C = C temp = temp
    # sts = @variables T(t) = V
    D = Differential(t)
    eqs = [
        D(T) ~ q / C,
    ]
    extend(ODESystem(eqs, t, [], ps; name = name, defaults = [T => temp]), oneport)
end

@component function ConstantVoltage(; name, temp = 1.0)
    @named oneport = OnePort()
    @unpack T = oneport
    ps = @parameters temp = temp
    eqs = [
        T ~ temp,
    ]
    extend(ODESystem(eqs, t, [], ps; name = name), oneport)
end

@component function VariableVoltage(Func; name, amp = 1.0, cst = 1.0)
    @named oneport = OnePort()
    @unpack T = oneport
    ps = @parameters amp=amp cst=cst
    sts = @variables F(t)
    eqs = [
       F ~ Func,
       T ~ F*amp + cst
    ]
    extend(ODESystem(eqs, t, sts, ps; name = name), oneport)
end

@component function VariableCurrent_evolve_var_pins_1(heat, nb_pin, T_soil,sign,  timer_,Velocity; name)
          

    
    @parameters x
    @parameters tspan
    @parameters DeltaT
    @variables t
    @variables T_1(..) cumuSum(..) T_2(..) 

    Velocity = Velocity
    R1 = Rcond
    R2 = Rgg
    R3 = Rsoil_1
    C1 = C_pipe
    C2 = C_grout
    
    heat=heat
    sign = sign 
    T_soil = T_soil
    tspan = timer_



    xmin = 0.0
    xmax = 1.0
    Dt = Differential(t)
    Dx = Differential(x)
    Dxx = Differential(x)^2
    nb_pin = nb_pin


    @named ground = Ground()
    @named soil_temp_1 = ConstantVoltage(temp=T_soil)

    component_ODE = Vector{ODESystem}()
    
    for index = 1:(nb_pin) 

        ODE_system_pin_i= Pin(name = Symbol(string("ODE_system_pin_1_",index)))
        push!(component_ODE, ODE_system_pin_i)
    
    end 
    
    component_ODE_2 = Vector{ODESystem}()
    
    for index = 1:(nb_pin) 

        ODE_system_pin_i= Pin(name = Symbol(string("ODE_system_pin_2_",index)))
        push!(component_ODE_2, ODE_system_pin_i)
    
    end 
    Heat_flow(x) = sum((component_ODE[k].q*f_c_d(x*nb_pin,k) for k=1:(nb_pin)))
    
    Heat_flow_2(x) = sum((component_ODE_2[k].q*f_c_d(x*nb_pin,k) for k=1:(nb_pin)))
    
    ## Q_pin vector with integral of pins voltage 
    
    
    
    Velocity = Velocity*sign
    Source_term = (1/((c_fluid*ro_fluid)))
    Diffusivity = k_f*Source_term
    #Diffusivity = 0 
    # x_min = 0.0:0.1:0.9
    # x_max = 0.1:0.1:1.0

    m = mass_density_fluid*volume_flow_rate
    Delta = heat/(m*c_fluid*ro_fluid)

    # Ix(xmin) = Integral(x in DomainSets.ClosedInterval(xmin, x))
    # @register_symbolic Ix(xmin)
    
    # f_IX_(x) = sum((Ix(x_min[k])(Heat_flow(x)))*f_c_d(x*nb_pin, k) for k=1:nb_pin)

    # @register_symbolic f_IX_(x)

    eq = [
        Dt(T_1(t, x)) - Diffusivity*Dxx(T_1(t, x)) + Velocity*Dx(T_1(t, x)) + Source_term*(Heat_flow(x)) ~ 0,
        Dt(T_2(t, x)) - Diffusivity*Dxx(T_2(t, x)) - Velocity*Dx(T_2(t, x)) + Source_term*(Heat_flow_2(x)) ~ 0,
       
    ]
    #We can't access to the variables used for discretization but we can use T_wall 
    # because cumusum and u aren't declare as the same as after discretization
    bcs = [T_1(t, xmax) ~ T_2(t, xmax), T_1(t, xmin) ~ 1 + T_2(t, xmin), T_1(0.0, x) ~ T_soil, T_2(0.0, x) ~ T_soil]

    domains = [t ∈ Interval(0.0, tspan), x ∈ Interval(xmin, xmax)]

    @named pde_system = PDESystem(eq, bcs, domains, [t, x], [T_1(t, x), T_2(t, x)], [DeltaT => Delta])
    discretization = MOLFiniteDifference([x => nb_pin], t)

    
    molsys, tspan = symbolic_discretize(pde_system,discretization)
    
    full_equations(molsys)

    @unpack T_1 = molsys
    @unpack T_2 = molsys

    equations_sys = Vector{Equation}()
    for k =1:(nb_pin)
        
        push!(equations_sys, component_ODE[k].q ~ c_fluid*ro_fluid*(((Symbolics.scalarize(T_1)[k]))-component_ODE[k].T))
        push!(equations_sys, component_ODE_2[k].q ~ c_fluid*ro_fluid*(((Symbolics.scalarize(T_2)[k]))-component_ODE_2[k].T))
        
    end 

    
    # sys_connection = ODESystem(eq, t; name =:Discretisation)
    
    sys_eq_1 = ODESystem(equations_sys, t; name=:Discretisation_1)
    sys_comp = compose(sys_eq_1,[component_ODE..., component_ODE_2...])
    
    sys = extend(sys_comp, molsys)
end
@component function VariableCurrent_evolve_var_pins_2(nb_pin, T_soil, sign, timer_, Velocity;  name)

    @parameters x
    @parameters tspan
    @variables t
    @variables u(..) cumuSum(..)

    Velocity = Velocity
    R1 = Rcond
    R2 = Rgg
    R3 = Rsoil_1
    C1 = C_pipe
    C2 = C_grout
    sign = sign 
    
    T_soil = T_soil
    tspan = timer_



    xmin = 0.0
    xmax = 1.0
    Dt = Differential(t)
    Dx = Differential(x)
    Dxx = Differential(x)^2
    nb_pin = nb_pin


    @named ground = Ground()
    @named soil_temp_1 = ConstantVoltage(V=T_soil)

    component_ODE = Vector{ODESystem}()
    
    for index = 1:(nb_pin) 

        ODE_system_pin_i= Pin(name = Symbol(string("ODE_system_pin_2_",index)))
        push!(component_ODE, ODE_system_pin_i)
    
    end 
    
    Heat_flow(x) = sum((component_ODE[k].i*f_c_d(x*nb_pin,k) for k=1:(nb_pin)))
    
    ## Q_pin vector with integral of pins voltage 
    
    
    
    Velocity = Velocity*sign
    Source_term = (1/(c_fluid*ro_fluid))
    # Diffusivity = k_f*Source_term
    Diffusivity = 0 

    # x_min = 0.0:0.1:0.9
    # x_max = 0.1:0.1:1.0


    # Ix(xmin) = Integral(x in DomainSets.ClosedInterval(xmin, x))
    # @register_symbolic Ix(xmin)
    
    # f_IX_(x) = sum((Ix(x_min[k])(Heat_flow(x)))*f_c_d(x*nb_pin, k) for k=1:nb_pin)

    # @register_symbolic f_IX_(x)

    eq = [
        Dt(T_1(t, x)) - Diffusivity*Dxx(T_1(t, x)) + Velocity*Dx(T_1(t, x)) - Source_term*(Heat_flow(x)) ~ 0,
    ]
    #We can't access to the variables used for discretization but we can use T_wall 
    # because cumusum and u aren't declare as the same as after discretization
    bcs = [u(0.0, x) ~ T_soil]

    domains = [t ∈ Interval(0.0, tspan), x ∈ Interval(xmin, xmax)]

    @named pde_system = PDESystem(eq, bcs, domains, [t, x], [T_1(t, x)])
    discretization = MOLFiniteDifference([x => nb_pin], t)
    # ds = DiscreteSpace(domains, [u(t,x).Tal], [x.Tal], discretization)
    
    molsys, tspan = symbolic_discretize(pde_system,discretization)
    

    @unpack u = molsys


    equations_sys = Vector{Equation}()
    for k =1:(nb_pin)

        push!(equations_sys, component_ODE[k].q ~ k_f*((Symbolics.scalarize(u)[k])-component_ODE[k].T))
        
    end 

    
    # sys_connection = ODESystem(eq, t; name =:Discretisation)
    
    sys_eq_1 = ODESystem(equations_sys, t; name=:Discretisation_2)
    sys_comp = compose(sys_eq_1,[component_ODE...])
    
    sys = extend(sys_comp, molsys)

end

@component function Fluid_2(nb_pin, sign, timer_;  name, u_1=1.0)

    @parameters x
    @parameters tspan
    @variables t
    @variables T_2(..) cumuSum(..)
    
    ps= @parameters u_1 = u_1
    
    R1 = Rcond
    R2 = Rgg
    R3 = Rsoil_1
    C1 = C_pipe
    C2 = C_grout
    sign = sign 
    
    tspan = timer_



    xmin = 0.0
    xmax = 1.0
    Dt = Differential(t)
    Dx = Differential(x)
    Dxx = Differential(x)^2
    nb_pin = nb_pin


    @named ground = Ground()
    @named soil_temp_1 = ConstantVoltage(V=T_soil)

    component_ODE = Vector{ODESystem}()
    
    for index = 1:(nb_pin) 

        ODE_system_pin_i= Pin(name = Symbol(string("ODE_system_pin_2_",index)))
        push!(component_ODE, ODE_system_pin_i)
    
    end 
    
    Heat_flow(x) = sum((component_ODE[k].i*f_c_d(x*nb_pin,k) for k=1:(nb_pin)))
    
    ## Q_pin vector with integral of pins voltage 
    
    
    
    Velocity = 0.0001*sign
    Source_term = (1/((4159*1000)))
    Diffusivity = 100*Source_term

    # x_min = 0.0:0.1:0.9
    # x_max = 0.1:0.1:1.0


    # Ix(xmin) = Integral(x in DomainSets.ClosedInterval(xmin, x))
    # @register_symbolic Ix(xmin)
    
    # f_IX_(x) = sum((Ix(x_min[k])(Heat_flow(x)))*f_c_d(x*nb_pin, k) for k=1:nb_pin)

    # @register_symbolic f_IX_(x)

    eq = [
        Dt(T_2(t, x)) - Diffusivity*Dxx(T_2(t, x)) + Velocity*Dx(T_2(t, x)) - Source_term*(Heat_flow(x)) ~ 0,
    ]
    #We can't access to the variables used for discretization but we can use T_wall 
    # because cumusum and u aren't declare as the same as after discretization
    bcs = [T_2(0.0, x) ~ u_1+5, T_2(t, 0.0) ~ u_1+5, Dx(T_2(t, xmin)) ~ 0.0]

    domains = [t ∈ Interval(0.0, tspan), x ∈ Interval(xmin, xmax)]

    @named pde_system = PDESystem(eq, bcs, domains, [t, x], [T_2(t, x)])
    discretization = MOLFiniteDifference([x => nb_pin], t)
    # ds = DiscreteSpace(domains, [u(t,x).Tal], [x.Tal], discretization)
    
    molsys, tspan = symbolic_discretize(pde_system,discretization)
    

    @unpack T_2 = molsys


    equations_sys = Vector{Equation}()
    for k =1:(nb_pin)

        push!(equations_sys, component_ODE[k].q ~ Velocity*(1/Source_term)*(((Symbolics.scalarize(T_2)[k]))-component_ODE[k].T))
        
    end 

    # push!()
    # sys_connection = ODESystem(eq, t; name =:Discretisation)
    
    sys_eq_1 = ODESystem(equations_sys, t, [], ps; name=:Discretisation_2)
    sys_comp = compose(sys_eq_1,[component_ODE...])
    
    sys = extend(sys_comp, molsys)

end

@component function VariableCurrent_evolve(sign;  name, tspan=1.0)
        
    ps = @parameters t, x tspan=tspan
    @variables u(..) cumuSum(..) T_wall(t) = 1.0
    #  q_pin(..)[1:dx]

        
        
        # Integral limits are defined with DomainSets.ClosedInterval
    
        # current = Symbolics.scalarize(current)
        # typeof(q_pin)
        xmin = 0.0
        xmax = 1.0
        Ix = Integral(x in DomainSets.ClosedInterval(xmin, xmax)) # basically cumulative sum from 0 to x


        Dt = Differential(t)
        Dx = Differential(x)
        Dxx = Differential(x)^2

        @named ODE_system_pin = Pin()
        # @named p = Pin()
        # @named n = Pin()
        sts = @variables T(t)=1.0 q(t)=1.0  
        
        
        Diffusivity = 100
        Velocity = sign*0.0001
        Source_term = (1/((4159*1000)))
        
        eq = [
            cumuSum(t, x) ~ Ix(T_1(t, x)) - T_wall, 

            Dt(T_1(t, x)) - Diffusivity*Dxx(T_1(t, x)) + Velocity*Dx(T_1(t, x)) - Source_term*(cumuSum(t,x)) ~ 0,

        ]
        #We can't access to the variables used for discretization but we can use T_wall 
        # because cumusum and u aren't declare as the same as after discretization
        bcs = [u(0.0, x) ~ 0.0, Dx(u(t, xmin)) ~ 0.0, Dx(u(t, xmax)) ~ 0]

        domains = [t ∈ Interval(0.0, tspan), x ∈ Interval(xmin, xmax)]

        @named pde_system = PDESystem(eq, bcs, domains, [t, x], [T_1(t, x), cumuSum(t,x)])
        discretization = MOLFiniteDifference([x => 10], t)
        molsys, tspan = symbolic_discretize(pde_system,discretization)
        @unpack cumuSum = molsys
        
        connected = extend(molsys,ODESystem([q ~ Symbolics.scalarize(cumuSum)[1], T_wall ~ T, T ~ ODE_system_pin.T, q ~ ODE_system_pin.q], t, sts, []; name = :connected))

        compose(connected,ODE_system_pin; name = name)
    end

@component function ConstantCurrent(; name, Heat = 1.0)
        @named oneport = OnePort()
        @unpack q = oneport
        ps = @parameters Heat = Heat
        eqs = [
            q ~ Heat,
        ]
        extend(ODESystem(eqs, t, [], ps; name = name), oneport)
end



@component function Wall_variable(; name,R1 =1.0, R2=1.0 , R3=1.0 ,C1=1.0 , C2=1.0)
    # @named oneport = OnePort()
    # @unpack T,q = oneport
    ps = @parameters R1=R1 R2=R2 R3=R3 C1=C1 C2=C2 

    @named resistor_1 = Resistor(R =R1)
    
    @named resistor_2 = Resistor(R =R2)
    
    @named resistor_3 = Resistor(R = R3)

    @named capacitor = Capacitor(C =C1)

    @named capacitor_2 = Capacitor(C=C2)

    @named ground_Wall = Ground()

    rc_eqs_2 = [
        # #oneport pin
        # connect(oneport.n, resistor_3.n),
        # connect(oneport.p, resistor_1.p),

        #internal connections
        connect(resistor_1.n, capacitor.p, resistor_2.p),
        connect(resistor_2.n, capacitor_2.p, resistor_3.p),
        connect(capacitor.n, capacitor_2.n, ground_Wall.g)
    ]
 
    compose(ODESystem(rc_eqs_2, t, [], ps; name = name), ground_Wall, capacitor,capacitor_2, resistor_1, resistor_2, resistor_3)

end

@component function Imp_n_var(n, T_soil; name, Rp = 1.0, Cp=1.0)
    # @named oneport = OnePort()
    # @unpack T, q = oneport
    ps = @parameters R=Rp C=Cp 
    T_soil = T_soil
    number = n 
    index_comp = 0
    index = 0
    @named ground= Ground()
    
    @named resistor_0 = Resistor(R = Rp)
    
    @named capacitor_0 = Capacitor(C = Cp, temp = T_soil)
    
    rc_eqs_2 = [
        
    # connect(oneport.p, resistor_0.p)
    connect(resistor_0.n, capacitor_0.p)
    connect(capacitor_0.n, ground.g)
    
    ]
    component = [ground, resistor_0, capacitor_0]

    if number >= 1

        index = 1
        index_comp= 3
        for index = 1:number
            ## Create the R 
            
            resistor_i= Resistor(name = Symbol(string("R",index)); R=Rp)
            
            # Create the C 
            capacitor_i= Capacitor(name = Symbol(string("C", index)); C=Cp, temp = T_soil)
            
            # add them to the component Array 
            push!(component, resistor_i)
            push!(component, capacitor_i)
            index_comp += 2
            # connect the the rest depending of the index and the n size 
            
            if (index_comp-2) != 0
                push!(rc_eqs_2, connect(resistor_i.p, component[index_comp-2].p, component[index_comp-3].n))
                
            end

            push!(rc_eqs_2, connect(capacitor_i.p, resistor_i.n))
            push!(rc_eqs_2, connect(capacitor_i.n, ground.g))
    
            
            # if index == number 
            # end 
            # connect them 
                
            # push!(ps,ps_i)
        end        
       
        
        
    else
        #error("n can't be less than 2 or higher than 10 ");
    end

    # order of the ODEsystem systems matter

    prob =  compose(ODESystem(rc_eqs_2, t, [], ps; name = name), component)
    return prob
end 




@component function Imp_n_var_RC_var(n, T_soil; name, Rp = [], Cp=[])
    # @named oneport = OnePort()
    # @unpack T, q = oneport
    ps = @parameters R=Rp C=Cp 
    T_soil = T_soil
    number = n 
    index_comp = 0
    index = 0
    @named ground= Ground()
    
    @named resistor_0 = Resistor(R = Rp[1])
    
    @named capacitor_0 = Capacitor(C = Cp[1], temp = T_soil)
    
    rc_eqs_2 = [
        
    # connect(oneport.p, resistor_0.p)
    connect(resistor_0.n, capacitor_0.p)
    connect(capacitor_0.n, ground.g)
    
    ]
    component = [ground, resistor_0, capacitor_0]

    if number >= 1

        index = 1
        index_comp =3
        for index = 1:number
            ## Create the R 
            
            resistor_i= Resistor(name = Symbol(string("R",index)); R=Rp[index+1])
            
            # Create the C 
            capacitor_i= Capacitor(name = Symbol(string("C", index)); C=Cp[index+1], temp= T_soil)
            
            # add them to the component Array 
            push!(component, resistor_i)
            push!(component, capacitor_i)
            index_comp += 2
            # connect the the rest depending of the index and the n size 
            
            if (index_comp-2) != 0
                push!(rc_eqs_2, connect(resistor_i.p, component[index_comp-2].p, component[index_comp-3].n))
                
            end

            push!(rc_eqs_2, connect(capacitor_i.p, resistor_i.n))
            push!(rc_eqs_2, connect(capacitor_i.n, ground.g))
    
            
            # if index == number 
            # end 
            # connect them 
                
            # push!(ps,ps_i)
        end        
       
        
        
    else
        #error("n can't be less than 2 or higher than 10 ");
    end

    # order of the ODEsystem systems matter

    prob =  compose(ODESystem(rc_eqs_2, t, [], ps; name = name), component)
    return prob
end 

@component function soil_MTRCM(T_soil;name)

    ### INSTANCIATION OF COMPONENT ______________________________________________
    @named fluid_1 = Imp_n_var(1, T_soil; Rp=Rf, Cp=C_fluid)
    @named fluid_2 = Imp_n_var(1, T_soil; Rp=Rf, Cp=C_fluid)

    @named source_g = ConstantVoltage(temp=T_soil)

    @named pipe_1 = Imp_n_var(0, T_soil; Rp=Rp, Cp=C_pipe)
    @named pipe_2 = Imp_n_var(0, T_soil; Rp=Rp, Cp=C_pipe)
    @named res_pipe_1 = Resistor(R = Rp)
    @named res_pipe_2 = Resistor(R = Rp)

    @named grout_grout = Imp_n_var(1, T_soil; Rp = Rgg, Cp= C_gg)
    @named res_gg = Resistor(R = Rgg)

    @named cap_joint_1 = Capacitor(C= C_g_p_gg, temp = T_soil)
    @named cap_joint_2 = Capacitor(C= C_g_p_gg, temp = T_soil)
    @named cap_joint_3 = Capacitor(C= C_g_s, temp = T_soil)

    @named grout_1 = Imp_n_var(1, T_soil; Rp = Rg, Cp=C_grout)
    @named grout_2 = Imp_n_var(1, T_soil; Rp = Rg, Cp=C_grout)
    @named res_grout_1= Resistor(R = Rg)
    @named res_grout_2= Resistor(R = Rg)

    @named soil = Imp_n_var_RC_var(1, T_soil; Rp = [Rsoil_1 Rsoil_2], Cp = [C_soil_2 C_soil_3])
    @named res_soil= Resistor(R = Rsoil_3)

    @named ground = Ground()
    ### connection of component ______________________________________________
    rc_eqs_2 = [
    # # connection branch 1 

            connect(fluid_1.R1.n, fluid_1.C1.p, pipe_1.resistor_0.p),
            connect(pipe_1.resistor_0.n, pipe_1.capacitor_0.p, res_pipe_1.p),
            connect(res_pipe_1.n, cap_joint_1.p, grout_grout.resistor_0.p, grout_1.resistor_0.p),

    # connection branch 2 and grout-grout

            
            connect(fluid_2.R1.n, fluid_2.C1.p, pipe_2.resistor_0.p),
            connect(pipe_2.resistor_0.n, pipe_2.capacitor_0.p, res_pipe_2.p),

            connect(res_gg.p, grout_grout.R1.n, grout_grout.C1.p),
            connect(res_pipe_2.n, cap_joint_2.p, res_gg.n, grout_2.resistor_0.p),

    #connection branch 1 grout         
            connect(grout_1.R1.n, grout_1.C1.p, res_grout_1.p),
            
    #connection branch 2 grout         
            connect(grout_2.R1.n, grout_2.C1.p, res_grout_2.p), 
            connect(res_grout_1.n, cap_joint_3.p, res_grout_2.n, soil.resistor_0.p),
            
    #connection soil     
            connect(soil.R1.n, soil.C1.p, res_soil.p),
            connect(res_soil.n, source_g.p),

    #Connection to the ground 
            connect(source_g.n, ground.g, fluid_1.ground.g,fluid_2.ground.g,pipe_1.ground.g,pipe_2.ground.g, grout_grout.ground.g, grout_1.ground.g, grout_2.ground.g, soil.ground.g, cap_joint_1.n, cap_joint_2.n, cap_joint_3.n),
        
    ]

    # @named model = OnePort()
    @named _rc_model_2 = ODESystem(rc_eqs_2, t)
    rc_model_2 =compose(_rc_model_2, [source_g,fluid_1, fluid_2, pipe_1, pipe_2, grout_grout, grout_1, grout_2, soil, res_gg,res_grout_1,res_grout_2,res_pipe_1,res_pipe_2,res_soil,cap_joint_1,cap_joint_2,cap_joint_3, ground]; name=name)

end

@component function soil_MTRCM_var(n,T_soil;name)
    


    @named Source_Fluid_1 = VariableCurrent_evolve(1)

    @named Source_Fluid_2 = VariableCurrent_evolve(-1)

    # @named pump = Distributor(Const_, L)

    @named layer_0 = soil_MTRCM(T_soil)
    
    number = n
    component = [layer_0, Source_Fluid_1,Source_Fluid_2]
    rc_eqs_2 = [
        
    # connect(Source_Fluid_1.n, pump.p)
    # connect(Source_Fluid_2.n, pump.n)

    connect(layer_0.fluid_2.resistor_0.p, Source_Fluid_2.ODE_system_pin)
    connect(layer_0.fluid_1.resistor_0.p, Source_Fluid_1.ODE_system_pin)
    ]
    #Instantiation of the other soil Impedances 

    if number >= 1

        index = 1
        for index = 1:number
            ## Create the R 

            layer_q = soil_MTRCM(T_soil;name = Symbol(string("layer_",index)))

            push!(component, layer_i)

            
            push!(rc_eqs_2, connect(layer_i.fluid_2.resistor_0.p, Source_Fluid_2.ODE_system_pin))
            push!(rc_eqs_2, connect(layer_i.fluid_1.resistor_0.p, Source_Fluid_1.ODE_system_pin))
           
        
        end        

    else
        #error("n can't be less than 2 or higher than 10 ");
    end
    compose(ODESystem(rc_eqs_2, t; name = name), component)
            
end

@component function soil_MTRCM_var_pin_var(n,T_soil, tspan;name)
    


    @named Source_Fluid_1 = VariableCurrent_evolve_var_pins_1(10,1;tspan=tspan)

    @named Source_Fluid_2 = VariableCurrent_evolve_var_pins_2(10,-1;tspan=tspan)

    # @named pump = Distributor(Const_, L)

    @named layer_1 = soil_MTRCM(T_soil)

    number = n
    component = [layer_1, ModelingToolkit.get_systems(Source_Fluid_1)... ,ModelingToolkit.get_systems(Source_Fluid_2)... ]
    rc_eqs_2 = [
        
    # connect(Source_Fluid_1.n, pump.p)
    # connect(Source_Fluid_2.n, pump.n)

    connect(layer_1.fluid_2.resistor_0.p, ModelingToolkit.get_systems(Source_Fluid_2)[1])
    connect(layer_1.fluid_1.resistor_0.p, ModelingToolkit.get_systems(Source_Fluid_1)[1])
    ]
    #Instantiation of the other soil Impedances 

    if number > 1

        index = 2
        for index = 2:number
            ## Create the R 

            layer_q = soil_MTRCM(T_soil;name = Symbol(string("layer_",index)))

            push!(component, layer_i)

            
            push!(rc_eqs_2, connect(layer_i.fluid_2.resistor_0.p, ModelingToolkit.get_systems(Source_Fluid_2)[index]))
            push!(rc_eqs_2, connect(layer_i.fluid_1.resistor_0.p, ModelingToolkit.get_systems(Source_Fluid_1)[index]))
            
        end        

    else
        #error("n can't be less than 2 or higher than 10 ");
    end
    compose(ODESystem(rc_eqs_2, t; name = name), component)
            
end

@component function Resistor_conv(V; name)

    @named oneport = OnePort()
    @unpack T,q = oneport
    sts = @variables R(t)
    B = 0.0033 
    visc = 16.2*10^6 
    Pr =  0.71

    eqs = [
        #R ~ R_therm_conv(h_convection(Nusselt(Grashof(10, B, abs(v), 0.5,visc),Pr),k_soil,0.5), V),
        R ~ 0.3, 
        T ~ q * R
        
    ]
    extend(ODESystem(eqs, t, sts, []; name = name), oneport)
end

@component function soil_MTRCM_var_pin_var_ver_2(heat, n, T_soil, tspan, Velocity; name)
    heat = heat
    Velocity = Velocity
    n = n 
    T_soil = T_soil
    tspan = tspan 
    @named Source_Fluid_1 = VariableCurrent_evolve_var_pins_1(heat, n, T_soil, 1, tspan, Velocity)
 
    # LAYER ----------------------------------------------------------------------------

    
    @named layer_1 = soil_MTRCM(T_soil)


    layers = Vector{ODESystem}()
    Rth_layers = Vector{ODESystem}()
    rc_eqs_2 = Vector{Equation}()

    layers = [
        layer_1
    ]
    rc_eqs_2 = [

    
    #connect(layer_1.fluid_2.resistor_0.p, getproperty(Source_Fluid_1, Symbol(string("ODE_system_pin_2_",n))))
    
    connect(layer_1.fluid_2.resistor_0.p, Source_Fluid_1.ODE_system_pin_2_1)
    connect(layer_1.fluid_1.resistor_0.p, Source_Fluid_1.ODE_system_pin_1_1)

    ]

    #  n= 10
    #  number = 9
    number = n
    ## Create the R 
    
    T_air = ConstantVoltage(;temp=20, name = Symbol(string("T_air_ambient")))
    push!(Rth_layers, T_air)
    #computing Res conT and R cond for air ambient 
    # considering 0,5 meter deepth between air and the first layer of soil (MRCTM)
    #tqking a 1 meter area 

    #Create resistor conT 1

    
    Rth_conv_1 = Resistor_conv(Vsoil_1;name = Symbol(string("Rth_conv_1")))

    Rth_cond_1 = Resistor(;R = R_therm_soil(0.5, k_soil), name = Symbol(string("Rth_cond_1")))
    
    push!(Rth_layers, Rth_conv_1)
    push!(Rth_layers, Rth_cond_1)

 
   #  push!(rc_eqs_2, Delta_T ~ layer_1.cap_joint_3.p.T - T_air.T)

    push!(rc_eqs_2, connect(Rth_conv_1.p, T_air.p))
    push!(rc_eqs_2, connect(T_air.n, layer_1.ground.g))
    
    push!(rc_eqs_2, connect(Rth_cond_1.p, Rth_conv_1.n))
    push!(rc_eqs_2, connect(Rth_cond_1.n, layer_1.cap_joint_3.p))

    #Create resistor conT 2

    
    Rth_conv_2 = Resistor_conv(Vsoil_2; name = Symbol(string("Rth_conv_2")))

    Rth_cond_2= Resistor(;R = R_therm_soil(0.5, k_soil), name = Symbol(string("Rth_cond_2")))
    
    push!(Rth_layers, Rth_conv_2)
    push!(Rth_layers, Rth_cond_2)

 
   #  push!(rc_eqs_2, Delta_T ~ layer_1.cap_joint_3.p.T - T_air.T)

    push!(rc_eqs_2, connect(Rth_conv_2.p, T_air.p))
    push!(rc_eqs_2, connect(T_air.n, layer_1.ground.g))
    
    push!(rc_eqs_2, connect(Rth_cond_2.p, Rth_conv_2.n))
    push!(rc_eqs_2, connect(Rth_cond_2.n, layer_1.soil.capacitor_0.p))

     
    if number > 1
        # 2 to 10
        index = 2

            
        for index = 2:number
         


            
            layer_i = soil_MTRCM(T_soil;name = Symbol(string("layer_",index)))

            push!(layers, layer_i)
            # connect pin 9 (10+1 -2 = 9)..2 to layer_.._fluid_2
            #push!(rc_eqs_2, connect(layer_i.fluid_2.resistor_0.p, getproperty(Source_Fluid_1, Symbol(string("ODE_system_pin_2_",((n+1)-index))))))
            push!(rc_eqs_2, connect(layer_i.fluid_2.resistor_0.p, getproperty(Source_Fluid_1, Symbol(string("ODE_system_pin_2_",index)))))
            # connect pin 2..10 to layer_.._fluid_1
            push!(rc_eqs_2, connect(layer_i.fluid_1.resistor_0.p, getproperty(Source_Fluid_1, Symbol(string("ODE_system_pin_1_",index)))))
            
            
            
            # connect between layer, circle one
            Rth_i_1 = Resistor(;R = R_therm_soil(L, k_soil), name = Symbol(string("Rth_",index)))

            push!(Rth_layers, Rth_i_1)
            push!(rc_eqs_2, connect(layer_i.cap_joint_3.p, Rth_i_1.n))
            
            push!(rc_eqs_2, connect(layers[index-1].cap_joint_3.p, Rth_i_1.p))

            # connect between layer, circle one
            Rth_i_2 = Resistor(;R = R_therm_soil(L, k_soil), name = Symbol(string("Rth_2_",index)))

            push!(Rth_layers, Rth_i_2)
            push!(rc_eqs_2, connect(layer_i.soil.capacitor_0.p, Rth_i_2.n))
            
            push!(rc_eqs_2, connect(layers[index-1].soil.capacitor_0.p, Rth_i_2.p))
            
            # connect between layer, circle one
            Rth_i_3 = Resistor(;R = R_therm_soil(L, k_soil), name = Symbol(string("Rth_3",index)))

            push!(Rth_layers, Rth_i_3)
            push!(rc_eqs_2, connect(layer_i.soil.C1.p, Rth_i_3.n))
            
            push!(rc_eqs_2, connect(layers[index-1].soil.C1.p, Rth_i_3.p))
            
        end        

    else
        #error("n can't be less than 2 or higher than 10 ");
    end

    # # SOLVING ----------------------------------------------------------------------------

    @named sys_2 = ODESystem(rc_eqs_2, t)
    
    sys_3= compose(sys_2, [Source_Fluid_1, layers..., Rth_layers...]; name =name)
    

end

@component function soil_n_var_simple(n, timer_; name, Rp = 1.0, Cp = 1.0, T_soil=1.0)

    ps = @parameters R=Rp C=Cp T_soil=T_soil

    # @parameters x
    
    number = n
    index = 0 
    
    @named soil_temp_1 = ConstantVoltage(temp=T_soil)
    
    @named Source_Fluid = VariableCurrent_evolve_var_pins_2(number,1, timer_)
    
    # @named test = VariableCurrent_evolve(1)
    pins = [Source_Fluid.ODE_system_pin_2_1, Source_Fluid.ODE_system_pin_2_2, Source_Fluid.ODE_system_pin_2_3, Source_Fluid.ODE_system_pin_2_4, Source_Fluid.ODE_system_pin_2_5, Source_Fluid.ODE_system_pin_2_6]
    # pins = ModelingToolkit.get_systems(Source_Fluid)
    # full_equations(Source_Fluid)
    @named ground= Ground()
    
    @named layer_1 = Resistor(R = Rp)

    
    layers = Vector{ODESystem}()
    soil_source = Vector{ODESystem}()
    
    connections = Vector{Equation}()
    
    rc_eqs_2 = Vector{Equation}()
    soil_source= [
        soil_temp_1

    ]
    layers = [
        layer_1
    ]
    rc_eqs_2 = [
        
    connect(layer_1.n, soil_temp_1.p)
    connect(soil_temp_1.n, ground.g)

    ]
    
    
    #Instantiation of the other soil Impedances 
    
    
    connections = [

        connect(layer_1.p, pins[1])
    ]

    if number > 1
    
        for index = 2:number
            ## Create the R 
    
    
            
            layer_q = Resistor(name = Symbol(string("layer_",index));R = Rp)
            soil_temp_q = ConstantVoltage(name = Symbol(string("soil_temp_",index)); V=T_soil)
            
            push!(layers, layer_i)
            push!(soil_source, soil_temp_i)

            
        
            push!(rc_eqs_2, connect(layer_i.n,soil_temp_i.p))
            push!(rc_eqs_2, connect(ground.g,soil_temp_i.n))

    
            push!(connections, connect(layer_i.p, pins[index]))
        
        end        
    
    else
        #error("n can't be less than 2 or higher than 10 ");
    end
    

    
    
    compose(ODESystem([connections..., rc_eqs_2...], t; name = :rc_model), [layers...,soil_source..., pins..., Source_Fluid, ground]; name=name)
    
end

@component function soil_n_var_layer(n,timer_, T_soil_temp;name)

    @parameters x
    @parameters tspan
    @variables t
    @variables u(..) cumuSum(..)
    
    T_soil = 20
    
    R1 = Rcond
    R2 = Rgg
    R3 = Rsoil_1
    C1 = C_pipe
    C2 = C_grout

    tspan = timer_



    xmin = 0.0
    xmax = 1.0
    Dt = Differential(t)
    Dx = Differential(x)
    Dxx = Differential(x)^2
    nb_pin = n


    @named ground = Ground()
    @named soil_temp_1 = ConstantVoltage(V=T_soil_temp)

    component_ODE = Vector{ODESystem}()
    
    for index = 1:(nb_pin) 

        ODE_system_pin_i= Pin(name = Symbol(string("ODE_system_pin_2_",index)))
        push!(component_ODE, ODE_system_pin_i)
    
    end 
    
    Heat_flow(x) = sum((component_ODE[k].i*f_c_d(x*nb_pin,k) for k=1:(nb_pin)))
    
    ## Q_pin vector with integral of pins voltage 
    
    
    
    Velocity = 0.0001
    Source_term = (1/((4159*1000)))
    Diffusivity = 100*Source_term

    # x_min = 0.0:0.1:0.9
    # x_max = 0.1:0.1:1.0


    # Ix(xmin) = Integral(x in DomainSets.ClosedInterval(xmin, x))
    # @register_symbolic Ix(xmin)
    
    # f_IX_(x) = sum((Ix(x_min[k])(Heat_flow(x)))*f_c_d(x*nb_pin, k) for k=1:nb_pin)

    # @register_symbolic f_IX_(x)

    eq = [
        Dt(T_1(t, x)) - Diffusivity*Dxx(T_1(t, x)) + Velocity*Dx(T_1(t, x)) - Source_term*(Heat_flow(x)) ~ 0,
    ]
    #We can't access to the variables used for discretization but we can use T_wall 
    # because cumusum and u aren't declare as the same as after discretization
    bcs = [u(0.0, x) ~ T_soil, Dx(u(t, xmin)) ~ 0.0, Dx(u(t, xmax)) ~ 0.0]

    domains = [t ∈ Interval(0.0, tspan), x ∈ Interval(xmin, xmax)]

    @named pde_system = PDESystem(eq, bcs, domains, [t, x], [T_1(t, x)])
    discretization = MOLFiniteDifference([x => nb_pin], t)
    # ds = DiscreteSpace(domains, [u(t,x).Tal], [x.Tal], discretization)
    
    molsys, tspan = symbolic_discretize(pde_system,discretization)
    

    @unpack u = molsys


    equations_sys = Vector{Equation}()
    for k =1:(nb_pin)

        push!(equations_sys, component_ODE[k].q ~ Velocity*(1/Source_term)*(((Symbolics.scalarize(u)[k]))-component_ODE[k].T))
        
    end 

    
    # sys_connection = ODESystem(eq, t; name =:Discretisation)
    
    sys_eq_1 = ODESystem(equations_sys, t; name=:Discretisation)
    sys_comp = compose(sys_eq_1,[component_ODE...])
    
    sys = extend(sys_comp, molsys)

    # LAYER ----------------------------------------------------------------------------

    
    @named layer_1 = Wall_variable(; R1 =R1, R2= R2, R3=R3,C1=C1, C2=C2)

    layers = Vector{ODESystem}()
    rc_eqs_2 = Vector{Equation}()

    layers = [
        layer_1
    ]
    rc_eqs_2 = [
    connect(layer_1.resistor_3.n, soil_temp_1.p)
    connect(soil_temp_1.n, ground.g)
    connect(layer_1.resistor_1.p, sys.ODE_system_pin_2_1)
    ]

    

    number = nb_pin
    if number > 1

        for index = 2:number
            ## Create the R 

            
            layer_q = Wall_variable(name = Symbol(string("layer_",index));R1 =R1, R2=R2, R3=R3,C1=C1 , C2=C2)
            # soil_temp_q = ConstantVoltage(name = Symbol(string("soil_temp_",index)); V=T_soil)
            
            push!(layers, layer_i)
            # push!(soil_source, soil_temp_i)

            
        
            push!(rc_eqs_2, connect(layer_i.resistor_3.n,soil_temp_1.p))
            # push!(rc_eqs_2, connect(soil_temp_i.n,ground.g))
            pin_q = string("ODE_system_pin_2_",index)

            # fieldname_q = fieldnames(sys)

            push!(rc_eqs_2, connect(layer_i.resistor_1.p, getproperty(sys, Symbol(pin_i))))


            # push!(rc_eqs_2, connect(layer_i.resistor_1.p, component_ODE[index]))
        end        

    else
        #error("n can't be less than 2 or higher than 10 ");
    end

    # push!(rc_eqs_2, connect(layers[2].resistor_1.p, sys.ODE_system_pin_2_2))

    # push!(rc_eqs_2, connect(layers[3].resistor_1.p, sys.ODE_system_pin_2_3))

    # push!(rc_eqs_2, connect(layers[4].resistor_1.p, sys.ODE_system_pin_2_4))

    # push!(rc_eqs_2, connect(layers[5].resistor_1.p, sys.ODE_system_pin_2_5))

    # push!(rc_eqs_2, connect(layers[6].resistor_1.p, sys.ODE_system_pin_2_6))

    # SOLVING ----------------------------------------------------------------------------

    @named sys_2 = ODESystem(rc_eqs_2, t)
    
    sys_3= compose(sys_2, [sys, ground, soil_temp_1, layers...]; name =name)
    
end
