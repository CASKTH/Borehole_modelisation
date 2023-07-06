using DataDrivenDiffEq
using Pkg
using LinearAlgebra
using OrdinaryDiffEq
using DataDrivenDMD
using Plots
using DomainSets
using MethodOfLines
using DifferentialEquations
include("../include_files/functions_component.jl")
include("../include_files/datas_code.jl")
# 


##Test 

@parameters x tspan
@variables t u(..) cumuSum(..)

# @register_symbolic Heat_flow(x)
T_soil = 20
T_fluid = 20

R1 = Rp
R2 = Rgg
R3 =Rsoil_1
C1=C_pipe
C2 = C_grout

tspan = timer_



xmin = 0.0
xmax = 1.0
Dt = Differential(t)
Dx = Differential(x)
Dxx = Differential(x)^2
nb_pin = 6

component_ODE = Vector{ODESystem}()
    
for index = 1:(nb_pin) 

    ODE_system_pin_i= Pin(name = Symbol(string("ODE_system_pin_2_",index)))
    push!(component_ODE, ODE_system_pin_i)

end 

Heat_flow(x) = sum((component_ODE[k].i*f_c_d(x*nb_pin,k) for k=1:(nb_pin)))
## Q_pin vector with integral of pins voltage 

Temperature_pins(x) = sum((component_ODE[k].v*f_c_d(x*nb_pin,k) for k=1:(nb_pin)))

Diffusivity = 100
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
    Dt(u(t, x)) - Diffusivity*Dxx(u(t, x)) + Velocity*Dx(u(t, x)) - Source_term*(Heat_flow(x)) ~ 0
]
#We can't access to the variables used for discretization but we can use T_wall 
# because cumusum and u aren't declare as the same as after discretization
bcs = [u(0.0, x) ~ T_soil, Dx(u(t, xmin)) ~ 0.0, Dx(u(t, xmax)) ~ 0.0]

domains = [t ∈ Interval(0.0, tspan), x ∈ Interval(xmin, xmax)]

@named pde_system = PDESystem(eq, bcs, domains, [t, x], [u(t, x)])
discretization = MOLFiniteDifference([x => nb_pin], t)
# ds = DiscreteSpace(domains, [u(t,x).val], [x.val], discretization)

molsys, tspan = symbolic_discretize(pde_system,discretization)


@unpack u = molsys


equations_sys = Vector{Equation}()
for k =1:(nb_pin)

    push!(equations_sys, component_ODE[k].i ~ Velocity*(1/Source_term)*(((Symbolics.scalarize(u)[k]))-component_ODE[k].v))
    
end 






# sys_connection = ODESystem(eq, t; name =:Discretisation)

sys_eq_1 = ODESystem(equations_sys, t; name=:Discretisation)
sys_comp = compose(sys_eq_1,[component_ODE...])

sys = compose(molsys, sys_comp)

# level1 = ODESystem(Equation[], t, [], []; name = :Discretisation)


# Pins1 = ModelingToolkit.get_systems(sys)[1] ∘ level1



@named layer_1 = Wall_variable(; R1 =R1, R2= R2, R3=R3,C1=C1 , C2=C2)
eq= [


    connect(layer_1.resistor_1.p, Pins1)
]

test = ODESystem(eq, t, [], []; name = :test)









@parameters t a b c d e f
p = [a #a is a local variable
    ParentScope(b) # b is a variable that belongs to one level up in the hierarchy
    ParentScope(ParentScope(c))# ParentScope can be nested
    DelayParentScope(d) # skips one level before applying ParentScope
    DelayParentScope(e, 2) # second argument allows skipping N levels
    GlobalScope(f)]

level0 = ODESystem(Equation[], t, [], p; name = :level0)
level1 = ODESystem(Equation[], t, [], []; name = :level1) ∘ level0
parameters(level1)

level2 = ODESystem(Equation[], t, [], []; name = :level2) ∘ level1
parameters(level2)
































### TEST SOIL LAYER WITH RESISTOR 
@named test_2 = soil_n_var_layer(6, 86400, 15)


@named test_4 = soil_n_var_layer(6, 86400, 25)


    sys_4 = structural_simplify(test_2)
    
    u0 = [
        # sys_4.Discretisation₊cumuSum => 0.0,
        # sys_4.soil_temp_1₊V => 0.0,


        sys_4.layer_1₊capacitor₊v => T_soil,
        
        sys_4.layer_2₊capacitor₊v => T_soil,
        sys_4.layer_3₊capacitor₊v => T_soil,
        sys_4.layer_4₊capacitor₊v => T_soil,
        
        sys_4.layer_5₊capacitor₊v => T_soil,
        sys_4.layer_6₊capacitor₊v => T_soil,

        sys_4.layer_1₊capacitor_2₊v => T_soil,
        
        sys_4.layer_2₊capacitor_2₊v => T_soil,
        sys_4.layer_3₊capacitor_2₊v => T_soil,
        sys_4.layer_4₊capacitor_2₊v => T_soil,
        
        sys_4.layer_5₊capacitor_2₊v => T_soil,
        sys_4.layer_6₊capacitor_2₊v => T_soil,
        T_soil => 20.0,
        # sys_4.soil_temp_1₊v => 0,
        Rp => 1.0,
        Cp => 1.0
        ]

    prob = ODEProblem(sys_4, u0, timer_)

    sol = solve(prob,saveat = 3600)

    
    sys_4 = structural_simplify(test_4)

    prob_2 = ODEProblem(sys_4, u0, timer_)

    sol_2 = solve(prob_2,saveat = 3600)
    
    solu = sol[prob.f.sys.Discretisation₊u]
    solu_2 = sol_2[prob_2.f.sys.Discretisation₊u]

    discrete_x = 1:6
    discrete_t = sol[t]
    
    using WGLMakie, ColorSchemes 

    x = discrete_x
    fig = Figure(resolution = (800, 600))
    ax = Axis(fig[1, 1], xlabel = "depth",xlabelsize = 1)
    obj = lines!(x, solu[1]; color = :red, label ="t=$(1)", colormap = :heat, linewidth = 1)
    obj = lines!(x, solu_2[1]; color = :blue, label ="t=$(1)", colormap = :heat, linewidth = 1)
    
    for i in 2:(length(discrete_t))
            
           lines!(x, solu[i]; color = :red, label ="t=$(i)",colormap = :heat, linewidth = 1)
           lines!(x, solu_2[i]; color = :blue, label ="t=$(i)",colormap = :heat, linewidth = 1)
    
    end
    # axislegend(ax)
    # Colorbar(fig[1, 2], obj, label = "temperature")
    
    # colgap!(fig.layout, 10)
    fig


 