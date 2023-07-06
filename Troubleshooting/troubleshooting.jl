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


Rp = 1.0
Cp = 1.0

timer_= 86400
T_soil = 20
n = 4
ps = @parameters R=Rp C=Cp T_soil=T_soil
# @parameters x

number = n
index = 2

@named soil_temp_1 = ConstantVoltage(V=T_soil)
@named ground = Ground()
soil_source = Vector{ODESystem}()

soil_source= [
    soil_temp_1,
    ground
]

@parameters x, cumuSum(x)
@parameters tspan
tspan = timer_
@variables t
# sts = @variables q(t) = 1.0 
@variables u(..)
xmin = 0.0
xmax = 1.0
Dt = Differential(t)
Dx = Differential(x)
Dxx = Differential(x)^2
nb_pin = n

component_ODE = Vector{ODESystem}()
equations_sys = Vector{Equation}()

for index = 1:(nb_pin) 

    ## Create the R 
    
    ODE_system_pin_i= Pin(name = Symbol(string("ODE_system_pin_2_",index)))
    push!(component_ODE, ODE_system_pin_i)

end 

## Q_pin vector with integral of pins voltage 


Diffusivity = 100
Velocity = 0.0001
Source_term = (1/((4159*1000)))


# f_IX_(x_1,x_2) = Integral(x in DomainSets.ClosedInterval(x_1...,x_2...)) 

    eq = [
        Dt(u(t, x)) - Diffusivity*Dxx(u(t, x)) + Velocity*Dx(u(t, x)) - Source_term*(cumuSum...) ~ 0,
    ]
    #We can't access to the variables used for discretization but we can use T_wall 
    # because cumusum and u aren't declare as the same as after discretization
    bcs = [u(0.0, x) ~ 0.0, Dx(u(t, xmin)) ~ 0.0, Dx(u(t, xmax)) ~ 0]

    domains = [t ∈ Interval(0.0, tspan), x ∈ Interval(xmin, xmax)]

    @named pde_system = PDESystem(eq, bcs, domains, [t, x], [u(t, x)], [cumuSum => 0.0])
    discretization = MOLFiniteDifference([x => nb_pin], t)
    # ds = DiscreteSpace(domains, [u(t,x).val], [x.val], discretization)
    
    molsys, tspan = symbolic_discretize(pde_system,discretization)
    
    @unpack cumuSum = molsys
    @unpack u = molsys

    # @parameters cumuSum
    # push!(component_ODE,molsys)
    for k =1:(nb_pin)

        push!(equations_sys, component_ODE[k].i ~ Velocity*Source_term*(((Symbolics.scalarize(u)[k]))-component_ODE[k].v))
       
    end 
    # push!(equations_sys, full_equations(molsys)...)
    # push!(equations, (Symbolics.scalarize(u)))
    push!(equations_sys,  cumuSum[1] ~ sum((component_ODE[k].i*f_c_d(x,k) for k=1:(nb_pin))))
    
    sys_eq_1 = ODESystem(equations_sys, t; name=:Discretisation)
    sys_comp=compose(sys_eq_1,[component_ODE...])
    sys = extend(sys_comp, molsys)
    
    

@named layer_1 = Wall_variable(; R1 =Rp, R2= Rp, R3=Rp,C1=Cp , C2=Cp)
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



if number > 1

    for index = 2:number
        ## Create the R 


        layer_i = Wall_variable(name = Symbol(string("layer_",index));R1 =Rp, R2= Rp, R3=Rp,C1=Cp , C2=Cp)
        soil_temp_i = ConstantVoltage(name = Symbol(string("soil_temp_",index)); V=T_soil)
        
        push!(layers, layer_i)
        push!(soil_source, soil_temp_i)

        
    
        push!(rc_eqs_2, connect(layer_i.resistor_3.n,soil_temp_i.p))
        push!(rc_eqs_2, connect(soil_temp_i.n,ground.g))


        # push!(rc_eqs_2, connect(layer_i.resistor_1.p, component_ODE[index]))
    end        

else
    #error("n can't be less than 2 or higher than 10 ");
end

        push!(rc_eqs_2, connect(layers[2].resistor_1.p,  sys.ODE_system_pin_2_2))
        push!(rc_eqs_2, connect(layers[3].resistor_1.p,  sys.ODE_system_pin_2_3))
        push!(rc_eqs_2, connect(layers[4].resistor_1.p,  sys.ODE_system_pin_2_4))

  ## Après l'ajout des liaisons
#   sys = extend(ODESystem(rc_eqs_2, t; name = :connected), sys_comp_pins)
#   sys = expand_connections(sys)
#   sys= alias_elimination(sys)

    sys_eq_2 = ODESystem(rc_eqs_2, t; name =:equation_sys)


  dis_sys = compose(sys_eq_2, [sys, layers..., soil_source...])
  

  inspect = expand_connections(dis_sys)
  inspect = alias_elimination(inspect)
  full_equations(expand_connections(inspect))
#   connected = compose(dis_sys, component_ODE...)
#   sys = extend(connected,sys_extend_molsys)
# sys_extend = ODESystem([rc_eqs_2...], t; name = :connected)


# sys = compose(sys_extend, [layers..., soil_source...,discrete_pins]; name=:final)
# sys = extend(ODESystem(rc_eqs_2, t; name = :connected),sys)


sys = expand_connections(sys)
sys= alias_elimination(sys)

structural_simplify(dis_sys)