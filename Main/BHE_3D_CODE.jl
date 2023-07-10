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

T_fluid = 20
T_grout = 20

T_soil= 20
timer_= 86400
n= 10
@named complex_soil = soil_MTRCM_var_pin_var_ver_2(n, T_soil, timer_)
# 
# @named ODE_Disc_var = VariableCurrent_evolve_var_pins_2(10,20, 1, timer_)
# @named simple_layers = soil_n_var_simple(6, timer_; Rp=1.0, Cp=1.0, T_soil = 20)
# @named wall_layers = soil_n_var_layer(4, timer_, T_soil)
# @unpack cumuSum = simple_layers.connected
# @unpack pins = simple_layers.connected

# sys= wall_layers
# sys = expand_connections(sys)
# sys= alias_elimination(sys)
# cumuSum[1] ~ sum((component_ODE[k].i*f_c_d(x,k) for k=1:nb_pin))


sys_2 = structural_simplify(complex_soil)
# @named complex_soil = soil_MTRCM_var_pin_var_ver_2(10,10,timer_)

# # traced_sys = modelingtoolkitize(complex_soil)
# sys_2 = structural_simplify(complex_soil)

# @named soil_test = soil_MTRCM_var_pin_var(1, T_soil)
# @named ODE_Disc = VariableCurrent_evolve(1)
# @named Layer_1 = soil_MTRCM_var(1,T_soil)
# # # @named Layer_1 = soil_n_var_simple(0)
# # sys_2 = structural_simplify(soil_test)
# sys_2 = structural_simplify(Layer_1)

@nonamespace u0= [
    # sys_2.pump₊q => 0,
    # sys_2.soil_temp_1₊V => T_soil,
    # sys_2.soil_temp_2₊V => T_soil,
    # sys_2.soil_temp_3₊V => T_soil,
    # sys_2.soil_temp_4₊V => T_soil,
    # sys_2.soil_temp_5₊V => T_soil,
    # sys_2.soil_temp_6₊V => T_soil,
    # sys_2.connected₊x => 0.01,
    # sys_2.connected₊cumuSum=> 0.1
    # (sys_2.Source_Fluid_1₊u[k] => T_fluid for k=1:10)...,
    
    # (sys_2.Source_Fluid_2₊u[k] => T_fluid for k=1:10)...,
    # sys_2.Source_Fluid_1₊cumuSum[1] => T_fluid,   
    # sys_2.Source_Fluid_2₊cumuSum[1] => T_fluid,   
    
    # sys_2.layer_1₊capacitor₊v =>T_soil,
    
    # sys_2.layer_2₊capacitor₊v =>T_soil,
    # sys_2.layer_3₊capacitor₊v =>T_soil,
    
        # HEAT_transfer_borehole => 0.0
        sys_2.layer_1₊soil₊capacitor_0₊v =>T_soil,
        sys_2.layer_1₊soil₊C1₊v => T_soil,
        sys_2.layer_1₊cap_joint_3₊v=> T_soil,
    
        sys_2.layer_1₊grout_2₊C1₊v => T_grout,
        sys_2.layer_1₊grout_2₊capacitor_0₊v =>T_grout,
        sys_2.layer_1₊grout_1₊capacitor_0₊v =>T_grout,
        sys_2.layer_1₊grout_1₊C1₊v =>T_grout,
        
        sys_2.layer_1₊grout_grout₊capacitor_0₊v => T_grout,
        sys_2.layer_1₊grout_grout₊C1₊v => T_grout,
        sys_2.layer_1₊fluid_1₊capacitor_0₊v => T_fluid,
        sys_2.layer_1₊fluid_1₊C1₊v => T_fluid,
    
        sys_2.layer_1₊fluid_2₊capacitor_0₊v => T_fluid,
        sys_2.layer_1₊fluid_2₊C1₊v => T_fluid,
    
        sys_2.layer_1₊pipe_1₊capacitor_0₊v => T_grout,
        sys_2.layer_1₊pipe_2₊capacitor_0₊v => T_grout,
    
        sys_2.layer_1₊cap_joint_1₊v=> T_grout,
        sys_2.layer_1₊cap_joint_2₊v=> T_grout,


        sys_2.layer_2₊soil₊capacitor_0₊v =>T_soil,
        sys_2.layer_2₊soil₊C1₊v => T_soil,
        sys_2.layer_2₊cap_joint_3₊v=> T_soil,

        sys_2.layer_2₊grout_2₊C1₊v => T_grout,
        sys_2.layer_2₊grout_2₊capacitor_0₊v =>T_grout,
        sys_2.layer_2₊grout_1₊capacitor_0₊v =>T_grout,
        sys_2.layer_2₊grout_1₊C1₊v =>T_grout,
        
        sys_2.layer_2₊grout_grout₊capacitor_0₊v => T_grout,
        sys_2.layer_2₊grout_grout₊C1₊v => T_grout,
        sys_2.layer_2₊fluid_1₊capacitor_0₊v => T_fluid,
        sys_2.layer_2₊fluid_1₊C1₊v => T_fluid,

        sys_2.layer_2₊fluid_2₊capacitor_0₊v => T_fluid,
        sys_2.layer_2₊fluid_2₊C1₊v => T_fluid,
    
        sys_2.layer_2₊pipe_1₊capacitor_0₊v => T_grout,
        sys_2.layer_2₊pipe_2₊capacitor_0₊v => T_grout,
    
        sys_2.layer_2₊cap_joint_1₊v=> T_grout,
        sys_2.layer_2₊cap_joint_2₊v=> T_grout,


        
        sys_2.layer_3₊soil₊capacitor_0₊v =>T_soil,
        sys_2.layer_3₊soil₊C1₊v => T_soil,
        sys_2.layer_3₊cap_joint_3₊v=> T_soil,

        sys_2.layer_3₊grout_2₊C1₊v => T_grout,
        sys_2.layer_3₊grout_2₊capacitor_0₊v =>T_grout,
        sys_2.layer_3₊grout_1₊capacitor_0₊v =>T_grout,
        sys_2.layer_3₊grout_1₊C1₊v =>T_grout,
        
        sys_2.layer_3₊grout_grout₊capacitor_0₊v => T_grout,
        sys_2.layer_3₊grout_grout₊C1₊v => T_grout,
        sys_2.layer_3₊fluid_1₊capacitor_0₊v => T_fluid,
        sys_2.layer_3₊fluid_1₊C1₊v => T_fluid,

        sys_2.layer_3₊fluid_2₊capacitor_0₊v => T_fluid,
        sys_2.layer_3₊fluid_2₊C1₊v => T_fluid,

        sys_2.layer_3₊pipe_1₊capacitor_0₊v => T_grout,
        sys_2.layer_3₊pipe_2₊capacitor_0₊v => T_grout,
    
        sys_2.layer_3₊cap_joint_1₊v=> T_grout,
        sys_2.layer_3₊cap_joint_2₊v=> T_grout,

        


        sys_2.layer_4₊soil₊capacitor_0₊v =>T_soil,
        sys_2.layer_4₊soil₊C1₊v => T_soil,
        sys_2.layer_4₊cap_joint_3₊v=> T_soil,

        sys_2.layer_4₊grout_2₊C1₊v => T_grout,
        sys_2.layer_4₊grout_2₊capacitor_0₊v =>T_grout,
        sys_2.layer_4₊grout_1₊capacitor_0₊v =>T_grout,
        sys_2.layer_4₊grout_1₊C1₊v =>T_grout,
        
        sys_2.layer_4₊grout_grout₊capacitor_0₊v => T_grout,
        sys_2.layer_4₊grout_grout₊C1₊v => T_grout,
        sys_2.layer_4₊fluid_1₊capacitor_0₊v => T_fluid,
        sys_2.layer_4₊fluid_1₊C1₊v => T_fluid,

        sys_2.layer_4₊fluid_2₊capacitor_0₊v => T_fluid,
        sys_2.layer_4₊fluid_2₊C1₊v => T_fluid,

        sys_2.layer_4₊pipe_1₊capacitor_0₊v => T_grout,
        sys_2.layer_4₊pipe_2₊capacitor_0₊v => T_grout,
    
        sys_2.layer_4₊cap_joint_1₊v=> T_grout,
        sys_2.layer_4₊cap_joint_2₊v=> T_grout,
            # HEAT_transfer_borehole => 0.0
    ]

prob = ODEProblem(sys_2, u0, timer_)

sol = solve(prob, saveat = 3600)



discrete_x = 1:10
discrete_t = sol[t]

solu = sol[prob.f.sys.Discretisation_1₊u]
solu_2 = sol[prob.f.sys.Discretisation_2₊u]

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

