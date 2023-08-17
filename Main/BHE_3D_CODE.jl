using DataDrivenDiffEq
using Pkg
using LinearAlgebra
using OrdinaryDiffEq
using DataDrivenDMD
using Plots
using DomainSets
using MethodOfLines
using ModelingToolkitDesigner
using DifferentialEquations
using WGLMakie, ColorSchemes 

# using ModelingToolkitStandardLibrary.Electrical
# using ModelingToolkitStandardLibrary.Blocks: Constant

include("../include_files/functions_component_thermal.jl")
include("../include_files/datas_code.jl")

T_fluid = 10
T_grout = 10

T_soil= 20


timer_= 86400*41
n = 10
Velocity = 0.0001
Delta = 1
@named complex_soil = soil_MTRCM_var_pin_var_ver_2(Delta, n, T_soil, timer_, Velocity)

# @named ODE_Disc_var = VariableCurrent_evolve_var_pins_1(10, 20, 1, timer_, Velocity)

# path = joinpath(@__DIR__, "design") # folder where visualization info is saved and retrieved
# design = ODESystemDesign(complex_soil, path);
# ModelingToolkitDesigner.view(design)
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
    # sys_2.soil_temp_1₊T => T_soil,
    # sys_2.soil_temp_2₊T => T_soil,
    # sys_2.soil_temp_3₊T => T_soil,
    # sys_2.soil_temp_4₊T => T_soil,
    # sys_2.soil_temp_5₊T => T_soil,
    # sys_2.soil_temp_6₊T => T_soil,
    # sys_2.connected₊x => 0.01,
    # sys_2.connected₊cumuSum=> 0.1
    # (sys_2.Source_Fluid_1₊u[k] => T_fluid for k=1:10)...,
    
    # (sys_2.Source_Fluid_2₊u[k] => T_fluid for k=1:10)...,
    # sys_2.Discretisation_1₊u[2] => T_fluid,   
    # sys_2.Discretisation_1₊u[3] => T_fluid,
    # sys_2.Source_Fluid_2₊cumuSum[1] => T_fluid,   

    # sys_2.layer_1₊capacitor₊T =>T_soil,
    
    # sys_2.layer_2₊capacitor₊T =>T_soil,
    # sys_2.layer_3₊capacitor₊T =>T_soil,
    
        # # # # HEAT_transfer_borehole => 0.0
        sys_2.Rth_conv_1₊R => 0
        sys_2.Rth_conv_2₊R => 0
        #=
     sys_2.layer_1₊soil₊capacitor_0₊T =>T_grout,
        sys_2.layer_1₊soil₊C1₊T => T_grout,
        sys_2.layer_1₊cap_joint_3₊T=> T_grout,
    
        sys_2.layer_1₊grout_2₊C1₊T => T_grout,
        sys_2.layer_1₊grout_2₊capacitor_0₊T =>T_grout,
        sys_2.layer_1₊grout_1₊capacitor_0₊T =>T_grout,
        sys_2.layer_1₊grout_1₊C1₊T =>T_grout,
        
        sys_2.layer_1₊grout_grout₊capacitor_0₊T => T_grout,
        sys_2.layer_1₊grout_grout₊C1₊T => T_grout,
        sys_2.layer_1₊fluid_1₊capacitor_0₊T => T_fluid,
        sys_2.layer_1₊fluid_1₊C1₊T => T_fluid,
    
        sys_2.layer_1₊fluid_2₊capacitor_0₊T => T_fluid,
        sys_2.layer_1₊fluid_2₊C1₊T => T_fluid,
    
        sys_2.layer_1₊pipe_1₊capacitor_0₊T => T_grout,
        sys_2.layer_1₊pipe_2₊capacitor_0₊T => T_grout,
    
        sys_2.layer_1₊cap_joint_1₊T=> T_grout,
        sys_2.layer_1₊cap_joint_2₊T=> T_grout,


        sys_2.layer_2₊soil₊capacitor_0₊T =>T_grout,
        sys_2.layer_2₊soil₊C1₊T => T_grout,
        sys_2.layer_2₊cap_joint_3₊T=> T_grout,

        sys_2.layer_2₊grout_2₊C1₊T => T_grout,
        sys_2.layer_2₊grout_2₊capacitor_0₊T =>T_grout,
        sys_2.layer_2₊grout_1₊capacitor_0₊T =>T_grout,
        sys_2.layer_2₊grout_1₊C1₊T =>T_grout,
        
        sys_2.layer_2₊grout_grout₊capacitor_0₊T => T_grout,
        sys_2.layer_2₊grout_grout₊C1₊T => T_grout,
        sys_2.layer_2₊fluid_1₊capacitor_0₊T => T_fluid,
        sys_2.layer_2₊fluid_1₊C1₊T => T_fluid,

        sys_2.layer_2₊fluid_2₊capacitor_0₊T => T_fluid,
        sys_2.layer_2₊fluid_2₊C1₊T => T_fluid,
    
        sys_2.layer_2₊pipe_1₊capacitor_0₊T => T_grout,
        sys_2.layer_2₊pipe_2₊capacitor_0₊T => T_grout,
    
        sys_2.layer_2₊cap_joint_1₊T=> T_grout,
        sys_2.layer_2₊cap_joint_2₊T=> T_grout,


        
        sys_2.layer_3₊soil₊capacitor_0₊T =>T_grout,
        sys_2.layer_3₊soil₊C1₊T => T_grout,
        sys_2.layer_3₊cap_joint_3₊T=> T_grout,

        sys_2.layer_3₊grout_2₊C1₊T => T_grout,
        sys_2.layer_3₊grout_2₊capacitor_0₊T =>T_grout,
        sys_2.layer_3₊grout_1₊capacitor_0₊T =>T_grout,
        sys_2.layer_3₊grout_1₊C1₊T =>T_grout,
        
        sys_2.layer_3₊grout_grout₊capacitor_0₊T => T_grout,
        sys_2.layer_3₊grout_grout₊C1₊T => T_grout,
        sys_2.layer_3₊fluid_1₊capacitor_0₊T => T_fluid,
        sys_2.layer_3₊fluid_1₊C1₊T => T_fluid,

        sys_2.layer_3₊fluid_2₊capacitor_0₊T => T_fluid,
        sys_2.layer_3₊fluid_2₊C1₊T => T_fluid,

        sys_2.layer_3₊pipe_1₊capacitor_0₊T => T_grout,
        sys_2.layer_3₊pipe_2₊capacitor_0₊T => T_grout,
    
        sys_2.layer_3₊cap_joint_1₊T=> T_grout,
        sys_2.layer_3₊cap_joint_2₊T=> T_grout,

        


        sys_2.layer_4₊soil₊capacitor_0₊T =>T_grout,
        sys_2.layer_4₊soil₊C1₊T => T_grout,
        sys_2.layer_4₊cap_joint_3₊T=> T_grout,

        sys_2.layer_4₊grout_2₊C1₊T => T_grout,
        sys_2.layer_4₊grout_2₊capacitor_0₊T =>T_grout,
        sys_2.layer_4₊grout_1₊capacitor_0₊T =>T_grout,
        sys_2.layer_4₊grout_1₊C1₊T =>T_grout,
        
        sys_2.layer_4₊grout_grout₊capacitor_0₊T => T_grout,
        sys_2.layer_4₊grout_grout₊C1₊T => T_grout,
        sys_2.layer_4₊fluid_1₊capacitor_0₊T => T_fluid,
        sys_2.layer_4₊fluid_1₊C1₊T => T_fluid,

        sys_2.layer_4₊fluid_2₊capacitor_0₊T => T_fluid,
        sys_2.layer_4₊fluid_2₊C1₊T => T_fluid,

        sys_2.layer_4₊pipe_1₊capacitor_0₊T => T_grout,
        sys_2.layer_4₊pipe_2₊capacitor_0₊T => T_grout,
    
        sys_2.layer_4₊cap_joint_1₊T=> T_grout,
        sys_2.layer_4₊cap_joint_2₊T=> T_grout,

        =#

        # sys_2.layer_10₊soil₊capacitor_0₊T =>T_soil,
        # sys_2.layer_10₊soil₊C1₊T => T_soil,
        # sys_2.layer_10₊cap_joint_3₊T=> T_soil,

        # sys_2.layer_10₊grout_2₊C1₊T => T_grout,
        # sys_2.layer_10₊grout_2₊capacitor_0₊T =>T_grout,
        # sys_2.layer_10₊grout_1₊capacitor_0₊T =>T_grout,
        # sys_2.layer_10₊grout_1₊C1₊T =>T_grout,
        
        # sys_2.layer_10₊grout_grout₊capacitor_0₊T => T_grout,
        # sys_2.layer_10₊grout_grout₊C1₊T => T_grout,
        # sys_2.layer_10₊fluid_1₊capacitor_0₊T => T_fluid,
        # sys_2.layer_10₊fluid_1₊C1₊T => T_fluid,

        # sys_2.layer_10₊fluid_2₊capacitor_0₊T => T_fluid,
        # sys_2.layer_10₊fluid_2₊C1₊T => T_fluid,

        # sys_2.layer_10₊pipe_1₊capacitor_0₊T => T_grout,
        # sys_2.layer_10₊pipe_2₊capacitor_0₊T => T_grout,
    
        # sys_2.layer_10₊cap_joint_1₊T=> T_grout,
        # sys_2.layer_10₊cap_joint_2₊T=> T_grout,
        #     # HEAT_transfer_borehole => 0.0
    ]

prob = ODEProblem(sys_2, u0, timer_)

# # Resolution of DTLessThanMin error on solving  => Automate index lowering
# prob_2 = modelingtoolkitize(prob)
# pendulum_sys = structural_simplify(prob_2)
# prob_3 = ODEProblem(pendulum_sys, Pair[], timer_)

sol = solve(prob, saveat = 3600)



discrete_x = 1:10

x = discrete_x
discrete_x_2 = reverse(discrete_x)
discrete_t = sol[t]

solu = sol[prob.f.sys.Discretisation_1₊T_1]
red_sol = reduce(hcat, solu)

solu_2 = sol[prob.f.sys.Discretisation_1₊T_2]
red_sol_2 = reduce(hcat, solu_2)


# fig = Figure(resolution = (800, 600))
# ax = Axis(fig[1, 1], xlabel = "depth",xlabelsize = 1)
# obj = lines!(x, solu[1]; color = :red, label ="t=$(1)", colormap = :heat, linewidth = 1)
# obj = lines!(discrete_x_2, solu_2[1]; color = :blue, label ="t=$(1)", colormap = :heat, linewidth = 1)

# obj = lines!(x, solu[2]; color = :brown, label ="t=$(2)", colormap = :heat, linewidth = 1)
# obj = lines!(discrete_x_2, solu_2[2]; color = :green, label ="t=$(2)", colormap = :heat, linewidth = 1)

# axislegend(ax)
# for i in 2:(length(discrete_t))
        
#        lines!(x, solu[i]; color = :red, label ="t=$(i)",colormap = :heat, linewidth = 1)
#        lines!(discrete_x_2, solu_2[i]; color = :blue, label ="t=$(i)",colormap = :heat, linewidth = 1)

# end

# # Colorbar(fig[1, 2], obj, label = "temperature")

# # colgap!(fig.layout, 10)
# fig

# WGLMakie.activate!()

zmin, zmax = minimum(red_sol), maximum(red_sol)
cmap = :viridis
cmap_2 = :viridis

fig = Figure(resolution = (800, 1400), fontsize = 22)
ax = Axis3(fig[1, 1], aspect = :data, perspectiveness = 0.5, elevation = π/9,
    xzpanelcolor = (:black, 0.75), yzpanelcolor = (:black, 0.75),
    zgridcolor = :grey, ygridcolor = :grey, xgridcolor = :grey)


sm_1 = WGLMakie.surface!(ax, x, discrete_t*0.000001, red_sol; colormap = cmap, colorrange = (zmin, zmax),
    transparency = true)
    #WGLMakie.wireframe!(ax, x, discrete_t*0.000001, red_sol; overdraw = true, transparency = true, color = (:black, 0.1))

ax2 = Axis3(fig[2, 1], aspect = :data, perspectiveness = 0.5, elevation = π / 9,
xzpanelcolor = (:black, 0.75), yzpanelcolor = (:black, 0.75),
zgridcolor = :grey, ygridcolor = :grey, xgridcolor = :grey)

zmin, zmax = minimum(red_sol_2), maximum(red_sol_2)
sm_2 = WGLMakie.surface!(ax2, discrete_x_2, discrete_t*0.000001, red_sol_2; colormap = cmap_2, colorrange = (zmin, zmax),
    transparency = true)

xm, ym, zm = minimum(ax.finallimits[])



Colorbar(fig[1, 2], sm_1, height = Relative(0.5))
colsize!(fig.layout, 1, Aspect(1, 1.0))

Colorbar(fig[2, 2], sm_2, height = Relative(0.5))
colsize!(fig.layout, 1, Aspect(1, 1.0))
fig

Sol_cap = Vector{Vector}()
v_n_tr = Vector{Matrix{Float64}}()

for n_z=1:n


    Sol_cap = Vector{Vector{Float64}}()
    push!(Sol_cap, sol[getproperty(prob.f.sys, Symbol(string("layer_"),n_z,string("₊cap_joint_3₊T")))])
    push!(Sol_cap, sol[getproperty(prob.f.sys, Symbol(string("layer_"),n_z,string("₊soil₊capacitor_0₊T")))])
    push!(Sol_cap, sol[getproperty(prob.f.sys, Symbol(string("layer_"),n_z,string("₊soil₊C1₊T")))])
    
    n_tr = reduce(hcat, Sol_cap)
    push!(v_n_tr, n_tr)


end


zcat(v) = cat(v...,dims=3)
M = zcat(v_n_tr[1:2])

discrete_t

Plots.contourf(v_n_tr[1],  levels=40, color=:viridis, linewidth=0)
min_z = minimum(M)
max_z = maximum(M)
anim = @animate for i ∈ 1:length(discrete_t)
    Plots.contourf(M[i,:,:]',levels=10, color=:viridis, linewidth=0, clims=(min_z, max_z), xlabel="Discrete radius", ylabel="Discrete depth");
end

gif(anim, "anim_fps15.gif", fps = 20)