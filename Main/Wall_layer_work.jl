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
# Rp


##Test 

T_soil = 20 
timer_ = 86400

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
        # Rp => Rcond,

        ]

    prob = ODEProblem(sys_4, u0, 86400)

    sol = solve(prob, saveat = 3600)

    
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


 