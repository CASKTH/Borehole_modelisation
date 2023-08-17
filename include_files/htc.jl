
include("../src/prop_th.jl")


"""
Function that return the heat transfer coefficient (h) in [W/m²/K]
It is to be noted that the tmperature should be the film temperature (average of bulk and wall) and not the fluid temperature.
Also, this quickly put together function works with refrigerants as a "fluid" key. The list is here:  http://www.coolprop.org/fluid_properties/PurePseudoPure.html

The correlation used has a certain validity range, but no error message is triggered if we are computing outside of it.

Moreover, two BCs are possible: constant heat flux of constant wall temperature. The constant heat flux is used here.

In the same way, two conditions are possible:
-Heating (Tf<Tw)
-Cooling (Tf>Tw)

The Dittus Boelter correlation is used for the turbulent case. (eq. 8.60 from Incropeira 2007)

(The flow is also assumed to be fully developped.)
"""
function h_circular(;ri, T, p, v, fluid, heating=true)
    ρ = ρ_(T=T, p=p, fluid=fluid) # density [kg/m³]
    μ = μ_(T=T, p=p, fluid=fluid) # viscosity [N⋅s/m²]
    Pr = Pr_(T=T, p=p, fluid=fluid) # Prandtl number []:ratio of momentum to thermal diffusivity
    λ = λ_(T=T, p=p, fluid=fluid) # thermal conductivity
    Di = ri * 2 #diameter
    Re = ρ * v * Di / μ # Reynolds number []: ratio of inertial forces to viscosity
    if Re < 5e3 #the limit is not well defined
        Nu = 4.36 #Nusselt number []: ratio of convection over conduction
        #constant HF hypothesis
    else
        heating ? n = 0.4 : n = 0.3
        Nu = 0.023*Re^(4/5)*Pr^n # see above
    end
    h = Nu * λ / Di
    return h
end

h_example = h_circular(
    ri = 0.08, # internal radius [m]
    T = 10+273.15, # temperature [K]
    p = 1e5, # pressure [Pa]
    v = 0.2, # velocity [m/s] not sure at all of this
    fluid = "water",
)
