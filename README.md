# Borehole_modelisation
Project dedicated to borehole modilisation and command using Julia 
---

_Understand the model_ 

The model is a Acausal component based model containing : 
- The discrete convection-diffusion component giving a temperature discretisation through time and space : u(t)[x]
- The layer instantiation, function of n size of u(t)[x]
- Ploting of the solution


_How to understand the code_ 

1 - Today, the working codes are called "BHE_3D_CODE" and "Wall_layer_work.jl", they contain the layer-discretisation model solving and plotting. And maybe some of the variables are "undefined" like _timer_ (or others) but you can just add a definition before instantiation of the model "soil_..._layer". The compilation is quite long for the first run. "BHE_3D_CODE" contain the computation with the variables layers of MRCTM. 

2 - You can have acces to the "include files" which contain _data's computation_ and the _components_ declaration. The "function_component.jl" contain the main features, the component declaration. So the focus have to be done on that. note : it uses the "wall_component" which is a serie of R-C-R-C-R (like a wall) and a component soil_temp which is a constant-voltage (temperature). 

3 - The component called "soil_MTRCM_var_pin_var_ver_2" contain the most updated work. So the focus have to be done on that. note : it uses the "soil_MRCTM" which is layers of MRCTM model, a component soil_temp which is a constant-voltage (temperature) and of course 2 fluid discretisation computation (with two different velocity sign).

![Alt Text](/anim_fps15.gif)

[deps]
CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
ColorSchemes = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
CompositeTypes = "b152e2b5-7a66-4b01-a709-34e65c35f657"
DataDrivenDiffEq = "2445eb08-9709-466a-b3fc-47e12bd697a2"
DifferentialEquations = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
DomainSets = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
IfElse = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
MethodOfLines = "94925ecb-adb7-4558-8ed8-f975c56a0bf4"
ModelingToolkit = "961ee093-0014-501f-94e3-6117800e7a78"
OrdinaryDiffEq = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
Setfield = "efcf1570-3423-57d1-acb7-fd33fddbac46"
SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"
StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
WGLMakie = "276b4fcb-3e11-5398-bf8b-a0c2d153d008"

