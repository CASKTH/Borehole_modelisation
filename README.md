# Borehole_modelisation
Project dedicated to borehole modilisation and command using Julia 
---

_Understand the model_ 

The model is a Acausal component based model containing : 
- The discrete convection-diffusion component giving a temperature discretisation through time and space : u(t)[x]
- The layer instantiation, function of n size of u(t)[x]
- Ploting of the solution


_How to understand the code_ 

1 - Today, the working code is called "Wall_layer_work.jl" and contain the layer-discretisation model solving and plotting. And maybe some of the variables are "undefined" like _timer_ (or others) but you can just add a definition before instantiation of the model "soil_n_var_layer". The compilation is quite long for the first run. 

2 - You can have acces to the "include files" which contain _data's computation_ and the _components_ declaration. The "function_component.jl" contain the main features, the component declaration. At the moment I'm writting this Readme.md the component called "soil_n_var_layer" contain the most updated work. So the focus have to be done on that. note : it uses the "wall_component" which is a serie of R-C-R-C-R (like a wall) and a component soil_temp which is a constant-voltage (temperature). 

3 - The code "function_component.jl" contain also @component which are prepared to be updated like the _"soil_MTRCM_var_pin_var"_, more updates will be done. 
