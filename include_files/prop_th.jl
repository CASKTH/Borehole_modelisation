using CoolProp
# using Plots
# using Optim
# using LinearAlgebra
# using Roots

# http://www.coolprop.org/coolprop/LowLevelAPI.html#generating-input-pairs
# http://www.coolprop.org/_static/doxygen/html/namespace_cool_prop.html#a58e7d98861406dedb48e07f551a61efb
valid_pairs = [
    "QT_INPUTS" # Molar quality, Temperature in K.
    "PQ_INPUTS" # Pressure in Pa, Molar qualit
    "QSmolar_INPUTS" # Molar quality, Entropy in J/mol/
    "QSmass_INPUTS" # Molar quality, Entropy in J/kg/
    "HmolarQ_INPUTS" # Enthalpy in J/mol, Molar qualit
    "HmassQ_INPUTS" # Enthalpy in J/kg, Molar qualit
    "DmolarQ_INPUTS" # Density in mol/m^3, Molar qualit
    "DmassQ_INPUTS" # Density in kg/m^3, Molar qualit
    "PT_INPUTS" # Pressure in Pa, Temperature in 
    "DmassT_INPUTS" # Mass density in kg/m^3, Temperature in 
    "DmolarT_INPUTS" # Molar density in mol/m^3, Temperature in 
    "HmolarT_INPUTS" # Enthalpy in J/mol, Temperature in 
    "HmassT_INPUTS" # Enthalpy in J/kg, Temperature in 
    "SmolarT_INPUTS" # Entropy in J/mol/K, Temperature in 
    "SmassT_INPUTS" # Entropy in J/kg/K, Temperature in 
    "TUmolar_INPUTS" # Temperature in K, Internal energy in J/mo
    "TUmass_INPUTS" # Temperature in K, Internal energy in J/k
    "DmassP_INPUTS" # Mass density in kg/m^3, Pressure in P
    "DmolarP_INPUTS" # Molar density in mol/m^3, Pressure in P
    "HmassP_INPUTS" # Enthalpy in J/kg, Pressure in P
    "HmolarP_INPUTS" # Enthalpy in J/mol, Pressure in P
    "PSmass_INPUTS" # Pressure in Pa, Entropy in J/kg/
    "PSmolar_INPUTS" # Pressure in Pa, Entropy in J/mol/
    "PUmass_INPUTS" # Pressure in Pa, Internal energy in J/k
    "PUmolar_INPUTS" # Pressure in Pa, Internal energy in J/mo
    "HmassSmass_INPUTS" # Enthalpy in J/kg, Entropy in J/kg/
    "HmolarSmolar_INPUTS" # Enthalpy in J/mol, Entropy in J/mol/
    "SmassUmass_INPUTS" # Entropy in J/kg/K, Internal energy in J/k
    "SmolarUmolar_INPUTS" # Entropy in J/mol/K, Internal energy in J/mo
    "DmassHmass_INPUTS" # Mass density in kg/m^3, Enthalpy in J/k
    "DmolarHmolar_INPUTS" # Molar density in mol/m^3, Enthalpy in J/mo
    "DmassSmass_INPUTS" # Mass density in kg/m^3, Entropy in J/kg/
    "DmolarSmolar_INPUTS" # Molar density in mol/m^3, Entropy in J/mol/
    "DmassUmass_INPUTS" # Mass density in kg/m^3, Internal energy in J/k
    "DmolarUmolar_INPUTS" # Molar density in mol/m^3, Internal energy in J/mol.
]

function decompose_pairs_symbols(pairs)
    rad = map(x -> x[1:end-7], pairs)
    function pos_cap(st)
        p = fill(false, length(st))
        for i in 1:length(st)
            c = st[i]
            if c == uppercase(c)
                p[i] = true
            end
        end
        return p
    end
    function split_cap(st)
        p = pos_cap(st)
        c2 = findlast(p)
        return [st[1:c2-1] st[c2:end]]
    end
    combos = reduce(vcat, split_cap.(rad))
    not_molar = .!occursin.("molar", combos[:, 1]) .&& .!occursin.("molar", combos[:, 2])
    useful_pairs = combos[not_molar, :]
    # values = unique(useful_pairs)

    key_to_symbol = Dict("Q" => :χ, "P" => :p, "Hmass" => :h, "Umass" => :U, "T" => :T, "Smass" => :S, "Dmass" => :ρ)

    pair_symbols_in = map(x -> key_to_symbol[x], useful_pairs)
    return pair_symbols_in
end
ins = decompose_pairs_symbols(valid_pairs)
symbol_to_in_key = Dict(:ρ => "Dmass", :T => "T", :p => "P", :χ => "Q", :h => "Hmass", :U => "Umass", :S => "Smass")

outs = vcat(unique(ins), [:μ; :λ; :σ; :Pr; :p; :dTdh_P])

symbol_to_out_key = Dict(:ρ => "Dmass", :T => "T", :p => "P", :Pr => "Prandtl", :μ => "viscosity", :λ => "conductivity", :σ => "surface_tension", :χ => "Q", :h => "Hmass", :U => "Umass", :S => "Smass", :cp => "Cpmass")
out_keys = [
    "igas_constant"  #Ideal-gas constant.
    "imolar_mass"  #Molar mass.
    "iacentric_factor"  #Acentric factor.
    "irhomolar_reducing"  #Molar density used for the reducing state.
    "irhomolar_critical"  #Molar density used for the critical point.
    "iT_reducing"  #Temperature at the reducing state.
    "iT_critical"  #Temperature at the critical point.
    "irhomass_reducing"  #Mass density at the reducing state.
    "irhomass_critical"  #Mass density at the critical point.
    "iP_critical"  #Pressure at the critical point.
    "iP_reducing"  #Pressure at the reducing point.
    "iT_triple"  #Triple point temperature.
    "iP_triple"  #Triple point pressure.
    "iT_min"  #Minimum temperature.
    "iT_max"  #Maximum temperature.
    "iP_max"  #Maximum pressure.
    "iP_min"  #Minimum pressure.
    "idipole_moment"  #Dipole moment.
    "iT"  #Temperature.
    "iP"  #Pressure.
    "iQ"  #Vapor quality.
    "iTau"  #Reciprocal reduced temperature.
    "iDelta"  #Reduced density.
    "iDmolar"  #Mole-based density.
    "iHmolar"  #Mole-based enthalpy.
    "iSmolar"  #Mole-based entropy.
    "iCpmolar"  #Mole-based constant-pressure specific heat.
    "iCp0molar"  #Mole-based ideal-gas constant-pressure specific heat.
    "iCvmolar"  #Mole-based constant-volume specific heat.
    "iUmolar"  #Mole-based internal energy.
    "iGmolar"  #Mole-based Gibbs energy.
    "iHelmholtzmolar"  #Mole-based Helmholtz energy.
    "iHmolar_residual"  #The residual molar enthalpy.
    "iSmolar_residual"  #The residual molar entropy (as a function of temperature and density)
    "iGmolar_residual"  #The residual molar Gibbs energy.
    "iDmass"  #Mass-based density.
    "iHmass"  #Mass-based enthalpy.
    "iSmass"  #Mass-based entropy.
    "iCpmass"  #Mass-based constant-pressure specific heat.
    "iCp0mass"  #Mass-based ideal-gas specific heat.
    "iCvmass"  #Mass-based constant-volume specific heat.
    "iUmass"  #Mass-based internal energy.
    "iGmass"  #Mass-based Gibbs energy.
    "iHelmholtzmass"  #Mass-based Helmholtz energy.
    "iviscosity"  #Viscosity.
    "iconductivity"  #Thermal conductivity.
    "isurface_tension"  #Surface tension.
    "iPrandtl"  #The Prandtl number.
    "ispeed_sound"  #Speed of sound.
    "iisothermal_compressibility"  #Isothermal compressibility.
    "iisobaric_expansion_coefficient"  #Isobaric expansion coefficient.
    "iisentropic_expansion_coefficient"  #Isentropic expansion coefficient.
    "ifundamental_derivative_of_gas_dynamics"  #The fundamental derivative of gas dynamics.
    "ialphar"
    "idalphar_dtau_constdelta"
    "idalphar_ddelta_consttau"
    "ialpha0"
    "idalpha0_dtau_constdelta"
    "idalpha0_ddelta_consttau"
    "iBvirial"  #Second virial coefficient.
    "iCvirial"  #Third virial coefficient.
    "idBvirial_dT"  #Derivative of second virial coefficient with temperature.
    "idCvirial_dT"  #Derivative of third virial coefficient with temperature.
    "iZ"  #The compressibility factor Z = p*v/(R*T)
    "iPIP"  #The phase identification parameter of Venkatarathnam and Oellrich.
    "ifraction_min"  #The minimum fraction (mole, mass, volume) for incompressibles.
    "ifraction_max"  #The maximum fraction (mole,mass,volume) for incompressibles.
    "iT_freeze"  #The freezing temperature for incompressibles.
    "iGWP20"  #The 20-year global warming potential.
    "iGWP100"  #The 100-year global warming potential.
    "iGWP500"  #The 500-year global warming potential.
    "iFH"  #Fire hazard index.
    "iHH"  #Health hazard index.
    "iPH"  #Physical hazard index.
    "iODP"  #Ozone depletion potential (R-11 = 1.0)
    "iPhase"  #The phase index of the given state.
    "iundefined_parameter"  #The last parameter, so we can check that all parameters are described in DataStructures.cpp. 
]


refs = ["water", "R404A"]
handles = Dict([r => CoolProp.AbstractState_factory("HEOS", r) for r in refs])


#improve again: do not update if the handle point didn't change
# apply the same logic for the partial derivative!!!
for output in outs
    if [c == 'd' for c in string(output)] |> sum == 2 # partial derivative

        exp = quote #une seule fonction par grandeur, other wise it erases
            function $(Symbol(string(output) * "_"))(; fluid, kwargs...)
                try
                    var = Dict(kwargs)
                    i1 = collect(keys(var))[1]
                    v1 = var[i1]
                    i2 = collect(keys(var))[2]
                    v2 = var[i2]
                    # key_in = symbol_to_in_key[i1]*symbol_to_in_key[i2]*"_INPUTS"
                    key_in1 = symbol_to_in_key[i1]
                    key_in2 = symbol_to_in_key[i2]
                    key_in1 * key_in2 * "_INPUTS" ∈ valid_pairs ? key_in = key_in1 * key_in2 * "_INPUTS" : key_in = key_in2 * key_in1 * "_INPUTS"
                    hand = handles[fluid] |> Int32

                    basis = $(string(output))
                    posd1 = findfirst('d', basis)
                    posd2 = findlast('d', basis)
                    pos_ = findlast('_', basis)
                    var1 = basis[posd1+1:posd2-1]
                    var2 = basis[posd2+1:pos_-1]
                    var3 = basis[pos_+1:end]
                    #manuel thingy, cannot get 2p thingy
                    key_out1 = symbol_to_out_key[Symbol(var1)]
                    key_out2 = symbol_to_out_key[Symbol(var2)]
                    key_out3 = symbol_to_out_key[Symbol(var3)]

                    nested_exp = quote
                        try
                            if !(CoolProp.AbstractState_output($hand, $key_in1) == $v1 && CoolProp.AbstractState_output($hand, $key_in2) == $v2)
                                CoolProp.AbstractState_update($hand, $key_in, $v1, $v2)
                            end
                        catch
                            CoolProp.AbstractState_update($hand, $key_in, $v1, $v2)
                        end
                        return CoolProp.AbstractState_first_partial_deriv($hand, get_param_index($key_out1), get_param_index($key_out2), get_param_index($key_out3))
                    end
                    eval(nested_exp)
                catch
                    throw("Pair not supported. Keys : ")
                end
            end
        end
        eval(exp)
    else #normal value
        key_call_out = symbol_to_out_key[output]
        exp = quote #only one function per variable, otherwise it erases
            function $(Symbol(string(output) * "_"))(; fluid, kwargs...)
                try
                    var = Dict(kwargs)
                    ia = collect(keys(var))[1]
                    va = var[ia]
                    ib = collect(keys(var))[2]
                    vb = var[ib]

                    key_ina = symbol_to_in_key[ia]
                    key_inb = symbol_to_in_key[ib]

                    if key_ina * key_inb * "_INPUTS" ∈ valid_pairs
                        key_in = key_ina * key_inb * "_INPUTS"
                        key_in1 = key_ina
                        v1 = va
                        key_in2 = key_inb
                        v2 = vb
                    else
                        key_in = key_inb * key_ina * "_INPUTS"
                        key_in1 = key_inb
                        v2 = va
                        key_in2 = key_ina
                        v1 = vb
                    end
                    hand = handles[fluid] |> Int32

                    try
                        if !(CoolProp.AbstractState_output(hand, key_in1) ≈ v1 && CoolProp.AbstractState_output(hand, key_in2) ≈ v2)
                            CoolProp.AbstractState_update(hand, key_in, v1, v2)
                        else
                        end
                        return CoolProp.AbstractState_output(hand, $key_call_out)
                    catch
                        CoolProp.AbstractState_update(hand, key_in, v1, v2)
                        return CoolProp.AbstractState_output(hand, $key_call_out)
                    end
                catch
                    # println($output)
                    throw("Pair not supported")
                end
            end
        end
        eval(exp)
    end
end
