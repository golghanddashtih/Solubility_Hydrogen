
####Packages
using Pkg
Pkg.add("Clapeyron")
using Clapeyron
####
#######Solubility_High_Pressure
function compute_solubility(model::Clapeyron.EoSModel,T_values::AbstractVector{<:Real}, P_values::AbstractVector{<:Real}, 
    composition_feed::Vector{Vector{Float64}}; 
    Salinity::AbstractVector{<:Real}= zeros(length(T_values)),  
    Pref::Real    = 0.1e6, 
   )
# --- Extract component names from model ---
    component_names = model.components 
# --- Normalize inputs to vectors (supports vector, range, or single number) ---
    T_vec = (T_values isa AbstractVector) ? collect(T_values) :
            (T_values isa Base.AbstractRange) ? collect(T_values) : [T_values]
    P_vec = (P_values isa AbstractVector) ? collect(P_values) :
            (P_values isa Base.AbstractRange) ? collect(P_values) : [P_values]

    nT = length(T_vec); nP = length(P_vec)

    nc=length(component_names)
    nfc = size(composition_feed)[1]
    nfr = size(composition_feed[1])[1]
    # println("→ Number of components: $nc")
    # println("→ Number of rows in composition_feed: $nfr")
    # println("→ Number of columns in composition_feed: $nfc")
    #### Size Checks ####
    @assert(nT == nP, "Error: nT and nP must be equal.")
    @assert(nc == nfr, "Error: Number of components in component_names and composition_feed must match.")
    @assert(nfc==nT, "Error: Number of columns in composition_feed must match number of temperature points.")
    #####
    #### Activity Coefficient due to Salinity Calculation
    ActSalt= zeros(nc, nT)
    A=0.3027532742783516
    B=0.00025474496618199763
    c=1.711e-5  

    row = exp.((A .- B .* T_vec) .* Salinity .* c)   # length nT
    ActSalt .= repeat(row', nc, 1)                   # nc × nT

    # println("→ Activity coefficient of H2 due to salinity calculated.",ActSalt )
    ### Preallocate fields ---  
    FugCoefV   = zeros(nc, nT)
    kH_T_Perry    = zeros(nc, nT)
    x_Perry_correct= zeros(nc, nT)
    PF_calc = zeros(nc, nT)
    x_tot_species= zeros(nT)
    ######
    for (i, comp) in enumerate(component_names)
        if comp=="water"||comp=="H2O"||comp=="H₂O"
            continue  # Skip water for now
        else
            Vbar= molar_liquid_volume(comp)
        
            # --- Build model and constants ---
            # model = _build_model_water(model_name)
            R     = Rgas(model)
            for j in 1:nT
                P = P_vec[j]  # Pa
                P_MPa = P * 1e-6    # MPa for the Li PF correlation
                T = T_vec[j]    # K
                z_feed = composition_feed[j]
                # println("Composition is ", z_feed)
                # Henry constant (Perry) — depends on T only; fill per cell for simplicity
                kH_T_Perry[i, j] = henry_constant(comp, T)

                # TP flash to get phase compositions (defensive try/catch)
                
                flash = try
                    Clapeyron.tp_flash2(model, P, T, z_feed )
                catch
                    
                    # Clapeyron.tp_flash(model, (p = P, T = T, z = z_feed))
                    Clapeyron.tp_flash(model,  P, T , z_feed )
                end

                x = flash.compositions[1]   # liquid
                y = flash.compositions[2]   # vapor
                # println("Composition is ", y)
                # Vapor-phase fugacity coefficients
                phi_v = fugacity_coefficient(model, P, T, y; phase = :v)
                FugCoefV[i,j]  = phi_v[i]
            # Perry-based x_H2 estimate and correction
                PF_calc[i,j]= Vbar*(P-Pref)/(Rgas(model) * T)
                x_Perry         = y[i] .* P./ kH_T_Perry[i,j] ./ 1e6
                x_Perry_correct[i,j] = x_Perry .* FugCoefV[i,j] ./ exp.(PF_calc[i,j]) ./ ActSalt[i,j]
            # Updating all species mole fraction except water
                x_tot_species[j]+= x_Perry_correct[i,j]
            end
        end
    end
    #Updating Water Mole Fraction in water phase 
    for (i, comp) in enumerate(component_names)
        if comp=="water"||comp=="H2O"||comp=="H₂O"
            for j in 1:nT
             x_Perry_correct[i,j]= 1 - x_tot_species[j]
            end
        else
                continue  # Skip
        end
    end
    # --- Assemble output (NamedTuple of matrices) ---
            mat_surface = (
                Pressure        = P_vec ,   # Pa
                Temperature     = T_vec,           # K
                FugCoefV     = FugCoefV,
                PF_calc         = PF_calc,
                PoyntingFactor  = exp.(PF_calc),
                kH_T_Perry      = kH_T_Perry,
                X    = x_Perry_correct
            )
        println("→ Solubility calculation completed and is equal to", mat_surface.X)
        return mat_surface
end
#######
using Clapeyron
###############Building Models
 # --- Build the model from name ---
    function _build_model_water(name)
        s = Symbol(name)
        lname = lowercase(String(s))
        if lname == "pr3"
            return PR(["hydrogen", "water"]; translation = PenelouxTranslation)
        elseif lname == "pr"
            return PR(["hydrogen", "water"])
        elseif lname == "srk3"
            return SRK(["hydrogen", "water"]; translation = PenelouxTranslation)
        elseif lname == "srk"
            return SRK(["hydrogen", "water"])
        elseif lname == "vdw"
            return vdW(["hydrogen", "water"])
        elseif lname == "gerg2008" || lname == "gerg"
            return GERG2008(["hydrogen", "water"])
        else
            error("Unknown model_name = $(model_name). Use :PR, :PR3, :SRK, :SRK3, :vdW, or :GERG2008.")
        end
    end
    ########
##### Henry's Law Constant Function (Perry's Method)
function henry_constant(gasname::AbstractString, T::Float64)
    # Default lower temperature bound
    Tfloor = 273.0
    Tcap = 373.0  # default upper bound

    # Initialize coefficients
    A = B = C = D = 0.0

    # Assign coefficients based on gas name
    gas = lowercase(strip(gasname))

    if gas == "acetylene"
        A, B, C, D = -156.51, 8160.2, 21.403, 0.0
        Tfloor, Tcap = 274.0, 343.0

    elseif gas == "carbon dioxide"
        A, B, C, D = -159.854, 8741.68, 21.6694, -1.10261e-3
        Tcap = 353.0

    elseif gas == "carbon monoxide"
        A, B, C, D = -171.764, 8296.9, 23.3376, 0.0
        Tcap = 353.0

    elseif gas == "ethane"
        A, B, C, D = -250.812, 12695.6, 34.7413, 0.0
        Tfloor, Tcap = 275.0, 323.0

    elseif gas == "ethylene"
        A, B, C, D = -153.027, 7965.2, 20.5248, 0.0
        Tfloor, Tcap = 287.0, 346.0

    elseif gas == "helium"
        A, B, C, D = -105.9768, 4259.62, 14.0094, 0.0
        Tcap = 348.0

    elseif gas == "hydrogen"
        A, B, C, D = -125.939, 5528.45, 16.8893, 0.0
        Tcap = 345.0

    elseif gas == "methane"
        A, B, C, D = -338.217, 13282.1, 51.9144, -0.0425831
        Tcap = 523.0

    elseif gas == "nitrogen"
        A, B, C, D = -181.587, 8632.12, 24.7981, 0.0
        Tcap = 350.0

    elseif gas == "oxygen"
        A, B, C, D = -171.2542, 8319.24, 23.24323, 0.0
        Tcap = 333.0

    else
        error("Gas name '$gasname' not recognized.")
    end

    # Henry’s law constant (Perry's equation form)
    H = 1 / exp(A + (B / T) + C * log(T) + D * T)*0.101325 # in atm to MPa converiosn now MPa/Mole fraction

    # Temperature validity check
    # if T > Tfloor && T < Tcap
    #     println("→ Temperature within valid range ($(Tfloor)–$(Tcap) K).")
    # else
    #     println("⚠ Temperature outside recommended range ($(Tfloor)–$(Tcap) K).")
    # end

    return H
end
###########
##Molar Liquid Volume Function (m3/mol)
function molar_liquid_volume(component::AbstractString)
    comp = lowercase(strip(component))

    if comp == "hydrogen"
        return 1.9954570665858385e-5  # m3/mol
    elseif comp == "nitrogen"
        return 3.311e-5  # m3/mol
    elseif comp == "methane"
        return 3.2e-5  # m3/mol
    elseif comp == "carbon dioxide"
        return 3.39e-5  # m3/mol
    # elseif comp == "water"
    #     return 1.8031534561120732e-5  # m3/mol
    else
        error("Component '$component' not recognized for molar liquid volume.")
    end
end