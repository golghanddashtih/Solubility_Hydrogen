using Plots

# --------------------------------------------------
# Temperature range (K)
# --------------------------------------------------
Tvals = collect(275.0:1.0:350)
n = length(Tvals)

# --------------------------------------------------
# Fixed pressures (Pa)
# --------------------------------------------------
  pressures_MPa = [0.1, 1.0, 10.0, 50.0, 100.0]
#  pressures_MPa = [0.101325]   # MPa
pressures = pressures_MPa .* 1e6   # MPa → Pa

# --------------------------------------------------
# Fixed inputs (size based on T)
# --------------------------------------------------
z   = [ [0.5, 0.5] for _ in 1:n ]   # composition
Sal = fill(0.0, n)                 # salinity

# --------------------------------------------------
# Plot setup
# --------------------------------------------------
plt = plot(
    xlabel = "Temperature (K)",
    ylabel = "Mole Fraction of H₂ in Liquid Phase",
    title  = "Solubility of H₂ in Water",
    legend = :topright,
    lw = 2
)

# --------------------------------------------------
# Loop over pressures
# --------------------------------------------------
for Pval in pressures
    P = fill(Pval, n)

    solubility = compute_solubility(
        model,
        Tvals,
        P,
        z;
        Salinity = Sal
    )

    plot!(
        plt,
        Tvals,
        solubility.X[1, :],   # H₂ in liquid
        label = "P = $(Pval / 1e6) MPa"
    )
end

# --------------------------------------------------
# Display plot
# --------------------------------------------------
display(plt)