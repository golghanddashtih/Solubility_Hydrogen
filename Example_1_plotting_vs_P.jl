using Plots

# --------------------------------------------------
# Pressure range (Pa)
# --------------------------------------------------
p = collect(0.1e6:0.1e6:100e6)   # 0.1–50 MPa
n = length(p)

# --------------------------------------------------
# Fixed inputs
# --------------------------------------------------
z   = [ [0.5, 0.5] for _ in 1:n ]   # composition
Sal = fill(0.0, n)                  # salinity

# --------------------------------------------------
# Temperatures to evaluate (K)
# --------------------------------------------------
temperatures = [298.0, 320.0, 350.0, 370.0, 400.0]

# --------------------------------------------------
# Plot setup
# --------------------------------------------------
plt = plot(
    xlabel = "Pressure (MPa)",
    ylabel = "Mole Fraction of H₂ in Liquid Phase",
    title  = "Solubility of H₂ in Water",
    legend = :topright
)

# --------------------------------------------------
# Loop over temperatures
# --------------------------------------------------
for Tval in temperatures
    T = fill(Tval, n)

    solubility = compute_solubility(
        model,
        T,
        p,
        z;
        Salinity = Sal
    )

    plot!(
        plt,
        solubility.Pressure ./ 1e6,     # Pa → MPa
        solubility.X[1, :],
        label = "T = $(Int(Tval)) K",
        lw = 2
    )
end

# --------------------------------------------------
# Display plot
# --------------------------------------------------
display(plt)