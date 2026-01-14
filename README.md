
<p align="center">
  <img src="Logo.png" alt="Project Logo" width="200"/>
</p>

# Solubility_Hydrogen
This code is developed in Julia Programming Language and is to compute the extent of solubility of hydrogen in water and brine under varying conditions of pressure, temperature and salinity. The accuracy has been tested across existing experimental data.
To use this file you can use the following simple codes:
Units are : Temperature, K; Pressure, Pa; Salinity, ppm or mg/l
#########Include The code
include("Solubility.jl")
######Define the model
model=PR(["hydrogen","water"]; translation=PenelouxTranslation)
model=SRK(["hydrogen","water"]) 
model=GERG2008(["hydrogen","water"])
########## Compute the solubility
compute_solubility(model,[298.0], [1.0e6], [[0.5,0.5]]; Salinity=[100000])
##############

You can give arrays of conditions as well. The model is verified vs experimental data for Hydrogen solubility in water and brine.
