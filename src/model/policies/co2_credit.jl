"""
GenX: An Configurable Capacity Expansion Model
Copyright (C) 2021,  Massachusetts Institute of Technology
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
A complete copy of the GNU General Public License v2 (GPLv2) is available
in LICENSE.txt.  Users uncompressing this from an archive may not have
received this license file.  If not, see <http://www.gnu.org/licenses/>.
"""

@doc raw"""
	co2_cap(EP::Model, inputs::Dict, setup::Dict)

	hat the generator-side rate-based constraint can be used to represent a fee-rebate (``feebate'') system: the dirty generators that emit above the bar ($\epsilon_{z,p,gen}^{maxCO_2}$) have to buy emission allowances from the emission regulator in the region $z$ where they are located; in the same vein, the clean generators get rebates from the emission regulator at an emission allowance price being the dual variable of the emissions rate constraint.
"""
function co2_credit!(EP::Model, inputs::Dict, setup::Dict)

    println("C02 Credit Module")

    dfGen = inputs["dfGen"]
    SEG = inputs["SEG"]  # Number of lines
    G = inputs["G"]     # Number of resources (generators, storage, DR, and DERs)
    T = inputs["T"]     # Number of time steps (hours)
    Z = inputs["Z"]     # Number of zones
    ### Expressions ###


    # determine the CO2 credit by mutiplying the capture co2 and co2 credit. This is a negative value.
    @expression(EP, ePlantCCO2Credit[y = 1:G], inputs["dfCO2Credit"][dfGen[y, :Zone], :CO2Credit] * sum(inputs["omega"][t] * EP[:eEmissionsCaptureByPlant][y, t] for t in 1:T))

    @expression(EP, eZonalCCO2Credit[z = 1:Z], EP[:vZERO] + sum(EP[:ePlantCCO2Credit][y] for y in dfGen[(dfGen[!, :Zone].==z), :R_ID]))

    # sum cross the zones
    @expression(EP, eTotalCCO2Credit, sum(eZonalCCO2Credit[z] for z in 1:Z))
    # add to objective function

    if setup["FLECCS"] > 0
        G_F = inputs["G_F"] # Number of FLECCS generator
        dfGen_ccs = inputs["dfGen_ccs"]
        @expression(EP, ePlantCCO2CreditFLECCS[y = 1:G_F], sum(EP[:eCO2CaptureByPlantFLECCS][y, t] for t in 1:T) * inputs["dfCO2Credit"][dfGen[y, :Zone], "CO2Credit"])
        # Sum CO2 Credit to zonal level
        @expression(EP, eZonalCCO2CreditFLECCS[z = 1:Z],  sum(EP[:ePlantCCO2CreditFLECCS][y] for y in unique(dfGen_ccs[(dfGen_ccs[!, :Zone].==z), :R_ID])))
        # Sum CO2 Credit to system level
        @expression(EP, eTotalCCO2CreditFLECCS, sum(EP[:eZonalCCO2CreditFLECCS][z] for z in 1:Z))
        add_to_expression!(EP[:eObj],  -1, EP[:eTotalCCO2CreditFLECCS])
    end    

    add_to_expression!(EP[:eObj], -1, EP[:eTotalCCO2Credit])

end
