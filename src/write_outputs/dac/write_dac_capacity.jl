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
	write_capacity(path::AbstractString, inputs::Dict, setup::Dict, EP::Model))

Function for writing the diferent capacities for the different generation technologies (starting capacities or, existing capacities, retired capacities, and new-built capacities).
"""
function write_dac_capacity(path::AbstractString, inputs::Dict, setup::Dict, EP::Model)
    # scale factor
	scale_factor = setup["ParameterScale"] == 1 ? ModelScalingFactor : 1
	# Capacity decisions
	dfDac = inputs["dfDac"]
	#MultiStage = setup["MultiStage"]
	NEW_CAP = dfDac[dfDac[!,:NEW_CAP] .==1,:R_ID]
	DAC_COMMIT = dfDac[dfDac[!,:DAC_COMMIT] .==1,:R_ID]
	capdDac = zeros(size(dfDac[!,:Resource]))
	for i in  NEW_CAP
		capdDac[i] = value(EP[:vCAP_DAC][i])
	end


	dfCapDac = DataFrame(
		Resource = dfDac[!,:Resource],
		Zone = dfDac[!,:Zone],
		StartCap = dfDac[!,:Existing_Cap_CO2],
		NewCap = capdDac[:],
		EndCap = capdDac[:]
	)

	dfCapDac.StartCap = dfCapDac.StartCap * scale_factor
	dfCapDac.NewCap = dfCapDac.NewCap * scale_factor
	dfCapDac.EndCap = dfCapDac.EndCap * scale_factor

	#dfCapDac = vcat(dfCapDac, total)
	CSV.write(joinpath(path, "capacity_dac.csv"), dfCapDac)


    # write the dual capex cost if fix cost = 0
	CapexdDac = zeros(size(dfDac[!,:Resource]))

	if !isempty(dual.(EP[:cMaxCap_DAC]))
		for i in NEW_CAP
		    CapexdDac[i] = dual.(EP[:cMaxCap_DAC][i])*scale_factor
	    end
	end

	dfCapexDac = DataFrame(
		Resource = dfDac[!,:Resource],
		Zone = dfDac[!,:Zone],
		Capex = -CapexdDac[:]
	)




	if 0 in dfDac.Fix_Cost_per_CO2perHr_yr
		CSV.write(joinpath(path, "capex_dual_dac.csv"), dfCapexDac)
	end


	# write dac cost

	dfDac = inputs["dfDac"]
    # Number of time steps (hours)

	dfCost = DataFrame(Costs = ["cTotal", "cFix", "cVar", "cCO2_seq", "cCO2_tax"])
	cVar = (!isempty(inputs["dfDac"]) ? value(EP[:eCTotalVariableDAC])*scale_factor^2 : 0.0)
	cFix = (!isempty(inputs["dfDac"]) ? value(EP[:eTotalCFixedDAC])*scale_factor^2 : 0.0)
    cCO2_seq =  (!isempty(inputs["dfDac"]) ? value(EP[:eCTotalCO2TS])*scale_factor^2 : 0.0)
	cCO2_tax =  ((setup["CO2Tax"]  > 0)  ? value.(EP[:eTotalCCO2TaxDAC])*scale_factor^2 : 0)

	cDacTotal = 0 
	cDacTotal += (cVar + cFix + cCO2_seq + cCO2_tax)

	dfCost[!,Symbol("Total")] = [cDacTotal, cFix, cVar, cCO2_seq, cCO2_tax]

	CSV.write(joinpath(path, "Dac_costs.csv"), dfCost)




	
end
