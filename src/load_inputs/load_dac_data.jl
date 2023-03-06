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
	load_dac_data!(setup::Dict, path::AbstractString, inputs::Dict)

Read input parameters related to dac
"""

function load_dac_data!(setup::Dict, path::AbstractString, inputs::Dict)
    filename ="dac.csv"

    dac_in = load_dataframe(joinpath(path, filename))

	# Add Resource IDs after reading to prevent user errors
	dac_in[!,:DAC_ID] = 1:length(collect(skipmissing(dac_in[!,1])))

	# Store DataFrame of dac input data for use in model
	inputs["dfDac"] = dac_in

    # Number of resources
	inputs["D"] = length(collect(skipmissing(dac_in[!,:DAC_ID])))

	# Set indices for internal use
	D = inputs["D"]   # Number of DAC resources 

    # Zones resources are located in
	zones = collect(skipmissing(dac_in[!,:Zone][1:inputs["D"]]))
	# Resource identifiers by zone (just zones in resource order + resource and zone concatenated)
	inputs["D_ZONES"] = zones
	#inputs["dfDac"][!,:DAC_ZONES] = dac_in["RESOURCES"] .* "_z" .* string.(zones)



    # the capacity of DAC is expressed as t CO2/h 
    # the unit of fixed cost is $/t CO2/have
    # the unit of non-fuel VOM is $/t CO2
    # the unit of heat consumption is MMBTU/t CO2
    # the unit of electricity consumption is MWh/t CO2

    if setup["ParameterScale"] == 1  # Parameter scaling turned on - adjust values of subset of parameter values, t CO2 should be converted to kt CO2
        # I dont think we have exiting capacity for DAC..so I will just do the conversion for Fixed Cost, VOM, heat, electricity consumption
        inputs["dfDac"][!, :Fix_Cost_per_CO2perHr_yr] = dac_in[!, :Fix_Cost_per_CO2perHr_yr]/ModelScalingFactor

        inputs["dfDac"][!, :Var_OM_Cost_per_CO2] = dac_in[!, :Var_OM_Cost_per_CO2]/ModelScalingFactor

        inputs["dfDac"][!, :Heat_MMBTU_per_CO2_metric_ton] = dac_in[!, :Heat_MMBTU_per_CO2_metric_ton]/ModelScalingFactor

        inputs["dfDac"][!, :Electricity_MWh_per_CO2_metric_ton] = dac_in[!, :Electricity_MWh_per_CO2_metric_ton]/ModelScalingFactor
    end

    println(filename * " Successfully Read!")

   # return inputs
end




