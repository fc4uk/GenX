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
	write_energy_revenue(path::AbstractString, inputs::Dict, setup::Dict, EP::Model, dfPower_FLECCS::DataFrame, dfPrice::DataFrame, dfCharge::DataFrame)

Function for writing energy revenue from the different generation technologies.
"""
function write_energy_revenue_fleccs(path::AbstractString, inputs::Dict, setup::Dict, EP::Model, dfPower_FLECCS::DataFrame, dfPrice::DataFrame)
	
	dfGen_ccs = inputs["dfGen_ccs"]
	T = inputs["T"]     # Number of time steps (hours)
	# fleccs generators
	FLECCS_ALL = inputs["FLECCS_ALL"]
	G_F = inputs["G_F"]
	N_F = inputs["N_F"]
	N = length(N_F)

    # the index for BOP
	i = inputs["BOP_id"]
     # if we are testing mutiple zones or mutiple cases..
	j = inputs["BOP_id"]*FLECCS_ALL
	# dfEnergyRevenue = DataFrame(Resource = inputs["RESOURCES"], Zone = dfGen_ccs[!,:Zone], AnnualSum = Array{Union{Missing,Float32}}(undef, G))
	# the price is already US$/MWh, and dfPower_FLECCS and dfCharge is already in MW, so no scaling is needed
	dfEnergyRevenue = DataFrame( Resource = dfGen_ccs[!,"Resource"][j], Zone = dfGen_ccs[!,:Zone][j],R_ID = dfGen_ccs[!,:R_ID][j], AnnualSum = Array{Union{Missing,Float32}}(undef, G_F), )
	# initiation
	dfEnergyRevenue_FLECCS1 = (DataFrame([[names(dfPower_FLECCS)]; collect.(eachrow(dfPower_FLECCS))], [:column; Symbol.(axes(dfPower_FLECCS, 1))])[4:T+3,i+1] .*
	DataFrame([[names(dfPrice)]; collect.(eachrow(dfPrice))], [:column; Symbol.(axes(dfPrice, 1))])[2:T+1,dfPower_FLECCS[i,:][:Zone]+1].*
	inputs["omega"])
    # if we have mutiple NGCC-CCS with differernt parameters or setting
	if length(j) >1
		for i in 2:length(j)
			dfEnergyRevenue_FLECCS1 = hcat(dfEnergyRevenue_FLECCS1,(DataFrame([[names(dfPower_FLECCS)]; collect.(eachrow(dfPower_FLECCS))], [:column; Symbol.(axes(dfPower_FLECCS, 1))])[4:T+3,j[i]+1] .*
			DataFrame([[names(dfPrice)]; collect.(eachrow(dfPrice))], [:column; Symbol.(axes(dfPrice, 1))])[2:T+1,dfPower_FLECCS[j[i],:][:Zone]+1].*
			inputs["omega"]))
		end
	end
			

	dfEnergyRevenue_FLECCS = hcat(dfEnergyRevenue, DataFrame(dfEnergyRevenue_FLECCS1', :auto))


	for i in 1:length(j)
		dfEnergyRevenue_FLECCS[!,:AnnualSum][i] = sum(dfEnergyRevenue_FLECCS1[:,i])
	end




	dfEnergyRevenue_FLECCS_annualonly = dfEnergyRevenue_FLECCS[!,1:4]
	CSV.write(joinpath(path, "EnergyRevenue_FLECCS.csv"), dfEnergyRevenue_FLECCS_annualonly)
	return dfEnergyRevenue_FLECCS
end
