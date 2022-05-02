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

Function for writing the diferent capacities for the FLECCS technologies (starting capacities or, existing capacities, retired capacities, and new-built capacities).
"""
function write_capacity_fleccs(path::AbstractString, inputs::Dict, setup::Dict, EP::Model)
	# Capacity decisions
	dfGen_ccs = inputs["dfGen_ccs"]
	FLECCS_ALL = inputs["FLECCS_ALL"]
	N_F = inputs["N_F"]
	COMMIT_ccs  =inputs["COMMIT_CCS"]
	# the number of rows for FLECCS generator 
	#n = length(dfGen_ccs[!,"Resource"])/length(N_F)

    # the number of subcompoents 
	N = length(N_F)

	capFLECCS = zeros(size(dfGen_ccs[!,"Resource"]))
    #reshape(zeros(size(dfGen_ccs[!,"Resource"])),length(FLECCS_ALL),length(N_F))
	
	for y in FLECCS_ALL
	    for i in inputs["NEW_CAP_FLECCS"]
		    for i in COMMIT_ccs
			    capFLECCS[(y-1)*N + i] = value(EP[:vCAP_FLECCS][y, i])* dfGen_ccs[(dfGen_ccs[!,:FLECCS_NO].==i),:Cap_Size][y]
		    end
		end

		for i in setdiff(N_F,COMMIT_ccs)
			capFLECCS[(y-1)*N + i] = value(EP[:vCAP_FLECCS][y, i])
		end
	end

	retcapFLECCS = zeros(size(dfGen_ccs[!,"Resource"]))

	
	for y in FLECCS_ALL
	    for i in inputs["RET_CAP_FLECCS"]
		    for i in COMMIT_ccs
			     retcapFLECCS[(y-1)*N + i] = value(EP[:vRETCAP_FLECCS][y, i])* dfGen_ccs[(dfGen_ccs[!,:FLECCS_NO].==i),:Cap_Size][y]
		    end
		end

		for i in setdiff(N_F,COMMIT_ccs)
			retcapFLECCS[(y-1)*N + i] = value(EP[:vRETCAP_FLECCS][y, i])
		end
	end




	EndCapFLECCS = zeros(size(dfGen_ccs[!,"Resource"]))
	for y in FLECCS_ALL
		for i in N_F
		    EndCapFLECCS[(y-1)*N+i] = value.(EP[:eTotalCapFLECCS])[y,i]
		end
	end




	dfCapFLECCS = DataFrame(
		Resource = dfGen_ccs[!,"Resource"],
		Zone = dfGen_ccs[!,:Zone], R_ID = dfGen_ccs[!,:R_ID],
		StartCap = dfGen_ccs[!,:Existing_Cap_Unit],
		RetCap = retcapFLECCS[:],
		NewCap = capFLECCS[:],
		EndCap = EndCapFLECCS,
	)
	if setup["ParameterScale"] ==1
		dfCapFLECCS.StartCap = dfCapFLECCS.StartCap * ModelScalingFactor
		dfCapFLECCS.RetCap = dfCapFLECCS.RetCap * ModelScalingFactor
		dfCapFLECCS.NewCap = dfCapFLECCS.NewCap * ModelScalingFactor
		dfCapFLECCS.EndCap = dfCapFLECCS.EndCap * ModelScalingFactor
	end

	total = DataFrame(
			Resource = "Total", Zone = "n/a", R_ID = "n/a", 
			StartCap = sum(dfCapFLECCS[!,:StartCap]), RetCap = sum(dfCapFLECCS[!,:RetCap]),
			NewCap = sum(dfCapFLECCS[!,:NewCap]), EndCap = sum(dfCapFLECCS[!,:EndCap]),
		)

	dfCap_FLECCS = vcat(dfCapFLECCS, total)
	CSV.write(joinpath(path,"capacity_FLECCS.csv"), dfCap_FLECCS)

####################################### plot net capacity


	capFLECCS_net = zeros(size(FLECCS_ALL))
    #reshape(zeros(size(dfGen_ccs[!,"Resource"])),length(FLECCS_ALL),length(N_F))
	
	for y in FLECCS_ALL
		capFLECCS_net[y] = value(EP[:vCAP_FLECCS_tx][y])
	end

    BOP_id = inputs["BOP_id"]
	j = BOP_id*FLECCS_ALL

	dfCapFLECCS_net = DataFrame(
		Resource = dfGen_ccs[!,"Resource"][j],
		Zone = dfGen_ccs[!,:Zone][j], R_ID = dfGen_ccs[!,:R_ID][j],
		EndCap = capFLECCS_net,
	)
	if setup["ParameterScale"] ==1
		dfCapFLECCS_net.EndCap = dfCapFLECCS.EndCap * ModelScalingFactor
	end

	total_net = DataFrame(
			Resource = "Total", Zone = "n/a", R_ID = "n/a", 
			EndCap = sum(dfCapFLECCS_net[!,:EndCap]),
		)

	dfCapFLECCS_net = vcat(dfCapFLECCS_net, total_net)

	CSV.write(joinpath(path,"Netcapacity_FLECCS.csv"), dfCapFLECCS_net)





	return dfCap_FLECCS
end
