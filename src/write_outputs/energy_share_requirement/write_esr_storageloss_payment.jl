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
	write_esr_storagelosspayment(path::AbstractString, inputs::Dict, setup::Dict, EP::Model)

"""
function write_esr_storagelosspayment(path::AbstractString, inputs::Dict, setup::Dict, EP::Model)
    Z = inputs["Z"]
    nESR = inputs["nESR"]
    dfESRStoragelosspayment = DataFrame(Zone=1:Z, AnnualSum=zeros(Z))
    for i in 1:inputs["nESR"]
	    tempesrpayment = zeros(Z)
	    tempesrpayment = (value.(EP[:eStorageLossByZone])) .* inputs["dfESR"][:,Symbol("ESR_$i")] .* dual.(EP[:cESRShare][i])
	    if setup["ParameterScale"] == 1
	        tempesrpayment *= (ModelScalingFactor^2)
	    end
	    dfESRStoragelosspayment.AnnualSum .+= tempesrpayment
	    dfESRStoragelosspayment = hcat(dfESRStoragelosspayment, DataFrame([tempesrpayment], [Symbol("ESR_$i")]))
	end
    CSV.write(joinpath(path, "ESR_StoragelossPayment.csv"), dfESRStoragelosspayment)
    return dfESRStoragelosspayment
end
