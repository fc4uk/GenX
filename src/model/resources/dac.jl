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
	DAC
"""

function dac!(EP::Model, inputs::Dict, setup::Dict)

	println("DAC module")

	dfDac = inputs["dfDac"]

	T = inputs["T"]     # Number of time steps (hours)
	Z = inputs["Z"]     # Number of zones
	D = inputs["D"] 	# Number of DAC

    START_SUBPERIODS = inputs["START_SUBPERIODS"]
	INTERIOR_SUBPERIODS = inputs["INTERIOR_SUBPERIODS"]
	hours_per_subperiod = inputs["hours_per_subperiod"]

    # load dac input 
    dfDac = inputs["dfDac"]
    DAC_ID = dfDac[!,:DAC_ID]
    G_DAC = length(collect(skipmissing(dfDac[!,:R_ID])))

    # fuel cost
    fuel_cost = inputs["fuel_costs"] 
    fuel_CO2 = inputs["fuel_CO2"]



    #decision variables
    # vCO2_DAC: the amount of capture captured by DAC facility, metric ton CO2/h. 
    # vDCAP: the capacity of DAC facility, metric ton CO2/h
    @variables(EP, begin
    vCO2_DAC[y in DAC_ID,t = 1:T] >= 0         
    vCAP_DAC[y in DAC_ID] >= 0    
    end)


    # DAC facility will either take heat (from natural gas) and electricity (from grid) as input, and produce negative co2 as write_outputs
    # heat consumption from heat resources, MMBTU/t CO2
    @expression(EP, eDAC_heat_consumption[y in DAC_ID, t = 1:T], vCO2_DAC[y,t] * dfDac[y,:Heat_MMBTU_per_CO2_metric_ton])
    # the use of heat resources (e.g., natural gas) may result in additional CO2 emissions as DAC may not capture the CO2 from the heat resources, the co2 content from heat resouces times (1 - capture rate) = emitted co2 from heat resource
    @expression(EP, eDAC_heat_CO2[y in DAC_ID, t = 1:T], eDAC_heat_consumption[y,t]* fuel_CO2[dfDac[y,:DAC_heat_resource]] * (1 - dfDac[y, :DAC_Fuel_Capture_Rate] ))
    # the electricity consumption for DAC, MWh/t CO2
    @expression(EP, eDAC_power[y in DAC_ID, t = 1:T], vCO2_DAC[y,t] * dfDac[y,:Electricity_MWh_per_CO2_metric_ton])
    # the power used for DAC is draw from the grid, so we have a power balance equation
	@expression(EP, ePowerBalanceDAC[t=1:T, z=1:Z],
		sum(eDAC_power[y,t] for y in intersect(dfDac[dfDac[!,:Zone].==z,:][!,:DAC_ID])))
		
	EP[:ePowerBalance] = EP[:ePowerBalance] - ePowerBalanceDAC

    # finally, the CO2 captured by DAC during each hr should be less than the capacity
	@constraint(EP, [y in DAC_ID, t=1:T],vCO2_DAC[y,t] <= vCAP_DAC[y])

    #---------------------------------- add up cost ---------------------------------------
    # Fixed Cost
    # Combine CAPEX and FOM into annulized Fixed Cost
    # Fixed cost for a DAC y 
	@expression(EP, eCFixed_DAC[y in DAC_ID], dfDac[y,:Fix_Cost_per_CO2perHr_yr] * vCAP_DAC[y])
	# total fixed costs for all the DAC
	@expression(EP, eTotalCFixedDAC, sum(eCFixed_DAC[y] for y in DAC_ID))
	EP[:eObj] += eTotalCFixedDAC


    # Variable cost
    omega = inputs["omega"]
    # the total variable cost (heat cost + non-fuel vom cost) for DAC y at time t, $/t CO2  
    @expression(EP, eCDAC_Variable[y in DAC_ID, t = 1:T],  (eDAC_heat_consumption[y,t] * fuel_cost[dfDac[y,:DAC_heat_resource]][t] + dfDac[y,:Var_OM_Cost_per_CO2]*vCO2_DAC[y,t]))

    # Cost associated with variable cost for DAC for the whole year
    
    @expression(EP, eCTotalVariableDACT[y in DAC_ID], sum(omega[t] * eCDAC_Variable[y,t] for t in 1:T ))

    # Total vairable cost for all the DAC at all the time
    @expression(EP, eCTotalVariableDAC, sum(eCTotalVariableDACT[y] for y in DAC_ID ))

    EP[:eObj] += eCTotalVariableDAC

    # skip the unit commitment for DAC for now, will implement constraint formulation later...
    # unit commitment variables for DAC
    

    if setup["UCommit"] > 0
        DAC_COMMIT = dfDac[dfDac[!,:DAC_COMMIT] .== 1, :R_ID]
        @variables(EP, begin
            vCOMMIT_DAC[y in DAC_COMMIT, t=1:T] >= 0 # commitment status
            vSTART_DAC[y in DAC_COMMIT, t=1:T] >= 0 # startup
            vSHUT_DAC[y in DAC_COMMIT, t=1:T] >= 0 # shutdown
        end)

         # set the unit commitment variables to integer is UC = 1
        for y in DAC_COMMIT
		    if setup["UCommit"] == 1
			    set_integer.(vCOMMIT_DAC[y,:])
			    set_integer.(vSTART_DAC[y,:])
			    set_integer.(vSHUT_DAC[y,:])
			    set_integer.(EP[:vCAP_DAC][y])
		    end
	    end


        @constraints(EP, begin
		    [y in DAC_COMMIT, t=1:T], vCOMMIT_DAC[y,t] <= EP[:vCAP_DAC][y] / dfDac[y,:Cap_Size]
		    [y in DAC_COMMIT, t=1:T], vSTART_DAC[y,t] <= EP[:vCAP_DAC][y] / dfDac[y,:Cap_Size]
		    [y in DAC_COMMIT, t=1:T], vSHUT_DAC[y,t] <= EP[:vCAP_DAC][y] / dfDac[y,:Cap_Size]
	    end)


        #max
        
        @constraints(EP, begin
        # Minimum negative CO2 per DAC "y" at hour "t" > Min stable CO2 capture
            [y in DAC_COMMIT, t=1:T], EP[:vCO2_DAC][y,t] >= dfDac[y,:Min_DAC]*dfDac[y,:Cap_Size]*EP[:vCOMMIT_DAC][y,t]

        # Maximum negative CO2 per DAC "y"  "y" at hour "t" < Max capacity
            [y in DAC_COMMIT, t=1:T], EP[:vCO2_DAC][y,t] <= dfDac[y,:Cap_Size]*EP[:vCOMMIT_DAC][y,t]

         end)


        	# Commitment state constraint linking startup and shutdown decisions (Constraint #4)
	    p = hours_per_subperiod
        @constraints(EP, begin
            [y in DAC_ID, t in 1:T], vCOMMIT_DAC[y,t] == vCOMMIT_DAC[y, hoursbefore(p, t, 1)] + vSTART_DAC[y,t] - vSHUT_DAC[y,t]
        end)

        	### Minimum up and down times (Constraints #9-10)
        
	    Up_Time = zeros(Int, nrow(dfDac))
	    Up_Time[DAC_COMMIT] .= Int.(floor.(dfDac[DAC_COMMIT,:Up_Time]))
	    @constraint(EP, [y in DAC_COMMIT, t in 1:T],
	    	EP[:vCOMMIT_DAC][y,t] >= sum(EP[:vSTART_DAC][y, hoursbefore(p, t, 0:(Up_Time[y] - 1))])
	    )

	    Down_Time = zeros(Int, nrow(dfDac))
	    Down_Time[DAC_COMMIT] .= Int.(floor.(dfDac[DAC_COMMIT,:Down_Time]))
	    @constraint(EP, [y in DAC_COMMIT, t in 1:T],
	    	EP[:vCAP_DAC][y]/dfDac[y,:Cap_Size]-EP[:vCOMMIT_DAC][y,t] >= sum(EP[:vSHUT_DAC][y, hoursbefore(p, t, 0:(Down_Time[y] - 1))])
	    )
    end

    #max capacity constraint for DAC
    	# Constraint on maximum capacity (if applicable) [set input to -1 if no constraint on maximum capacity]
	# DEV NOTE: This constraint may be violated in some cases where Existing_Cap_MW is >= Max_Cap_MW and lead to infeasabilty
    @constraint(EP, cMaxCap_DAC[y in intersect(dfDac[dfDac.Max_Cap_DAC.>0,:R_ID], 1:G_DAC)], vCAP_DAC[y] <= dfDac[y,:Max_Cap_DAC])

    # get the CO2 balance 
    # the net negative CO2 for each DAC y at each hour t, CO2 emissions from heat consumption minus CO2 captured by DAC = net negative emissions
    @expression(EP, eCO2_DAC_net[y in DAC_ID, t = 1:T], eDAC_heat_CO2[y,t] - vCO2_DAC[y,t] )
    # the net negative CO2 from all the DAC
    @expression(EP, eCO2_DAC_net_ByZoneT[z = 1:Z, t = 1:T], 
        sum(eCO2_DAC_net[y, t] for y in dfDac[(dfDac[!, :Zone].==z), :DAC_ID]))
    # the net negative CO2 from all the DAC during the whole year
    @expression(EP, eCO2_DAC_net_ByZone[z = 1:Z], 
        sum(eCO2_DAC_net_ByZoneT[z, t] for t in 1:T))
    # sum of net CO2 across all the zone
    @expression(EP, eCO2_ToT_DAC_net, sum(eCO2_DAC_net_ByZone[z] for z in 1:Z))


    # separately account for the amount of CO2 that is captured.
    # actually the eCO2_net should be the total sequestration carbon. since eCO2_DAC_net should be a negative value, put minus sign in front of it..
    # costs associated with co2 transport & storage ($/(t CO2/h)) = captured co2 (t CO2/h) * Co2 transport and storage cost ($/t CO2)
    @expression(EP, eCCO2_TS_ByPlant[y in DAC_ID, t = 1:T], -eCO2_DAC_net[y, t]* dfDac[y, :CO2_Transport_Storage_Per_t])
    # the sequestrated CO2 from all the DAC 
    @expression(EP, eCCO2_TS_ByZoneT[z = 1:Z, t = 1:T], 
        sum(eCCO2_TS_ByPlant[y, t] for y in dfDac[(dfDac[!, :Zone].==z), :DAC_ID]))
    # the sequestrated CO2 from all the DAC during the whole year ($/t CO2)
    @expression(EP, eCCO2_TS_ByZone[z = 1:Z], 
        sum(eCCO2_TS_ByZoneT[z, t] for t in 1:T))
    # sum of CO2 sequestration costs.
    @expression(EP, eCTotalCO2TS, sum(eCCO2_TS_ByZone[z] for z in 1:Z))


    EP[:eObj] += eCTotalCO2TS

    
    return EP
end

