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

function dac(EP::Model, inputs::Dict)

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
    # the use of heat resources (e.g., natural gas) may result in additional CO2 emissions as DAC may not capture the CO2 from the heat resources
    # CO2 from the use of heat resource
    @expression(EP, eDAC_heat_CO2[y in DAC_ID, t = 1:T], eDAC_heat_consumption[y,t]* fuel_CO2[dfDAC[y,:DAC_heat_resource]] )
    # the electricity consumption for DAC, MWh/t CO2
    @expression(EP, eDAC_power[y in DAC_ID, t = 1:T], vCO2_DAC[y,t] * dfDac[y,:Electricity_MWh_per_CO2_metric_ton])
    # the power used for DAC is draw from the grid, so we have a power balance equation
	@expression(EP, ePowerBalanceDAC[t=1:T, z=1:Z],
		sum(eDAC_power[y,t] for y in intersect(dfDac[dfDAC[!,:Zone].==z,:][!,:DAC_ID])))
		
	EP[:ePowerBalance] = EP[:ePowerBalance] - ePowerBalanceDAC

    # finally, the CO2 captured by DAC during each hr should be less than the capacity
	@constraint(EP, [y in DAC_ID, t=1:T],vCO2_DAC[y,t] <= vCAP_DAC[y])

    #---------------------------------- add up cost ---------------------------------------
    # Fixed Cost
    # Combine CAPEX and FOM into annulized Fixed Cost
    # Fixed cost for a DAC y 
	@expression(EP, eCFixed_DAC[y in TS], dfDac[y,:Fixed_Cost_per_MW_th] * vCAP_DAC[y])
	# total fixed costs for all the DAC
	@expression(EP, eTotalCFixedDAC, sum(eCFixed_DAC[y] for y in DAC_ID))
	EP[:eObj] += eTotalCFixedDAC

    # Variable cost
    omega = inputs["omega"]
    # the total variable cost (heat cost + non-fuel vom cost) for DAC y at time t, $/t CO2  
    @expression(EP, eCDAC_Variable[y in DAC_ID, t = 1:T], omega * (eDAC_heat_consumption[y,t] * fuel_cost[dfDAC[y,:DAC_heat_resource]][t] +dfDac[y,:Var_OM_Cost_per_CO2]*vCO2_DAC[y,t]))
    # Cost associated with variable cost for all the DAC at time T
    @expression(EP, eCTotalVariableDACT[ t = 1:T], eCDAC_Variable[y,t] for y in DAC_ID )
    # Total vairable cost for all the DAC at all the time
    @expression(EP, eCTotalVariableDAC, eCTotalVariableDACT[t] for t in 1:T )

    EP[:eObj] += eCTotalVariableDAC
    
    # skip the unit commitment for DAC for now, will implement constraint formulation later...
    # get the CO2 balance 
    # the net negative CO2 for each DAC y at each hour t, CO2 emissions from heat consumption minus CO2 captured by DAC = net negative emissions
    @expression(EP, eCO2_DAC_net[y in DAC_ID, t = 1:T], eDAC_heat_CO2[y,t] - vCO2_DAC[y,t] )

    expression(EP, eCO2_DAC_net_ByZone[z = 1:Z, t = 1:T], 
        sum(eCO2_DAC_net[y, t] for y in dfDac[(dfDac[!, :Zone].==z), :DAC_ID]))
    # the net negative CO2 from all the DAC 
    
    return EP
end