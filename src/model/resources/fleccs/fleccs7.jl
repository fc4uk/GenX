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
	FLECCS7(EP::Model, inputs::Dict, UCommit::Int, Reserves::Int)

The FLECCS7 module creates decision variables, expressions, and constraints related to FLECCS Process with DAC (MIT). In this module, we will write up all the constraints formulations associated with the power plant.
"""

function fleccs7(EP::Model, inputs::Dict)

	println("FLECCS7, FLECCS Process with DAC (MIT)")

	T = inputs["T"]     # Number of time steps (hours)
    Z = inputs["Z"]     # Number of zones
    G_F = inputs["G_F"] # Number of FLECCS generator
	FLECCS_ALL = inputs["FLECCS_ALL"] # set of FLECCS generator
	dfGen_ccs = inputs["dfGen_ccs"] # FLECCS general data
	#dfGen_ccs = inputs["dfGen_ccs"] # FLECCS specific parameters
	# get number of flexible subcompoents
	N_F = inputs["N_F"]
	n = length(N_F)
 


	NEW_CAP_ccs = inputs["NEW_CAP_FLECCS"] #allow for new capcity build
	RET_CAP_ccs = inputs["RET_CAP_FLECCS"] #allow for retirement

	START_SUBPERIODS = inputs["START_SUBPERIODS"] #start
    INTERIOR_SUBPERIODS = inputs["INTERIOR_SUBPERIODS"] #interiors

    hours_per_subperiod = inputs["hours_per_subperiod"]

	fuel_type = collect(skipmissing(dfGen_ccs[!,:Fuel]))

	fuel_CO2 = inputs["fuel_CO2"]
	fuel_costs = inputs["fuel_costs"]



	STARTS = 1:inputs["H"]:T
    # Then we record all time periods that do not begin a sub period
    # (these will be subject to normal time couping constraints, looking back one period)
    INTERIORS = setdiff(1:T,STARTS)

	# capacity decision variables


	# varibales related to power generation/consumption
    @variables(EP, begin
        # Continuous decision variables
        vP_NGCC[y in FLECCS_ALL, 1:T] >= 0  # Net power produced by the NGCC
		vFuel_use_NGCC[y in FLECCS_ALL, 1:T] >= 0 # Molar flow rate of fuel into the NGCC block
        vCaCO3_use_calciner[y in FLECCS_ALL, 1:T] >= 0  # Molar flow rate of fresh CaCO3 into the calciner block
		vCO2_atmosphere[y in FLECCS_ALL, 1:T] >=0  # mole of CO2 captured from atmosphere 
		nCaO_DAC[y in FLECCS_ALL, 1:T] >= 0  # mol of CaO stored in the DAC unit

    end)

	# variables related to CO2 and solvent
	# @variable(EP, NGCC_on[y in FLECCS_ALL, 1:T], Bin)

	# the order of those variables must follow the order of subcomponents in the "FLECCS_data.csv"
	# 1. NGCC
	# 2. CAL
	# 3. DAC


	# get the ID of each subcompoents 
	# NGCC
	NGCC_id = inputs["NGCC_id"]
	# CAL
	CAL_id = inputs["CAL_id"]
	# DAC
	DAC_id = inputs["DAC_id"]


	# Specific constraints for FLECCS system
	# Piecewise heat rate UC Equation 16 and 17
	@constraint(EP, [y in FLECCS_ALL, t = 1:T], 
	    vFuel_use_NGCC[y,t] >= vP_NGCC[y,t]*dfGen_ccs[!,:c1][y] + EP[:vCOMMIT_FLECCS][y, NGCC_id, t]*dfGen_ccs[!,:d1][y])
	@constraint(EP, [y in FLECCS_ALL, t = 1:T],
		vFuel_use_NGCC[y,t] >= vP_NGCC[y,t]*dfGen_ccs[!,:c2][y] + EP[:vCOMMIT_FLECCS][y, NGCC_id, t]*dfGen_ccs[!,:d2][y])
	#@expression(EP, vFuel_use_NGCC[y in FLECCS_ALL, t = 1:T], vP_NGCC[y,t] * 6.8)


    # (1) Molar flow rate of CO2 in the flue gas, eCO2use (kmol/hr) = 1:35eFuel_use_NGCCt (MMBtu/hr), p1 =  1.35
    @expression(EP, eCO2_flue[y in FLECCS_ALL,t=1:T], dfGen_ccs[!,:P1_1][1+n*(y-1)] * vFuel_use_NGCC[y,t])
	# (5) Molar flow rate of fuel into the calciner block, need to define eFuel_use_calciner first before using it in a formulation.
	@expression(EP, eFuel_use_calciner[y in FLECCS_ALL,t=1:T], dfGen_ccs[!,:P5_1][1+n*(y-1)] * vP_NGCC[y,t] + dfGen_ccs[!,:P5_2][1+n*(y-1)] * vCaCO3_use_calciner[y,t] +  EP[:vCOMMIT_FLECCS][y, NGCC_id, t]* dfGen_ccs[!,:P5_3][1+n*(y-1)])
	# (2) Molar flow rate of CO2 from the calciner
    @expression(EP, eCO2_calciner[y in FLECCS_ALL,t=1:T], dfGen_ccs[!,:P2_1][1+n*(y-1)] * vP_NGCC[y,t] + dfGen_ccs[!,:P2_2][1+n*(y-1)] *vCaCO3_use_calciner[y,t] + EP[:vCOMMIT_FLECCS][y, NGCC_id, t]*dfGen_ccs[!,:P2_3][1+n*(y-1)] )
	
    # (3) CO2 vented to the atmosphere
	@expression(EP, eCO2_vent[y in FLECCS_ALL,t=1:T], dfGen_ccs[!,:P3_1][1+n*(y-1)] * vP_NGCC[y,t] + dfGen_ccs[!,:P3_2][1+n*(y-1)] * vCaCO3_use_calciner[y,t] +  EP[:vCOMMIT_FLECCS][y, NGCC_id, t]*dfGen_ccs[!,:P3_3][1+n*(y-1)])
    
	# (4) Molar flow rate of CO2 exiting the separation block
	@expression(EP, CO2_liquified[y in FLECCS_ALL,t=1:T], eCO2_calciner[y,t] - eCO2_vent[y,t])
	
	# (6) Net power required by block with calcium loop
    @expression(EP, vP_calciner[y in FLECCS_ALL,t=1:T], dfGen_ccs[!,:P6_1][1+n*(y-1)] * vP_NGCC[y,t] + dfGen_ccs[!,:P6_2][1+n*(y-1)] * vCaCO3_use_calciner[y,t] + EP[:vCOMMIT_FLECCS][y, NGCC_id, t]* dfGen_ccs[!,:P6_3][1+n*(y-1)])
    
	# (7) Molar flow rate of CaO from the block with calcium loop
    @expression(EP, CaO_out_calciner[y in FLECCS_ALL,t=1:T], vCaCO3_use_calciner[y,t])

    # (8) Net power required by DAC
    @expression(EP, vP_use_DAC[y in FLECCS_ALL,t=1:T], dfGen_ccs[!,:P8_1][1+n*(y-1)] * vCaCO3_use_calciner[y,t])
   
	# (9) Net power produced by the overall plant (vP_out = eCCS_net) ///is it meant to be @expression(EP, eCCS_net[y in FLECCS_ALL,t=1:T], vP_NGCC[y,t] + vP_calciner[y,t] + vP_use_DAC[y,t])
    @expression(EP, eCCS_net[y in FLECCS_ALL,t=1:T], vP_NGCC[y,t] - vP_calciner[y,t] - vP_use_DAC[y,t])

	# (12) Used CaO is protional to the CO2 captured by atmosphere, need to reorder those equations.
	@expression(EP, CaO_use[y in FLECCS_ALL,t=1:T], vCO2_atmosphere[y,t]/dfGen_ccs[!,:alpha][1+n*(y-1)])
	# (11) Moles of fully converted CaO particles exported from DAC
	#@constraint(EP, [y in FLECCS_ALL,t=1:24], CaO_use[y,t] == CaO_out_calciner[y,t+8736])
	@constraint(EP, [y in FLECCS_ALL,t=1:24], CaO_use[y,t] == 0)
	@constraint(EP, [y in FLECCS_ALL,t=25:T], CaO_use[y,t] ==CaO_out_calciner[y,t-24])

	# (10) Moles of CaO in the DAC block
    # dynamic of DAC system, normal [tonne solvent/sorbent]
    @constraint(EP, [y in FLECCS_ALL, t in INTERIOR_SUBPERIODS],nCaO_DAC[y, t] == nCaO_DAC[y, t-1] +  (CaO_out_calciner[y,t] - CaO_use[y,t]))
    # dynamic of DAC system, wrapping [tonne solvent/sorbent]
	@constraint(EP, [y in FLECCS_ALL, t in START_SUBPERIODS], nCaO_DAC[y, t] == 0 )

    #@constraint(EP, [y in FLECCS_ALL, t in START_SUBPERIODS], nCaO_DAC[y, t] == nCaO_DAC[y,t+hours_per_subperiod-1] +  (CaO_out_calciner[y,t] - CaO_use[y,t]))


	# Power Balance /// I'm using eCCS_net to replace vP_out 
	@expression(EP, ePowerBalanceFLECCS[t=1:T, z=1:Z], sum(eCCS_net[y,t] for y in unique(dfGen_ccs[(dfGen_ccs[!,:Zone].==z),:R_ID])))
	EP[:ePowerBalance] += ePowerBalanceFLECCS

    #@constraint(EP,[y in FLECCS_ALL,t = 1:T],  vCaCO3_use_calciner[y,t] >= 3000)


	@constraint(EP,[y in FLECCS_ALL, t=1:T], vCaCO3_use_calciner[y,t] >= EP[:eTotalCapFLECCS][y, NGCC_id] * 3000/688)
	@constraint(EP,[y in FLECCS_ALL, t=1:T], vCaCO3_use_calciner[y,t] <= EP[:eTotalCapFLECCS][y, NGCC_id] * 6500/688)
	@constraint(EP,[y in FLECCS_ALL, t=1:T], nCaO_DAC[y,t] <= EP[:eTotalCapFLECCS][y, NGCC_id] * 156000/688)


	@constraints(EP, begin
	    [y in FLECCS_ALL, i in NGCC_id, t = 1:T], EP[:vFLECCS_output][y,i,t] == vP_NGCC[y,t]
		[y in FLECCS_ALL, i in CAL_id,t = 1:T],EP[:vFLECCS_output][y,i,t] == vCaCO3_use_calciner[y,t]
		[y in FLECCS_ALL, i in DAC_id,t = 1:T],EP[:vFLECCS_output][y,i,t] == nCaO_DAC[y,t]
	end)

	@constraint(EP, [y in FLECCS_ALL, t = 1:T], vCO2_atmosphere[y,t] <= EP[:eTotalCapFLECCS][y, DAC_id])



    ########### amount of CO2 sequestration 
	@expression(EP, eCO2_sequestration[y in FLECCS_ALL,t=1:T], CO2_liquified[y,t])

    ######################
	# Cost
	# Fuel: Cost of natural gas. inputs["omega"] should be added to all the variable cost related formulations
	@expression(EP, eCVar_fuel[y in FLECCS_ALL, t = 1:T], fuel_costs[fuel_type[1]][t] * (eFuel_use_calciner[y,t] + vFuel_use_NGCC[y,t]))
		# Sum to annual level
	@expression(EP, eCFuelFLECCS[y in FLECCS_ALL],inputs["omega"][t]*sum(eCVar_fuel[y,t] for t in 1:T))
    # Sum to zonal-annual level
	@expression(EP, eZonalCFuelFLECCS[z = 1:Z], (sum(EP[:eCFuelFLECCS][y] for y in unique(dfGen_ccs[(dfGen_ccs[!,:Zone].==z),:R_ID]))))
	# Sum to system level
    @expression(EP, eTotalCFuelFLECCS, sum(eZonalCFuelFLECCS[z] for z in 1:Z))


	
	
	# Cost of feed CaCO3
	@expression(EP, eCVar_CaCO3[y in FLECCS_ALL, t = 1:T],  inputs["omega"][t]*dfGen_ccs[!,:cost_limestone][1+n*(y-1)] * vCaCO3_use_calciner[y,t])
	# Cost of limestone transport
	@expression(EP, eCVar_CaCO_transport[y in FLECCS_ALL, t = 1:T], inputs["omega"][t]*dfGen_ccs[!,:cost_limestone_transport][1+n*(y-1)] * CaO_use[y,t])
	# Carbon tax/credit, this is implemented in fleccs.jl
	#@expression(EP, eCVar_CO2[y in FLECCS_ALL, t = 1:T], dfGen_ccs[!,:cost_carbon][y] * (eCO2_vent[y,t]-CO2_atmosphere[y,t]))
	# CO2 sequestration cost applied to sequestrated CO2
	@expression(EP, eCVar_CO2_sequestration[y in FLECCS_ALL, t = 1:T], inputs["omega"][t]*dfGen_ccs[!,:cost_seq][1+n*(y-1)] * CO2_liquified[y,t])

	#adding up variable costs
	@expression(EP, eCVOMFLECCS[y in FLECCS_ALL], sum(eCVar_CaCO3[y,t] +  eCVar_CaCO_transport[y,t] + eCVar_CO2_sequestration[y,t] for t in 1:T))
    # Sum to zonal-annual level
	@expression(EP, eZonalCVOMFLECCS[z = 1:Z], (sum(EP[:eCVOMFLECCS][y] for y in unique(dfGen_ccs[(dfGen_ccs[!,:Zone].==z),:R_ID]))))
	# Sum to system level
    @expression(EP, eTotalCVOMFLECCS, sum(eZonalCVOMFLECCS[z] for z in 1:Z))

	# total variable cost 
	@expression(EP, eCVarFLECCS[y in FLECCS_ALL],eCVOMFLECCS[y] + eCFuelFLECCS[y])
	@expression(EP, eZonalCVarFLECCS[z = 1:Z],  sum(EP[:eCVarFLECCS][y] for y in unique(dfGen_ccs[(dfGen_ccs[!,:Zone].==z),:R_ID])))
    @expression(EP, eTotalCVarFLECCS, sum(EP[:eZonalCVarFLECCS][z] for z in 1:Z))


	# Add total variable discharging cost contribution to the objective function
	EP[:eObj] += eTotalCVarFLECCS

	return EP
end
