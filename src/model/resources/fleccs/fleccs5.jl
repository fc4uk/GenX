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
	FLECCS5(EP::Model, inputs::Dict, UCommit::Int, Reserves::Int)
The FLECCS5 module creates decision variables, expressions, and constraints related to NGCC-CCS coupled with H2 generation and storage systems. 
"""

function fleccs5(EP::Model, inputs::Dict)

	println("FLECCS5, NGCC coupled with H2 generation and storage Resources Module")

	T = inputs["T"]     # Number of time steps (hours)
    Z = inputs["Z"]     # Number of zones
    G_F = inputs["G_F"] # Number of FLECCS generator
	FLECCS_ALL = inputs["FLECCS_ALL"] # set of FLECCS generator
	dfGen_ccs = inputs["dfGen_ccs"] # FLECCS general data
	# get number of flexible subcompoents
	N_F = inputs["N_F"]
	n = length(N_F)
 

	START_SUBPERIODS = inputs["START_SUBPERIODS"] #start
    INTERIOR_SUBPERIODS = inputs["INTERIOR_SUBPERIODS"] #interiors
    hours_per_subperiod = inputs["hours_per_subperiod"]

	fuel_type = collect(skipmissing(dfGen_ccs[!,:Fuel]))
	fuel_CO2 = inputs["fuel_CO2"]
	fuel_costs = inputs["fuel_costs"]




	# variales related to power generation/consumption
    @variables(EP, begin
        # Continuous decision variables
        vP_gt[y in FLECCS_ALL, 1:T]  >= 0 # generation from combustion TURBINE (gas TURBINE)
    end)

	# variales related to CO2, H2, and storage
	@variables(EP, begin
        vCAPTURE[y in FLECCS_ALL,1:T] >= 0 # captured CO2, tonne/h
        vSTORE_H2[y in FLECCS_ALL,1:T] >= 0 # h2 storage, MMBTU H2
		vH2_in[y in FLECCS_ALL,1:T] >= 0 # h2 in, MMBTU H2/h
		vH2_out[y in FLECCS_ALL,1:T] >= 0 # h2 out, MMBTU H2/h
		vH2_sale[y in FLECCS_ALL,1:T] >= 0 # h2 sold for credit,  MMBTU H2/h
	end)

	# get the ID of each subcompoents 
	# gas turbine 
	NGCT_id = inputs["NGCT_id"]
	# steam turbine
	NGST_id = inputs["NGST_id"]
	# absorber 
	PCC_id = inputs["PCC_id"]
	# compressor
	Comp_id = inputs["Comp_id"] 
	# electrolyzer
	Electrolyzer_id = inputs["Electrolyzer_id"] 
	# h2 compressor
	H2_Compressor_id = inputs["H2_Compressor_id"] 
	# h2 tank 
	H2_Tank_id = inputs["H2_Tank_id"]
	#BOP 
	BOP_id = inputs["BOP_id"]

	# Specific constraints for FLECCS system
    # eq1a: Thermal Energy input of combustion TURBINE (or oxyfuel power cycle) at hour "t" [MMBTU]
	# depending on how much H2 (vH2_out) is feed into the combustion turbine, the NG input is determined by the energy difference.
    @expression(EP, eFuel[y in FLECCS_ALL,t=1:T], dfGen_ccs[!,:pHeatRate_gt][1+n*(y-1)] * vP_gt[y,t])
    # eq1b: The thermal energy input (efuel_(y,t)) can be from H2 and NG, depending on how much H2 (vH2_out) is feed into the combustion turbine, the NG input is determined by the energy difference.
	@expression(EP, eFuel_NG[y in FLECCS_ALL,t=1:T], eFuel[y,t] - vH2_out[y,t])

    # eqn 2a - 2c
	# stema generated by high, mid, low, pressure turbine
	@expression(EP, eSteam_high[y in FLECCS_ALL,t=1:T], dfGen_ccs[!,:pSteamRate_high][1+n*(y-1)] * eFuel[y,t])
	# mid pressure steam
	@expression(EP, eSteam_mid[y in FLECCS_ALL,t=1:T], dfGen_ccs[!,:pSteamRate_mid][1+n*(y-1)] * eSteam_high[y,t] )
	# low pressure steam
	@expression(EP, eSteam_low[y in FLECCS_ALL,t=1:T], dfGen_ccs[!,:pSteamRate_low][1+n*(y-1)] * eSteam_mid[y,t])
	
	# eqn 3, CO2 generated by combustion TURBINE (or oxyfuel power cycle) at hour "t" [tonne]
    @expression(EP, eCO2_flue[y in FLECCS_ALL,t=1:T], inputs["CO2_per_MMBTU_FLECCS"][y,NGCT_id] * eFuel_NG[y,t])

	# eq4, vCAPTURE must less than eCO2_flue
	@constraint(EP, [y in FLECCS_ALL,t=1:T], vCAPTURE[y,t] <= dfGen_ccs[!,:pCO2CapRate][1+n*(y-1)]*eCO2_flue[y,t] )
	
	# eq5, CO2 vented at time "t" [tonne]
    @expression(EP, eCO2_vent[y in FLECCS_ALL,t=1:T], eCO2_flue[y,t] - vCAPTURE[y,t])

    # eq6, steam used by post-combustion carbon capture (PCC) unit [MMBTU], since steam generated by auxiliary boiler could be used to regenerate solvent, vSTEAM_in is incorporated into this equation.
    @expression(EP, eSteam_use_pcc[y in FLECCS_ALL,t=1:T], dfGen_ccs[!,:pSteamUseRate][1+n*(y-1)] * vCAPTURE[y,t])
    
	# eq7, power used by post-combustion carbon capture (PCC) unit [MWh]
    @expression(EP, ePower_use_pcc[y in FLECCS_ALL,t=1:T], dfGen_ccs[!,:pPowerUseRate][1+n*(y-1)]  * vCAPTURE[y,t])
   
    # eq8, power used by compressor unit [MWh]
    @expression(EP, ePower_use_comp[y in FLECCS_ALL,t=1:T], dfGen_ccs[!,:pCO2CompressRate][1+n*(y-1)] * vCAPTURE[y,t])
	
	# eq 9, power used by auxiliary [MWh]
	@expression(EP, ePower_use_other[y in FLECCS_ALL,t=1:T], dfGen_ccs[!,:pPowerUseRate_Other][1+n*(y-1)] * vP_gt[y,t])

	#eqn 10a - 10 d, Power generated by steam turbine [MWh]
	@expression(EP, ePower_st_high[y in FLECCS_ALL,t=1:T], eSteam_high[y,t]/dfGen_ccs[!,:pHeatRate_st_high][1+n*(y-1)] )

	@expression(EP, ePower_st_mid[y in FLECCS_ALL,t=1:T], eSteam_mid[y,t]/dfGen_ccs[!,:pHeatRate_st_mid][1+n*(y-1)])
	
	@expression(EP, ePower_st_low[y in FLECCS_ALL,t=1:T], (eSteam_low[y,t] - eSteam_use_pcc[y,t])/dfGen_ccs[!,:pHeatRate_st_low][1+n*(y-1)])
	
	@constraint(EP, [y in FLECCS_ALL,t=1:T], ePower_st_low[y,t] >= 0)

	@expression(EP, ePower_st[y in FLECCS_ALL,t=1:T], ePower_st_high[y,t] +ePower_st_mid[y,t] +ePower_st_low[y,t] )

    # eqn 11. power used by electrolyzer for H2 generation
	@expression(EP, ePower_use_h2[y in FLECCS_ALL,t=1:T], dfGen_ccs[!,:pPowerUseRate_H2][1+n*(y-1)]  * vH2_in[y,t])
	
	@expression(EP, ePower_use_h2_comp[y in FLECCS_ALL,t=1:T], dfGen_ccs[!,:pPowerUseRate_H2_comp][1+n*(y-1)]  * vH2_in[y,t])
  
	#eqn 12, maximum H2 ratio
	@constraint(EP, [y in FLECCS_ALL,t=1:T], vH2_out[y,t] <= eFuel[y,t] * dfGen_ccs[!,:pMaxH2Ratio][1+n*(y-1)])

	# eqn 13 - 14, H2 storage mass balance
	# dynamic of H2 storage system, normal
	@constraint(EP, cStore_H2[y in FLECCS_ALL, t in INTERIOR_SUBPERIODS],vSTORE_H2[y, t] == vSTORE_H2[y, t-1] + vH2_in[y,t] - vH2_out[y,t] - vH2_sale[y,t])
	# dynamic of H2, wrapping 
	@constraint(EP, cStore_H2wrap[y in FLECCS_ALL, t in START_SUBPERIODS],vSTORE_H2[y, t] == vSTORE_H2[y,t+hours_per_subperiod-1] +  vH2_in[y,t] - vH2_out[y,t]- vH2_sale[y,t])


	# eqn15, NGCC-CCS net power output 
	@expression(EP, eCCS_net[y in FLECCS_ALL,t=1:T], vP_gt[y,t] + ePower_st[y,t] - ePower_use_comp[y,t]- ePower_use_pcc[y,t] - ePower_use_other[y,t] - ePower_use_h2[y,t] -  ePower_use_h2_comp[y,t])
    
	# add to power balance
	@expression(EP, ePowerBalanceFLECCS[t=1:T, z=1:Z], sum(eCCS_net[y,t] for y in unique(dfGen_ccs[(dfGen_ccs[!,:Zone].==z),:R_ID])))

	## Power Balance##
	EP[:ePowerBalance] += ePowerBalanceFLECCS

	# create a container for FLECCS output.
	@constraints(EP, begin
	    [y in FLECCS_ALL, i in NGCT_id, t = 1:T],EP[:vFLECCS_output][y,i,t] == vP_gt[y,t]
		[y in FLECCS_ALL, i in NGST_id,t = 1:T],EP[:vFLECCS_output][y,i,t] == ePower_st[y,t]	
		[y in FLECCS_ALL, i in PCC_id,t = 1:T],EP[:vFLECCS_output][y,i,t] == vCAPTURE[y,t]
		[y in FLECCS_ALL, i in Comp_id, t =1:T],EP[:vFLECCS_output][y,i,t] == ePower_use_comp[y,t]	
		[y in FLECCS_ALL, i in Electrolyzer_id, t = 1:T],EP[:vFLECCS_output][y,i,t] == ePower_use_h2[y,t]
		[y in FLECCS_ALL, i in H2_Compressor_id, t =1:T],EP[:vFLECCS_output][y,i,t] == ePower_use_h2_comp[y,t]
		[y in FLECCS_ALL, i in H2_Tank_id, t =1:T],EP[:vFLECCS_output][y,i,t] == vSTORE_H2[y,t]	
		#[y in FLECCS_ALL, i in BOP_id, t =1:T],EP[:vFLECCS_output][y,i,t] == eCCS_net[y,t]			
	end)

	@constraint(EP, [y in FLECCS_ALL], EP[:eTotalCapFLECCS][y, BOP_id] == EP[:eTotalCapFLECCS][y, NGCT_id]+ EP[:eTotalCapFLECCS][y,NGST_id])


    ########### amount of CO2 sequestration 
	@expression(EP, eCO2_sequestration[y in FLECCS_ALL,t=1:T], vCAPTURE[y,t])




	###########variable cost
	#fuel, eFuel_NG is used
	@expression(EP, eCVar_fuel[y in FLECCS_ALL, t = 1:T],(fuel_costs[fuel_type[1]][t]*eFuel_NG[y,t]))
	# Sum to annual level
	@expression(EP, eCFuelFLECCS[y in FLECCS_ALL],inputs["omega"][t]*sum(eCVar_fuel[y,t] for t in 1:T))
    # Sum to zonal-annual level
	@expression(EP, eZonalCFuelFLECCS[z = 1:Z], (sum(EP[:eCFuelFLECCS][y] for y in unique(dfGen_ccs[(dfGen_ccs[!,:Zone].==z),:R_ID]))))
	# Sum to system level
    @expression(EP, eTotalCFuelFLECCS, sum(eZonalCFuelFLECCS[z] for z in 1:Z))

	#VOM
	# CO2 sequestration cost applied to sequestrated CO2
	@expression(EP, eCVar_CO2_sequestration[y in FLECCS_ALL, t = 1:T],(inputs["omega"][t]*eCO2_sequestration[y,t]*dfGen_ccs[!,:pCO2_sequestration][1+n*(y-1)]))
	# variable O&M for all the teams: combustion turbine (or oxfuel power cycle)
	@expression(EP,eCVar_gt[y in FLECCS_ALL, t = 1:T], inputs["omega"][t]*(dfGen_ccs[(dfGen_ccs[!,:FLECCS_NO].==NGCT_id) .& (dfGen_ccs[!,:R_ID].==y),:Var_OM_Cost_per_Unit][1])*vP_gt[y,t])
	# variable O&M for NGCC-based teams: VOM of steam turbine and co2 compressor
	# variable O&M for steam turbine
	@expression(EP,eCVar_st[y in FLECCS_ALL, t = 1:T], inputs["omega"][t]*(dfGen_ccs[(dfGen_ccs[!,:FLECCS_NO].==NGST_id) .& (dfGen_ccs[!,:R_ID].==y),:Var_OM_Cost_per_Unit][1])*ePower_st[y,t])
	 # variable O&M for compressor
	@expression(EP,eCVar_comp[y in FLECCS_ALL, t = 1:T], inputs["omega"][t]*(dfGen_ccs[(dfGen_ccs[!,:FLECCS_NO].== Comp_id) .& (dfGen_ccs[!,:R_ID].==y),:Var_OM_Cost_per_Unit][1])*(eCO2_sequestration[y,t])) 
	# assume H2 = $14/mmbtu, eCH2 = credit from H2 sales
	@expression(EP,eCH2[y in FLECCS_ALL, t = 1:T], inputs["omega"][t]*(vH2_sale[y,t]*dfGen_ccs[!,:H2_per_MMBTU][1+n*(y-1)]))
	

	#adding up variable cost
	@expression(EP, eCVOMFLECCS[y in FLECCS_ALL], sum(eCVar_CO2_sequestration[y,t] + eCVar_gt[y,t] + eCVar_st[y,t] + eCVar_comp[y,t] - eCH2[y,t]  for t in 1:T))
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
