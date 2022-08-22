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
"""

function fleccs_commit(EP::Model, inputs::Dict, UCommit::Int)

	println("FLECCS (Unit Commitment) Resources Module")

	dfGen_ccs = inputs["dfGen_ccs"]

	T = inputs["T"]     # Number of time steps (hours)
	Z = inputs["Z"]     # Number of zones
	G = inputs["G"]     # Number of resources


	FLECCS_ALL = inputs["FLECCS_ALL"] # set of FLECCS generator
	N_F = inputs["N_F"] 	# get number of flexible subcompoents
	n_F = length(N_F)
	COMMIT_ccs = inputs["COMMIT_CCS"] # CCS compoents subjected to UC
 

	hours_per_subperiod = inputs["hours_per_subperiod"] #total number of hours per subperiod
	START_SUBPERIODS = inputs["START_SUBPERIODS"]
	INTERIOR_SUBPERIODS = inputs["INTERIOR_SUBPERIODS"]



	### Variables ###

	## Decision variables for unit commitment
	# gas turbine and steam turbine are grouped into vCOMMIT_NGCC, for Allam cycle, we also use vCOMMIT_NGCC to represent the commitment status
	# commitment state variable
	@variable(EP, vCOMMIT_FLECCS[y in FLECCS_ALL, i in COMMIT_ccs, t=1:T] >= 0)
	# startup event variable
	@variable(EP, vSTART_FLECCS[y in FLECCS_ALL, i in COMMIT_ccs, t=1:T] >= 0)
	# shutdown event variable
	@variable(EP, vSHUT_FLECCS[y in FLECCS_ALL, i in COMMIT_ccs, t=1:T] >= 0)

	### Expressions ###
	## Objective Function Expressions ##
	# Startup costs of "generation" for resource "y" during hour "t"
	@expression(EP, eCStart_FLECCS[y in FLECCS_ALL, i in COMMIT_ccs, t=1:T],(inputs["omega"][t]*inputs["C_Start_FLECCS"][y,i,t]*vSTART_FLECCS[y,i,t]))
	# Julia is fastest when summing over one row one column at a time
	
	@expression(EP, eZonalCStartFLECCS[z = 1:Z], sum(eCStart_FLECCS[y,i,t] for y in  unique(dfGen_ccs[(dfGen_ccs[!, :Zone].==z), :R_ID]), i in COMMIT_ccs, t in 1:T))
	
    @expression(EP, eTotalCStartFLECCS, sum(eZonalCStartFLECCS[z] for z in 1:Z))


	EP[:eObj] +=  eTotalCStartFLECCS

	### Constratints ###
	## Declaration of integer/binary variables
	if UCommit == 1 # Integer UC constraints
		for y in FLECCS_ALL
		    for i in COMMIT_ccs
			    set_integer.(vCOMMIT_FLECCS[y,i,:])
	    		set_integer.(vSTART_FLECCS[y,i,:])
		    	set_integer.(vSHUT_FLECCS[y,i,:])
			    if y in inputs["RET_CAP_FLECCS"]
    				set_integer(EP[:vRETCAP_FLECCS][y,i])
	    		end
		    	if y in inputs["NEW_CAP_FLECCS"]
			     	set_integer(EP[:vRETCAP_FLECCS][y,i])
				end
			end
		end
	end 




	### Expressions ###

	### Constraints ###

	### Capacitated limits on unit commitment decision variables (Constraints #1-3)
	@constraints(EP, begin
		[y in FLECCS_ALL, i in COMMIT_ccs, t=1:T], vCOMMIT_FLECCS[y,i,t] <= EP[:eTotalCapFLECCS][y,i]/dfGen_ccs[(dfGen_ccs[!,:R_ID].==y),:Cap_Size][i]
		[y in FLECCS_ALL, i in COMMIT_ccs, t=1:T], vSTART_FLECCS[y,i,t] <= EP[:eTotalCapFLECCS][y,i]/dfGen_ccs[(dfGen_ccs[!,:R_ID].==y),:Cap_Size][i]
		[y in FLECCS_ALL, i in COMMIT_ccs, t=1:T], vSHUT_FLECCS[y,i,t] <= EP[:eTotalCapFLECCS][y,i]/dfGen_ccs[(dfGen_ccs[!,:R_ID].==y),:Cap_Size][i]
	end)

	# Commitment state constraint linking startup and shutdown decisions (Constraint #4)
	@constraints(EP, begin
		# For Start Hours, links first time step with last time step in subperiod
		[y in FLECCS_ALL, i in COMMIT_ccs, t in START_SUBPERIODS], vCOMMIT_FLECCS[y,i,t] == vCOMMIT_FLECCS[y,i,(t+hours_per_subperiod-1)] + vSTART_FLECCS[y,i,t] - vSHUT_FLECCS[y,i,t]
		# For all other hours, links commitment state in hour t with commitment state in prior hour + sum of start up and shut down in current hour
		[y in FLECCS_ALL, i in COMMIT_ccs, t in INTERIOR_SUBPERIODS], vCOMMIT_FLECCS[y,i,t] == vCOMMIT_FLECCS[y,i,t-1] + vSTART_FLECCS[y,i,t] - vSHUT_FLECCS[y,i,t]
	end)

	### Maximum ramp up and down between consecutive hours (Constraints #5-6)
    ## 
	## For Start Hours
	# Links last time step with first time step, ensuring position in hour 1 is within eligible ramp of final hour position
		# rampup constraints
	@constraint(EP,[y in FLECCS_ALL, i in COMMIT_ccs, t in START_SUBPERIODS],
		EP[:vFLECCS_output][y,i,t]-EP[:vFLECCS_output][y,i,(t+hours_per_subperiod-1)] <= dfGen_ccs[(dfGen_ccs[!,:R_ID].==y),:Ramp_Up_Percentage][i]*dfGen_ccs[(dfGen_ccs[!,:R_ID].==y),:Cap_Size][i]*(vCOMMIT_FLECCS[y,i,t]-vSTART_FLECCS[y,i,t])
			+ min(1,max(dfGen_ccs[(dfGen_ccs[!,:R_ID].==y),:Min_Power][i],dfGen_ccs[!,:Ramp_Up_Percentage][y]))*dfGen_ccs[(dfGen_ccs[!,:R_ID].==y),:Cap_Size][i]*vSTART_FLECCS[y,i,t]
			-dfGen_ccs[(dfGen_ccs[!,:R_ID].==y),:Min_Power][i]*dfGen_ccs[(dfGen_ccs[!,:R_ID].==y),:Cap_Size][i]*vSHUT_FLECCS[y,i,t])

		# rampdown constraints
	@constraint(EP,[y in FLECCS_ALL, i in COMMIT_ccs, t in START_SUBPERIODS],
		EP[:vFLECCS_output][y,i,(t+hours_per_subperiod-1)]-EP[:vFLECCS_output][y,i,t] <=dfGen_ccs[(dfGen_ccs[!,:R_ID].==y),:Ramp_Dn_Percentage][i]*dfGen_ccs[(dfGen_ccs[!,:R_ID].==y),:Cap_Size][i]*(vCOMMIT_FLECCS[y,i,t]-vSTART_FLECCS[y,i,t])
			-dfGen_ccs[(dfGen_ccs[!,:R_ID].==y),:Min_Power][i]*dfGen_ccs[(dfGen_ccs[!,:R_ID].==y),:Cap_Size][i]*vSTART_FLECCS[y,i,t]
			+ min(1,max(dfGen_ccs[(dfGen_ccs[!,:R_ID].==y),:Min_Power][i],dfGen_ccs[!,:Ramp_Dn_Percentage][y]))*dfGen_ccs[(dfGen_ccs[!,:R_ID].==y),:Cap_Size][i]*vSHUT_FLECCS[y,i,t])

	## For Interior Hours
		# rampup constraints
	@constraint(EP,[y in FLECCS_ALL, i in COMMIT_ccs, t in INTERIOR_SUBPERIODS],
		EP[:vFLECCS_output][y,i,t]-EP[:vFLECCS_output][y,i,t-1] <=dfGen_ccs[(dfGen_ccs[!,:R_ID].==y),:Ramp_Up_Percentage][i]*dfGen_ccs[(dfGen_ccs[!,:R_ID].==y),:Cap_Size][i]*(vCOMMIT_FLECCS[y,i,t]-vSTART_FLECCS[y,i,t])
			+ min(1,max(dfGen_ccs[(dfGen_ccs[!,:R_ID].==y),:Min_Power][i],dfGen_ccs[!,:Ramp_Up_Percentage][y]))*dfGen_ccs[(dfGen_ccs[!,:R_ID].==y),:Cap_Size][i]*vSTART_FLECCS[y,i,t]
			-dfGen_ccs[(dfGen_ccs[!,:R_ID].==y),:Min_Power][i]*dfGen_ccs[(dfGen_ccs[!,:R_ID].==y),:Cap_Size][i]*vSHUT_FLECCS[y,i,t])

		# rampdown constraints
	@constraint(EP,[y in FLECCS_ALL, i in COMMIT_ccs, t in INTERIOR_SUBPERIODS],
		EP[:vFLECCS_output][y,i,t-1]-EP[:vFLECCS_output][y,i,t] <=dfGen_ccs[(dfGen_ccs[!,:R_ID].==y),:Ramp_Dn_Percentage][i]*dfGen_ccs[(dfGen_ccs[!,:R_ID].==y),:Cap_Size][i]*(vCOMMIT_FLECCS[y,i,t]-vSTART_FLECCS[y,i,t])
			-dfGen_ccs[(dfGen_ccs[!,:R_ID].==y),:Min_Power][i]*dfGen_ccs[(dfGen_ccs[!,:R_ID].==y),:Cap_Size][i]*vSTART_FLECCS[y,i,t]
			+min(1,max(dfGen_ccs[(dfGen_ccs[!,:R_ID].==y),:Min_Power][i],dfGen_ccs[!,:Ramp_Dn_Percentage][y]))*dfGen_ccs[(dfGen_ccs[!,:R_ID].==y),:Cap_Size][i]*vSHUT_FLECCS[y,i,t])

	### Minimum and maximum power output constraints (Constraints #7-8)
	@constraints(EP, begin
		# Minimum stable power generated per technology "y" at hour "t" > Min power
		[y in FLECCS_ALL, i in COMMIT_ccs, t=1:T], EP[:vFLECCS_output][y,i,t] >= dfGen_ccs[(dfGen_ccs[!,:R_ID].==y),:Min_Power][i]*dfGen_ccs[(dfGen_ccs[!,:R_ID].==y),:Cap_Size][i]*vCOMMIT_FLECCS[y,i,t]

		# Maximum power generated per technology "y" at hour "t" < Max power
		[y in FLECCS_ALL, i in COMMIT_ccs, t=1:T], EP[:vFLECCS_output][y,i,t] <= dfGen_ccs[(dfGen_ccs[!,:R_ID].==y),:Cap_Size][i]*vCOMMIT_FLECCS[y,i,t]
	end)


    ####Fangwei
	### Minimum up and down times (Constraints #9-10)
	for y in FLECCS_ALL
		for i in COMMIT_ccs
		    ## up time
		    Up_Time = Int(floor(dfGen_ccs[(dfGen_ccs[!,:R_ID].==y),:Up_Time][i]))
		    Up_Time_HOURS = [] # Set of hours in the summation term of the maximum up time constraint for the first subperiod of each representative period
		    for s in START_SUBPERIODS
			    Up_Time_HOURS = union(Up_Time_HOURS, (s+1):(s+Up_Time-1))
		    end

		    @constraints(EP, begin
			    # cUpTimeInterior: Constraint looks back over last n hours, where n = dfGen_ccs[(dfGen_ccs[!,:R_ID].==y),:Up_Time][i])
			    [t in setdiff(INTERIOR_SUBPERIODS,Up_Time_HOURS)], vCOMMIT_FLECCS[y,i,t] >= sum(vSTART_FLECCS[y,i,e] for e=(t-dfGen_ccs[(dfGen_ccs[!,:R_ID].==y),:Up_Time][i]):t)

			    # cUpTimeWrap: If n is greater than the number of subperiods left in the period, constraint wraps around to first hour of time series
			    # cUpTimeWrap constraint equivalant to: sum(vSTART_FLECCS[y,e] for e=(t-((t%hours_per_subperiod)-1):t))+sum(vSTART_FLECCS[y,e] for e=(hours_per_subperiod_max-(dfGen_ccs[(dfGen_ccs[!,:R_ID].==y),:Up_Time][i])-(t%hours_per_subperiod))):hours_per_subperiod_max)
			    [t in Up_Time_HOURS], vCOMMIT_FLECCS[y,i,t] >= sum(vSTART_FLECCS[y,i,e] for e=(t-((t%hours_per_subperiod)-1):t))+sum(vSTART_FLECCS[y,i,e] for e=((t+hours_per_subperiod-(t%hours_per_subperiod))-(dfGen_ccs[(dfGen_ccs[!,:R_ID].==y),:Up_Time][i]-(t%hours_per_subperiod))):(t+hours_per_subperiod-(t%hours_per_subperiod)))

			    # cUpTimeStart:
			    # NOTE: Expression t+hours_per_subperiod-(t%hours_per_subperiod) is equivalant to "hours_per_subperiod_max"
			    [t in START_SUBPERIODS], vCOMMIT_FLECCS[y,i,t] >= vSTART_FLECCS[y,i,t]+sum(vSTART_FLECCS[y,i,e] for e=((t+hours_per_subperiod-1)-(dfGen_ccs[(dfGen_ccs[!,:R_ID].==y),:Up_Time][i]-1)):(t+hours_per_subperiod-1))
		    end)

		    ## down time
		    Down_Time = Int(floor(dfGen_ccs[(dfGen_ccs[!,:R_ID].==y),:Down_Time][i]))
		    Down_Time_HOURS = [] # Set of hours in the summation term of the maximum down time constraint for the first subperiod of each representative period
		    for s in START_SUBPERIODS
		    	Down_Time_HOURS = union(Down_Time_HOURS, (s+1):(s+Down_Time-1))
		    end

		    # Constraint looks back over last n hours, where n = dfGen_ccs[(dfGen_ccs[!,:R_ID].==y),:Down_Time][i]
		    # TODO: Replace LHS of constraints in this block with eNumPlantsOffline[y,t]
		    @constraints(EP, begin
		    	# cDownTimeInterior: Constraint looks back over last n hours, where n = inputs["pDMS_Time"][y]
		    	[t in setdiff(INTERIOR_SUBPERIODS,Down_Time_HOURS)], EP[:eTotalCapFLECCS][y,i]/dfGen_ccs[(dfGen_ccs[!,:R_ID].==y),:Cap_Size][i]-vCOMMIT_FLECCS[y,i,t] >= sum(vSHUT_FLECCS[y,i,e] for e=(t-dfGen_ccs[(dfGen_ccs[!,:R_ID].==y),:Down_Time][i]):t)

		    	# cDownTimeWrap: If n is greater than the number of subperiods left in the period, constraint wraps around to first hour of time series
		    	# cDownTimeWrap constraint equivalant to: EP[:eTotalCapFLECCS][y,i]/dfGen_ccs[(dfGen_ccs[!,:R_ID].==y),:Cap_Size][i]-vCOMMIT_FLECCS[y,t] >= sum(vSHUT_FLECCS[y,e] for e=(t-((t%hours_per_subperiod)-1):t))+sum(vSHUT_FLECCS[y,e] for e=(hours_per_subperiod_max-(dfGen_ccs[(dfGen_ccs[!,:R_ID].==y),:Down_Time][i]-(t%hours_per_subperiod))):hours_per_subperiod_max)
		    	[t in Down_Time_HOURS], EP[:eTotalCapFLECCS][y,i]/dfGen_ccs[(dfGen_ccs[!,:R_ID].==y),:Cap_Size][i]-vCOMMIT_FLECCS[y,i,t] >= sum(vSHUT_FLECCS[y,i,e] for e=(t-((t%hours_per_subperiod)-1):t))+sum(vSHUT_FLECCS[y,i,e] for e=((t+hours_per_subperiod-(t%hours_per_subperiod))-(dfGen_ccs[(dfGen_ccs[!,:R_ID].==y),:Down_Time][i]-(t%hours_per_subperiod))):(t+hours_per_subperiod-(t%hours_per_subperiod)))
    
		    	# cDownTimeStart:
		    	# NOTE: Expression t+hours_per_subperiod-(t%hours_per_subperiod) is equivalant to "hours_per_subperiod_max"
		    	[t in START_SUBPERIODS], EP[:eTotalCapFLECCS][y,i]/dfGen_ccs[(dfGen_ccs[!,:R_ID].==y),:Cap_Size][i]-vCOMMIT_FLECCS[y,i,t]  >= vSHUT_FLECCS[y,i,t]+sum(vSHUT_FLECCS[y,i,e] for e=((t+hours_per_subperiod-1)-(dfGen_ccs[(dfGen_ccs[!,:R_ID].==y),:Down_Time][i]-1)):(t+hours_per_subperiod-1))
		    end)
	    end
	end

	## END Constraints for thermal units subject to integer (discrete) unit commitment decisions



	return EP
end

