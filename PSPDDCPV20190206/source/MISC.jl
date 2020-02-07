# function to wrap around the unit commitment coupling constraints
function time_cycling(t:: Int64, Tstart:: Int64, Tend:: Int64)
	if t < Tstart
		t = Tend + rem(t-Tstart+1, Tend-Tstart+1)
	elseif t > Tend
        t = Tstart-1 + mod(t-Tstart+1, Tend-Tstart+1);
	end
    return t
end

# return month index  and scenario index based on subproblem index
function returnIDX(subproblemIdx::Int64, NSubhoriozns::Int64)
	if rem(subproblemIdx,NSubhoriozns) == 0
		MonthIdx = NSubhoriozns
		ScenarioIdx = div(subproblemIdx,NSubhoriozns)

	else
		MonthIdx = rem(subproblemIdx,NSubhoriozns)
		ScenarioIdx = div(subproblemIdx,NSubhoriozns) +1
	end
	return MonthIdx,ScenarioIdx
end

# projected reliability, assuming consumers follow system profile
function projectedReliability(u_consumer_sol, supplyR_sol)
    # subscription profile of each option
   subsProfileSystem = Array{Float64}(NOptions, Thours)
    for option = OptionsSet
        subsProfileSystem[option, :] = sum(
            u_consumer_sol[consumer, option]*Dl[consumer]*SystemProfile'
        for consumer = ConsumersSet)
    end
	r_projected =sum(sum(supplyR_X_sol[:,:,scenario],2)*ScenarioProb[scenario] for scenario = ScenariosSet)./sum(subsProfileSystem,2)
	r_projected = r_projected[:]
	return r_projected
end
