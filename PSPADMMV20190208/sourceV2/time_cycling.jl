## function to wrap around the unit commitment coupling constraints

function time_cycling(t:: Int64, Tstart:: Int64, Tend:: Int64)
	if t < Tstart
		t = Tend + rem(t-Tstart+1, Tend-Tstart+1)
	elseif t > Tend
        t = Tstart-1 + mod(t-Tstart+1, Tend-Tstart+1);
	end
    return t
end

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
