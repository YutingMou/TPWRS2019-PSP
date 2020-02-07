## cutting plane master problem

master = Model(solver = GurobiSolver(
            OutputFlag=0,
            #MIPGap=0.001,
            TimeLimit=3600
            )
    )


# first formulation of
@variables master begin
     z_master <= -5e7
    #-1e2 <= dual_vb_master <= 1e2
    dual_profits_master <=1
    dual_reliability_master[OptionsSet]
#    -1 <=  dual_profits_master <= 1    # otherwise SlaveX will be unbounded
#    -1e2 <= dual_reliability_master[OptionsSet] <= 1e2
end



#=
@constraints(master, begin
    [option = 2:NOptions],
        dual_reliability_master[option] == 0
end)
=#

@objective(master, Max, z_master)
