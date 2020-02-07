## Z update in ADMM, to project X in the ignored constraint sets in X update

# compile the common constraints for Slave Z, to save time
function DefineSlaveZ(u_consumer_initial)
    # define model slave Z, use the same solver as slave x
    slaveZ = Model(solver = SlaveXSolverChoice)
    # define decision variables, with suffix _Z
    @variables slaveZ begin
        p[OptionsSet] >= 0   # price in the menu
        r[OptionsSet] >= 0   # reliability in the menu
        z_consumer[ConsumersSet, OptionsSet] >= 0  # auxilary variable for r*u
        y_consumer[ConsumersSet, OptionsSet] >= 0  # auxilary variable for p*u
        surplus[ConsumersSet] >= 0  # gamma in the note
        supplyR_Z[OptionsSet, TimeSet,ScenariosSet] >= 0   # total supply to residential consuemrs in each option
        totalCost_Z[TimeSet,ScenariosSet] >= 0
    end

    # promised reliability is delivered, the aggregator is only aware of system profile
for option = OptionsSet
    Expression_1 = AffExpr()
    for consumer = ConsumersSet, t = TimeSet
        append!(Expression_1, z_consumer[consumer, option]/zSCALE*Dl[consumer]*SystemProfile[t])
    end
    Expression_2 = AffExpr()
    for scenario = ScenariosSet
        for t = TimeSet
            append!(Expression_2, supplyR_Z[option,t, scenario]/sSCALE[option]*ScenarioProb[scenario])
        end
    end
    @constraints(slaveZ, begin
        Expression_1 == Expression_2
    end)
end


ProfitsExpression_1 = AffExpr()
 for option = OptionsSet, consumer = ConsumersSet
     append!(ProfitsExpression_1, y_consumer[consumer,option]/ySCALE*Dl[consumer])
 end
 ProfitsExpression_2 = AffExpr()
 for scenario = ScenariosSet
     for t = TimeSet
         append!(ProfitsExpression_2, totalCost_Z[t,scenario]/cSCALE*ScenarioProb[scenario])
     end
 end
 @constraints(slaveZ, begin
     ProfitsExpression_1*sum(SystemProfile) -  ProfitsExpression_2 ==  Profits
 end)


    ## define common constraints
    @constraints(slaveZ, begin
        # optimality conditions
        [consumer = ConsumersSet], surplus[consumer] <=
            sum((Vl[consumer]*z_consumer[consumer, option]/zSCALE - y_consumer[consumer,option]/ySCALE)
                for option = OptionsSet)
    end)

    for option = OptionsSet
        for consumer = ConsumersSet
            @constraints(slaveZ, begin
                surplus[consumer] >= Vl[consumer]*r[option]-p[option]
                z_consumer[consumer, option]/zSCALE == r[option]*u_consumer_initial[consumer, option]
                y_consumer[consumer, option]/ySCALE == p[option]*u_consumer_initial[consumer, option]
                # McCormick envelopes
                #=
                z_consumer[consumer, option]/zSCALE <= u_consumer_initial[consumer, option]
                z_consumer[consumer, option]/zSCALE <= r[option]
                z_consumer[consumer, option]/zSCALE >= u_consumer_initial[consumer, option] + r[option] - 1
                y_consumer[consumer, option]/ySCALE <= u_consumer_initial[consumer, option]*PriceMax
                y_consumer[consumer, option]/ySCALE <= p[option]
                y_consumer[consumer, option]/ySCALE >= u_consumer_initial[consumer, option]*PriceMax + p[option] - PriceMax
                =#
            end)

        end
    end

    for consumer = 1:NConsumers-1
        for optionStart = OptionsSet
            @constraints(slaveZ, begin
                # optimalizty cuts on z_consumer and y_consumer
                sum(z_consumer[consumer, option] for option = optionStart:NOptions)  <=
                    sum(z_consumer[consumer+1, option] for option = optionStart:NOptions)
                sum(y_consumer[consumer, option] for option = optionStart:NOptions)  <=
                    sum(y_consumer[consumer+1, option] for option = optionStart:NOptions)
            end)
        end
        @constraints(slaveZ, begin
            # surplus of each consumer should also be increasing
            surplus[consumer] <= surplus[consumer+1]
        end)
    end

    @constraints(slaveZ, begin
        # price and reliability in the menu should be increasing
        [option = 1:NOptions-1],
           r[option] <= r[option+1]
        [option = 1:NOptions-1],
            p[option] <= p[option+1]
        r[NOptions] <= 1
    end)


    #println("Compiling Time of DefineSlaveZ: ", round(toq()))
    # return the prefined model, to be used in SolveSlaveZ
    return slaveZ, p, r, z_consumer, y_consumer, surplus, supplyR_Z, totalCost_Z

end



# change the objective function of slaveZ at each iterations
function SolveSlaveZ(slaveZ, p, r, z_consumer, y_consumer, surplus, supplyR_Z, totalCost_Z,
    totalCost_X_sol, supplyR_X_sol, totalCost_Z_sol, supplyR_Z_sol, niterationsADMM)

    #tic()
    # the objective needs to be updated at each iteration of ADMM
    # use append to add up terms, more efficient
    objExpression = QuadExpr()

    for scenario = ScenariosSet, t = TimeSet
        append!(objExpression, ScenarioProb[scenario]*(totalCost_Z[t,scenario]/cSCALE - totalCost_X_sol[t,scenario]/cSCALE -totalCost_U[t,scenario]/cSCALE)^2)
    end

    for scenario = ScenariosSet, t = TimeSet, option = OptionsSet
        #append!(objExpression, ScenarioProb[scenario]*(supplyR_Z[option,t,scenario]-
        #    supplyR_X_sol[option,t,scenario]- supplyR_U[option,t,scenario])^2)
            append!(objExpression, ScenarioProb[scenario]*(supplyR_Z[option,t,scenario]/sSCALE[option]-
                supplyR_X_sol[option,t,scenario]/sSCALE[option]- supplyR_U[option,t,scenario]/sSCALE[option])^2)
    end


    @objective(slaveZ, Min, 1000*objExpression)
    #println("Compiling Slave_Z takes: ", round(toq(),2))

    #tic()
    # solve it
    status_slaveZ = solve(slaveZ)
    #println("Solving Slave_Z takes: ", round(toq(),2) )
    #println("SlaveZ: ", status_slaveZ,
    #    "; ObjValue: ", round(getobjectivevalue(slaveZ)))
    println("*****Z Update is solved")
    # get solutions
    totalCost_Z_current_sol = getvalue(totalCost_Z[:,:])
    #println(totalCost_Z_current_sol)
    supplyR_Z_current_sol = getvalue(supplyR_Z[:,:,:])
    #println(supplyR_Z_current_sol)

    #r_sol = getvalue(r[:])
    #p_sol = getvalue(p[:])
    y_consumer_sol = getvalue(y_consumer[:,:])
    ## calculate dual residual -RHO(Z_k - Z_{k-1})
    dualResidual_temp= sum((totalCost_Z_current_sol - totalCost_Z_sol).^2) +
                    sum((supplyR_Z_current_sol- supplyR_Z_sol).^2)
    # not multiplied by RHO, to be used in relative dualResidual
    dualResidual = dualResidual_temp^0.5


    ## update values of Z
    totalCost_Z_sol = totalCost_Z_current_sol
    supplyR_Z_sol = supplyR_Z_current_sol

    ## calculate primal Residual X_k - Z_k
    primalResidual_temp= sum((totalCost_X_sol- totalCost_Z_sol).^2) +
                        sum((supplyR_X_sol - supplyR_Z_sol).^2)

    primalResidual = primalResidual_temp^0.5

    ## calculate relative primalResidual
    primalDtemp1 = sum(totalCost_X_sol.^2) + sum(supplyR_X_sol.^2)
    primalDtemp2 = sum(totalCost_Z_sol.^2) + sum(supplyR_Z_sol.^2)
    primalDtempMax = max(primalDtemp1^0.5, primalDtemp2^0.5)
    primalResidual_rel= primalResidual/primalDtempMax*100  # in percentage
    #println("**** primalResidual_rel: ", round(primalResidual_rel,4),"%")
    flush(STDOUT)
    return totalCost_Z_sol, supplyR_Z_sol,y_consumer_sol, primalResidual_rel, dualResidual

end
