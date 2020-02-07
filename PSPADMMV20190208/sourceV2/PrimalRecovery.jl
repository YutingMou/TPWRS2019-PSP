### used for primal recovery
### fix the solution from the producer model (cost and supply to each option)
function SlaveZPR(u_consumer_initial, totalCost_X_sol, supplyR_X_sol)

        # consumer Model
        slaveZPR = Model(solver = SlaveXSolverChoice)


        @variables slaveZPR begin
            0 <= p[OptionsSet] <= PriceMax  # price in the menu
            0 <= r[OptionsSet] <= 1   # reliability in the menu
            z_consumer[ConsumersSet, OptionsSet] >= 0  # auxilary variable for r*u
            y_consumer[ConsumersSet, OptionsSet] >= 0  # auxilary variable for p*u
            surplus[ConsumersSet] >= 0  # gamma in the note
            profit_slack   # slack variable to make the problem feasible forever
            profit_Abs >= 0   # absolute value of the slack variable
            r_slack[OptionsSet]
            r_Abs[OptionsSet] >= 0
        end


    # promised reliability is delivered, the aggregator is only aware of system profile
    for option = OptionsSet
        Expression_1 = AffExpr()
#=
        for consumer = ConsumersSet, t = TimeSet
            append!(Expression_1, z_consumer[consumer, option]/zSCALE*Dl[consumer]*SystemProfile[t])
        end
=#      # subscribed
        append!(Expression_1, r[option]*Di[option]*sum(SystemProfile))

        Expression_2 = AffExpr()
        for scenario = ScenariosSet
            for t = TimeSet
                append!(Expression_2, supplyR_X_sol[option,t, scenario]/sSCALE[option]*ScenarioProb[scenario])
            end
        end
        @constraints(slaveZPR, begin
            Expression_1 == Expression_2 + r_slack[option]
        end)
    end


    ProfitsExpression_1 = AffExpr()
    #=
     for option = OptionsSet, consumer = ConsumersSet
         append!(ProfitsExpression_1, y_consumer[consumer,option]/ySCALE*Dl[consumer])
     end
     =#
     for option = OptionsSet
         append!(ProfitsExpression_1, p[option]*Di[option])
     end
     ProfitsExpression_1 = ProfitsExpression_1*sum(SystemProfile)

     ProfitsExpression_2 = AffExpr()
     for scenario = ScenariosSet
         for t = TimeSet
             append!(ProfitsExpression_2, totalCost_X_sol[t,scenario]/cSCALE*ScenarioProb[scenario])
         end
     end
     @constraints(slaveZPR, begin
         ProfitsExpression_1 -  ProfitsExpression_2 ==  Profits + profit_slack
     end)


        ## define common constraints
        @constraints(slaveZPR, begin
            # optimality conditions
            [consumer = ConsumersSet], surplus[consumer] <=
                sum((Vl[consumer]*z_consumer[consumer, option]/zSCALE - y_consumer[consumer,option]/ySCALE)
                    for option = OptionsSet)
        end)

        for option = OptionsSet
            for consumer = ConsumersSet
                @constraints(slaveZPR, begin
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
                @constraints(slaveZPR, begin
                    # optimalizty cuts on z_consumer and y_consumer
                    sum(z_consumer[consumer, option] for option = optionStart:NOptions)  <=
                        sum(z_consumer[consumer+1, option] for option = optionStart:NOptions)
                    sum(y_consumer[consumer, option] for option = optionStart:NOptions)  <=
                        sum(y_consumer[consumer+1, option] for option = optionStart:NOptions)
                end)
            end
            @constraints(slaveZPR, begin
                # surplus of each consumer should also be increasing
                surplus[consumer] <= surplus[consumer+1]
            end)
        end

        @constraints(slaveZPR, begin
            # price and reliability in the menu should be increasing
            [option = 1:NOptions-1],
               r[option] <= r[option+1]
            [option = 1:NOptions-1],
                p[option] <= p[option+1]
            r[NOptions] <= 1
        end)

        @constraints(slaveZPR, begin
            profit_Abs >= profit_slack
            profit_Abs >= - profit_slack
        end)
    for option = OptionsSet
        @constraints(slaveZPR, begin
            r_Abs[option] >= r_slack[option]
            r_Abs[option] >= - r_slack[option]
        end)
    end


        objslaveZPR = Vmax*(profit_Abs + sum(r_Abs[option] for option = OptionsSet))

        @objective(slaveZPR, Min,  objslaveZPR)

        println("## solving slaveZPR")

        status_slaveZPR = solve(slaveZPR)
        if status_slaveZPR == :Optimal
            r_sol = getvalue(r[:])
            p_sol = getvalue(p[:])
            objvalue_slaveZPR = getobjectivevalue(slaveZPR)
            y_consumer_sol = getvalue(y_consumer[:,:])
            z_consumer_sol = getvalue(z_consumer[:,:])
            #println("r_sol: ",r_sol)
            #println("p_sol: ", p_sol)
            #println("u_consumer_sol_sum: ", u_consumer_sol_sum)
            #println("objvalue_slaveZPR: ", objvalue_slaveZPR)

        else
            println(status_slaveZPR)
        end

    return objvalue_slaveZPR, r_sol, p_sol#, y_consumer_sol, z_consumer_sol
end
