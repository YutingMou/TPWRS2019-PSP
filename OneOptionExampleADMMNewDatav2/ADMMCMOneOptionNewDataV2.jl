### dual decomposition cutting plane consumer subproblem
function ADMMCM(rho, supplyR_X_sol, totalCost_X_sol, supplyR_U_sol, totalCost_U_sol)

        # consumer Model
        CM = Model(solver = GurobiSolver(
                    OutputFlag=0,
                    MIPGap=0.001,
                    TimeLimit=3600
                    )
            )

        @variables CM begin
            p[OptionsSet] >= 0   # price in the menu
            r[OptionsSet] >= 0   # reliability in the menu
            u_consumer[ConsumersSet, OptionsSet], Bin   # consumers choice, binary
            z_consumer[ConsumersSet, OptionsSet] >= 0  # auxilary variable for r*u
            y_consumer[ConsumersSet, OptionsSet] >= 0  # auxilary variable for p*u
            surplus[ConsumersSet] >= 0  # gamma in the note
            supplyR[OptionsSet, TimeSet] >= 0   # total supply to residential consuemrs in each option
            totalCost[TimeSet] >= 0
        end

        # consumers choose at most one option
        @constraints(CM, begin
            [consumer = ConsumersSet],
                sum(u_consumer[consumer, option] for option = OptionsSet) <= 1
        end)


        # constraints on gamma for each consumer, i.e., lower level optimality condition
        for consumer = ConsumersSet, option = OptionsSet
            @constraints(CM, begin
                surplus[consumer] >= Vl[consumer]*r[option]-p[option]
            end)
        end

        for  consumer = ConsumersSet
            @constraints(CM, begin
            surplus[consumer] <= sum(Vl[consumer]*z_consumer[consumer,option] - y_consumer[consumer,option]
             for option = OptionsSet)
            end)
        end

        @constraints(CM, begin
            sum(Dl*u_consumer[consumer, 1] for consumer = ConsumersSet) == D
        end)

        #  McCormick envelope
        for option = OptionsSet, consumer = ConsumersSet
            @constraints(CM, begin
                z_consumer[consumer, option] <= u_consumer[consumer, option]
                z_consumer[consumer, option] <= r[option]
                z_consumer[consumer, option] >= u_consumer[consumer, option] + r[option] - 1
                y_consumer[consumer, option] <= u_consumer[consumer, option]*PriceMax
                y_consumer[consumer, option] <= p[option]
                y_consumer[consumer, option] >= u_consumer[consumer, option]*PriceMax + p[option] - PriceMax

            end)
        end


        # valid cuts
        for consumer = 1:NConsumers-1
            for optionStart = OptionsSet
                @constraints(CM, begin
                    # Valid cuts on u_consumer
                    sum(u_consumer[consumer, option] for option = optionStart:NOptions)  <=
            		    sum(u_consumer[consumer+1, option] for option = optionStart:NOptions)
                    sum(y_consumer[consumer, option] for option = optionStart:NOptions)  <=
                        sum(y_consumer[consumer+1, option] for option = optionStart:NOptions)
                end)
            end
            @constraints(CM, begin
                # surplus of each consumer should also be increasing
                surplus[consumer] <= surplus[consumer+1]
            end)
        end


        # promised reliablity must be delivered
        for option = OptionsSet
            @constraints(CM, begin
            sum(z_consumer[consumer, option]*Dl*DemandScaling[t]/2
                for consumer = ConsumersSet, t = TimeSet) ==
            sum(supplyR[option,t]*Omega[t] for t = TimeSet)
            end)
        end

        @constraints(CM, begin
            sum(y_consumer[consumer, option]  for consumer = ConsumersSet, option = OptionsSet) -
            sum( totalCost[t] *TotalCostScale * Omega[t] for t = TimeSet) == Profits
        end)


        #

        @constraints(CM, begin
            r[1]*VB[1] == p[1]
        end)


        objCM_expression1 = sum((totalCost[t] - totalCost_X_sol[t] - totalCost_U_sol[t])^2 for t = TimeSet)
        objCM_expression2 = sum((supplyR[option, t] - supplyR_X_sol[option,t] - supplyR_U_sol[option,t])^2
                for option = OptionsSet, t = TimeSet)
        objCM = rho[1]/2*objCM_expression1 + rho[2]/2*objCM_expression2


        @objective(CM, Min,  objCM)

        println("solving consumer model")

        status_CM = solve(CM)
        if status_CM == :Optimal
            supplyR_Z_sol = getvalue(supplyR[:,:])
            totalCost_Z_sol = getvalue(totalCost[:])
            r_sol = getvalue(r[:])
            p_sol = getvalue(p[:])
            objvalue_CM = getobjectivevalue(CM)
            #println("r_sol: ",r_sol)
            #println("p_sol: ", p_sol)
            #println("u_consumer_sol_sum: ", u_consumer_sol_sum)
            #println("objvalue_CM: ", objvalue_CM)

        else
            println(status_CM)
        end

    return objvalue_CM, r_sol, p_sol,supplyR_Z_sol,totalCost_Z_sol
end
