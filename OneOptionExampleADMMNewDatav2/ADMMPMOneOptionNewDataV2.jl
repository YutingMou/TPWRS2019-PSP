
function ADMMPM(niterationsADMM, rho, supplyR_Z_sol, totalCost_Z_sol, supplyR_U_sol, totalCost_U_sol)


        # producer model
        PM = Model(solver = GurobiSolver(
                    OutputFlag=0,
                    MIPGap=0.001,
                    TimeLimit=3600
                    )
            )

        # define variables
        @variables PM begin
            production[GeneratorsSet, TimeSet] >= 0
            productionCost[GeneratorsSet, TimeSet] >= 0
            supplyR[OptionsSet, TimeSet] >= 0   # total supply to residential consuemrs in each option
            totalCost[TimeSet] >= 0
        end

        # define constraints
        # productionCapacity limit
        for g = GeneratorsSet, t = TimeSet
            @constraints(PM, begin
                production[g,t] <= Pmax[g,t]
                productionCost[g,t] == production[g,t]*MC[g]
            end)
        end
        # totalCost is the sum of productionCost
        for t = TimeSet
            @constraints(PM, begin
                totalCost[t] * TotalCostScale == sum(productionCost[g,t] for g = GeneratorsSet)
            end)
        end
        # supply demand balance
        for t = TimeSet
            @constraints(PM, begin
                sum(supplyR[option, t] for option = OptionsSet) == sum(production[g, t] for g = GeneratorsSet)
            end)
        end

        # supply can not exceed subscription
        for t = TimeSet, option = 1
            @constraints(PM, begin
                supplyR[option, t] <=  D
            end)
        end

        objPM_expression1 = sum(productionCost[g, t]*Omega[t] for g = GeneratorsSet, t = TimeSet) -
            sum((supplyR[1,t]*(VB[1] + VB[2]))*0.5*Omega[t]   for t = TimeSet)

        objPM_expression2_1 = sum((totalCost[t] - totalCost_Z_sol[t] + totalCost_U_sol[t])^2 for t = TimeSet)
        objPM_expression2_2 = sum((supplyR[option, t] - supplyR_Z_sol[option,t] + supplyR_U_sol[option,t])^2
                for option = OptionsSet, t = TimeSet)
        objPM = objPM_expression1 + rho[1]/2*objPM_expression2_1 + rho[2]/2*objPM_expression2_2

        @objective(PM, Min,  objPM)

        println("solving producer model")

        status_PM = solve(PM)
        if status_PM == :Optimal

            supplyR_X_sol = getvalue(supplyR[:,:])
            totalCost_X_sol = getvalue(totalCost[:])
            objvalue_PM = getobjectivevalue(PM)
            production_X_sol = getvalue(production[:,:])
            objPM_expression1_sol = getvalue(objPM_expression1)
            objPM_expression2_1_sol = getvalue(objPM_expression2_1)*rho[1]/2
            objPM_expression2_2_sol = getvalue(objPM_expression2_2)*rho[2]/2

            #println("rho: ", rho)

            println("*********Output From the producer model")
            println("objvalue_PM: ", objvalue_PM)
            println("objPM_expression1_sol: ", objPM_expression1_sol)
            println("objPM_expression2_1_sol: ", objPM_expression2_1_sol)
            println("objPM_expression2_2_sol: ", objPM_expression2_2_sol)

            println("production_X_sol: ", production_X_sol)
            println("supplyR_Z_sol: ", supplyR_Z_sol)
            println("totalCost_Z_sol: ", totalCost_Z_sol)
            println("supplyR_U_sol: ", supplyR_U_sol)
            println("totalCost_U_sol: ", totalCost_U_sol)
            #println("supplyR_X_sol: ", supplyR_X_sol)
            #println("totalCost_X_sol: ", totalCost_X_sol)

            #costs = sum(MC[g] * production_sol[g, t] * Omega[t] for g = GeneratorsSet, t = TimeSet)
            #socialwelfare = benefits - costs
            #println("subR_sol: ", subR_sol)
            #println("supplyR_X_sol: ", supplyR_X_sol)
            #println("totalCost_X_sol: ", totalCost_X_sol)


            #println("profits: ", profits)
            #println("socialwelfare: ", socialwelfare)
        else
            println(status_PM)
        end

        return objvalue_PM, supplyR_X_sol, totalCost_X_sol, production_X_sol

    end
