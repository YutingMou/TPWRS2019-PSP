## this is used to solve the slave problem parallelly,
## each subproblem is UC
## decouple the time-coupled constraints (reliability and profits)


function SolveSlaveX(inputSlaveX)
    tic()

    # input data
    subproblemIdx = inputSlaveX[1]  # index of this subproblem
    MonthIdx = inputSlaveX[2]
    ScenarioIdx =  inputSlaveX[3]
    niterationsDD = inputSlaveX[4]
    dual_vb  = inputSlaveX[5]
    dual_reliability = inputSlaveX[6]
    dual_profits = inputSlaveX[7]

    TdaysX = DaysX[MonthIdx]
    ThoursX = 24*TdaysX
    TimeSetX = 1:ThoursX

    # define the model
    slaveX = Model(solver = SlaveXMIPSolverChoice)


    ## define decision varioables
    @variables slaveX begin
        # conventional generators
        u[GeneratorsSet, TimeSetX] , Bin   # on/off binary
        v[GeneratorsSet, TimeSetX] , Bin   # start-up binary
        z[GeneratorsSet, TimeSetX] , Bin   # shut-down binary
        production[GeneratorsSet, TimeSetX] >= 0
        productionCost[GeneratorsSet, TimeSetX] >= 0
        startupCost0[GeneratorsSet, TimeSetX] >=0    # startup cost + startup fuel cost
        miniLoadCost[GeneratorsSet, TimeSetX] >=0    # regardless of the load, this cost is applicable if the unit is on
        totalCost[TimeSetX] >= 0 # total cost of all generators at t + loadshedding cost
        # pumped hydro storage
        0 <= s[TimeSetX] <= PSEnergy             # energy in pumped storage in MWH
        0 <= p_pump[TimeSetX] <= PSPumpingMax    # power pumped into pumped storage
        0 <= p_prod[TimeSetX] <= PSProducingMax  # power production of pumped sorage
        0 <= loadShedding[TimeSetX] <= sum(Di)  # in case of IC demand is too high

        renewableShedding[TimeSetX] >= 0
        residentialD[TimeSetX] >= 0   # total supply to residential consuemrs
        supplyR[OptionsSet, TimeSetX] >= 0   # total supply to residential consuemrs in each option
    end


    println("***** variables are defined")



    @constraints(slaveX, begin
        # PS empty everyday at mightnight
        [d = 1:TdaysX], s[d*24] <= 1
    end)

    ## PS and UC
    @constraints(slaveX, begin
    # initial of pumped hydro
        s[1] == 0 + p_pump[1]*PSEff - p_prod[1]
    end)
    for t = 2:ThoursX
        @constraints(slaveX, begin
            s[t] == s[t-1] + p_pump[t]*PSEff - p_prod[t]
        end)
    end

    for t = TimeSetX
        for g = GeneratorsSet
            @constraints(slaveX, begin
            # cycling UC
                u[g, time_cycling(t-1, 1, ThoursX)]+v[g, t] - z[g, t] == u[g, t]
                sum(v[g, time_cycling(y, 1, ThoursX)]
                    for y = t - Generators[:UT][g] + 1:t ) <= u[g, t]
                sum(v[g,time_cycling(yy, 1, ThoursX)]
                    for yy = t + 1 : t + Generators[:DT][g]) <= 1 - u[g, t]
                # conventional gnerators' production bounds
                production[g,t] >= MinRunCapacityDerated[g, StartHourX[MonthIdx]+t-1]*u[g,t]
                production[g,t] <= MaxRunCapacityDerated[g, StartHourX[MonthIdx]+t-1]*u[g,t]
            end)
        end
    end



    for t = TimeSetX
        # supply cannot exceeds consumer Subscription
        for option = OptionsSet
            @constraints(slaveX, begin
                supplyR[option,t]/sSCALE[option] <=
                    Di[option]*SystemProfile[StartHourX[MonthIdx]+t-1]
            end)
        end

        @constraints(slaveX, begin
            # total supply to all options = residential demand
            residentialD[t] == sum(supplyR[option,t]/sSCALE[option] for option = OptionsSet)
            # supply demand balance
            sum(production[g, t] for g = GeneratorsSet)  + loadShedding[t] +
                WindProfiles[StartHourX[MonthIdx]+t-1,ScenarioIdx] - renewableShedding[t] +
			         SolarProfiles[StartHourX[MonthIdx]+t-1,ScenarioIdx] +
                     HistoricalProfiles[:Import][StartHourX[MonthIdx]+t-1] +
                     p_prod[t] - p_pump[t] + 77 ==
				residentialD[t] + HistoricalProfiles[:ICDemand][StartHourX[MonthIdx]+t-1]
        end)

    end

    ## model the costs of conventional generators
    for g = GeneratorsSet
        ## collect heat rate curve and fuel price data for generator g
    	hrc = view(HeatRateCurves, HeatRateCurves[:SeriesHRC] .==
    		Generators[:HRSGenerator][g])
    	fprows = .&(FuelPrice[:Fuel] .== Generators[:FuelGenerator][g],
    		FuelPrice[:FuelHub] .== Generators[:FHGenerator][g])
        ## cost curve is piecewise Linear
        FuelPrice_g = view(FuelPrice[:Price], fprows)[1]
        for t = TimeSetX
            #for i = 1:numhr
            @constraints(slaveX, begin
                # productionCost
                productionCost[g,t] .>=
                FuelPrice_g*(hrc[:SlopeHRC]*production[g,t] + hrc[:InterceptHRC]*u[g,t])
                # miniLoadCost
                miniLoadCost[g,t] ==
                Generators[:NoLoadConsumption][g]*MinRunCapacityDerated[g,StartHourX[MonthIdx]+t-1]*FuelPrice_g*u[g, t]
                # startup costs
                startupCost0[g,t] ==
                (Generators[:StartupCost][g] + Generators[:StartupFuel][g]*FuelPrice_g)*MinRunCapacityDerated[g,StartHourX[MonthIdx]+t-1]*v[g, t]
            end)
            #end
        end
    end

    @constraints(slaveX, begin
        [t = TimeSetX],
        totalCost[t] == sum( (productionCost[g, t] + startupCost0[g,t] + miniLoadCost[g,t])
            for g = GeneratorsSet) + loadShedding[t]*VOLL
    end)

    ### creat objective function in a more efficient way

    # costs
    objExpression_1 = AffExpr()
    for t =  TimeSetX
        append!(objExpression_1, totalCost[t])
    end

    # benefits
    objExpression_2 = AffExpr()
    for t = TimeSetX
        append!(objExpression_2,
         sum(supplyR[option,t]*(VB[option]+VB[option+1])*0.5/sSCALE[option] for option = OptionsSet)
         )
    end

    # the  term related to reliability
    objExpression_3 =  AffExpr()
    for t = TimeSetX,  option = OptionsSet
        append!(objExpression_3, dual_reliability[option]*supplyR[option,t]/sSCALE[option])
    end

    # the term related to profits
    # dual profits can not exceed 1, otherwise, the problem is unbounded
    # i.e., it becomes maximize totalcost
    objExpression_4 =  AffExpr()
    for  t = TimeSetX
        append!(objExpression_4, totalCost[t])
    end
    objExpression_4 = dual_profits*objExpression_4

    objExpression = objExpression_1 - objExpression_2 -
        objExpression_3 - objExpression_4
    @objective(slaveX, Min, objExpression)

    toq()
    #println("Compiling SlaveX ", subproblemIdx, " takes: ", round(toq(),2))
    #flush(STDOUT)
    tic()
    # solve it
    status_slaveX = solve(slaveX;  suppress_warnings=true, relaxation=true)
    if status_slaveX == :Optimal

        # get objvalue
        currentslaveXobjval = getobjectivevalue(slaveX)
        # get solutions
        totalCost_current_sol = getvalue(totalCost[:])
        supplyR_current_sol = getvalue(supplyR[:,:])
        println("SlaveX ", subproblemIdx, ": ", status_slaveX, "; ObjValue: ", round(currentslaveXobjval), "; It takes: ", round(toq(),2))
        flush(STDOUT)

        # prepare the results to be output at the last iteration of DD
        if false #niterationsDD < MaxIterationDD
            # return objective value and dual, to be used to construct Benders cuts
            return subproblemIdx, totalCost_current_sol, supplyR_current_sol
        else    # the last iteration, output the results
            u_sol = Int.(round.(getvalue(u[:,:])))
            v_sol = Int.(round.(getvalue(v[:,:])))
            z_sol = Int.(round.(getvalue(z[:,:])))
            production_sol = getvalue(production[:,:])
            productionCost_sol = getvalue(productionCost[:,:])
            startupCost0_sol = getvalue(startupCost0[:,:])
            miniLoadCost_sol = getvalue(miniLoadCost[:,:])
            s_sol = getvalue(s[:])
            p_pump_sol = getvalue(p_pump[:])
            p_prod_sol = getvalue(p_prod[:])
            loadShedding_sol = getvalue(loadShedding[:])
            renewableShedding_sol = getvalue(renewableShedding[:])

            return  subproblemIdx, totalCost_current_sol, supplyR_current_sol, u_sol,
                    v_sol, z_sol, production_sol, productionCost_sol,
                    startupCost0_sol, miniLoadCost_sol, s_sol, p_pump_sol,
                    p_prod_sol, loadShedding_sol, renewableShedding_sol, currentslaveXobjval

        end
    else
        println("SlaveX ", subproblemIdx, ": ", status_slaveX)
    end

end
