## The consumer model

# compile the common constraints for Slave Z, to save time
function DefineSlaveZ()
    #tic()
    # define model slave Z, use the same solver as slave x
    slaveZ = Model(solver = SlaveXMIPSolverChoice)
    # define decision variables, with suffix _Z
    @variables slaveZ begin
        p[OptionsSet] >= 0   # price in the menu
        r[OptionsSet] >= 0   # reliability in the menu
        u_consumer[ConsumersSet, OptionsSet], Bin
        z_consumer[ConsumersSet, OptionsSet] >= 0  # auxilary variable for r*u
        y_consumer[ConsumersSet, OptionsSet] >= 0  # auxilary variable for p*u
        surplus[ConsumersSet] >= 0  # gamma in the note
        0 <= vb_c <= Vmax
    end

#=
    @constraints(slaveZ, begin
        p[1] == r[1]*VB[1]
    end)
    for option = 2:NOptions
        @constraints(slaveZ, begin
            p[option] == p[option-1] + VB[option]*(r[option] - r[option-1])
        end)

    end
=#

#=
    # each consumer chooses at most one option
    @constraints(slaveZ, begin
        [consumer = ConsumersSet],
        sum(u_consumer[consumer, option] for option = OptionsSet) <= 1
    end)

    # total subscription to a certain option can  not exceed Di
    @constraints(slaveZ, begin
        [option = OptionsSet],
        sum(Dl[consumer]*u_consumer[consumer,option] for consumer = ConsumersSet) <= Di[option]
    end)


    @constraints(slaveZ, begin
        [option = 2:NOptions],
        sum(Dl[consumer]*u_consumer[consumer,option] for consumer = ConsumersSet) == Di[option]
    end)
=#
    ## define common constraints
    @constraints(slaveZ, begin
        # optimality conditions
        [consumer = ConsumersSet], surplus[consumer] <=
            sum((Vl[consumer]*z_consumer[consumer, option]/zSCALE - y_consumer[consumer,option]/ySCALE)
                for option = OptionsSet)
    end)

    for consumer = 1:NConsumers-1
        #=
        for optionStart = OptionsSet
            @constraints(slaveZ, begin
                # Valid cuts on z_consumer and y_consumer
                sum(u_consumer[consumer, option] for option = optionStart:NOptions)  <=
                    sum(u_consumer[consumer+1, option] for option = optionStart:NOptions)
                sum(z_consumer[consumer, option] for option = optionStart:NOptions)  <=
                    sum(z_consumer[consumer+1, option] for option = optionStart:NOptions)
                sum(y_consumer[consumer, option] for option = optionStart:NOptions)  <=
                    sum(y_consumer[consumer+1, option] for option = optionStart:NOptions)
            end)
        end
        =#
        @constraints(slaveZ, begin
            # surplus of each consumer should also be increasing
            surplus[consumer] <= surplus[consumer+1]
        end)
    end

    #println("Compiling Time of DefineSlaveZ: ", round(toq()))
    # return the prefined model, to be used in SolveSlaveZ
    return slaveZ, p, r, u_consumer, z_consumer, y_consumer, surplus

end



# change the objective function of slaveZ at each iterations
function SolveSlaveZ(u_consumer_initial, slaveZ, p, r,u_consumer, z_consumer, y_consumer, surplus,
    niterationsDD, dual_reliability, dual_profits)


    #  make u_consumer equal to the predefined u_consumer_initial
    @constraints(slaveZ, begin
        [consumer = ConsumersSet, option = OptionsSet],
            u_consumer[consumer, option] == u_consumer_initial[consumer, option]
    end)

    @constraints(slaveZ, begin
        # price and reliability in the menu should be increasing
        [option = 1:NOptions-1],
           r[option] <= r[option+1]
        [option = 1:NOptions-1],
            p[option] <= p[option+1]
        r[NOptions] <= 1
    end)

    for option = OptionsSet
        for consumer = ConsumersSet
            @constraints(slaveZ, begin
                surplus[consumer] >= Vl[consumer]*r[option]-p[option]
                # McCormick envelopes
                #=
                z_consumer[consumer, option]/zSCALE <= u_consumer[consumer, option]
                z_consumer[consumer, option]/zSCALE <= r[option]
                z_consumer[consumer, option]/zSCALE >= u_consumer[consumer, option] + r[option] - 1
                y_consumer[consumer, option]/ySCALE <= u_consumer[consumer, option]*PriceMax
                y_consumer[consumer, option]/ySCALE <= p[option]
                y_consumer[consumer, option]/ySCALE >= u_consumer[consumer, option]*PriceMax + p[option] - PriceMax
                =#
                z_consumer[consumer, option]/zSCALE == r[option]*u_consumer_initial[consumer, option]
                y_consumer[consumer, option]/ySCALE == p[option]*u_consumer_initial[consumer, option]
            end)
        end
    end

    ## multiple parts of the objective function
    # Company Revenue
    objExpression_1 = AffExpr()
    #=
    for option = OptionsSet, consumer = ConsumersSet
        append!(objExpression_1, y_consumer[consumer,option]/ySCALE*Dl[consumer])
    end
    =#
    for option = OptionsSet
        append!(objExpression_1, p[option]*Di[option])
    end
    objExpression_1 = objExpression_1*sum(SystemProfile)*dual_profits

    # subscription quantity
    objExpression_2 = AffExpr()
    #=
    for option = OptionsSet, consumer = ConsumersSet
        append!(objExpression_2, z_consumer[consumer,option]/zSCALE*Dl[consumer]*dual_reliability[option])
    end
    =#
    for option = OptionsSet
        append!(objExpression_2, r[option]*Di[option]*dual_reliability[option])
    end

    objExpression_2 = objExpression_2*sum(SystemProfile)
    @objective(slaveZ, Min, objExpression_1 + objExpression_2)
    tic()
    # solve it
    status_slaveZ = solve(slaveZ; relaxation=true)
    if status_slaveZ == :Optimal
        currentslaveZobjval = getobjectivevalue(slaveZ)

        println("Solving Slave_Z takes: ", round(toq(),2) )
        #println("SlaveZ: ", status_slaveZ,"; ObjValue: ", round(currentslaveZobjval))
        #println("***** consumer subproblem is solved")
        # get solutions
        r_sol = getvalue(r[:])
        p_sol = getvalue(p[:])
        u_consumer_sol = getvalue(u_consumer[:,:])
        y_consumer_sol = getvalue(y_consumer[:,:])
        z_consumer_sol = getvalue(z_consumer[:,:])

        #println("dual_reliability: ", dual_reliability)
        #println("r_sol: ", r_sol)
        #println("p_sol: ", p_sol)
        flush(STDOUT)


    end

    return  r_sol, p_sol, u_consumer_sol, y_consumer_sol, currentslaveZobjval

end
