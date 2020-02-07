# Benders is not used because  the first iteration only prints some info
tic()   # record time to load moduels and data
HOST = "LENOVO"   #LENOVO, LEMAITRE

if HOST == "LEMAITRE"
    println("Date/time: ", now(), "\n\n")
    println("N.O. of proceses: ", nprocs())
    flush(STDOUT)
    @everywhere begin
        global  test_5D = 0
        global  SubHorizonChoice = "MONTH"  # "DAY", "WEEK", "MONTH"
        global  MIPGapSlave = 0.001
        global  MIPGapMaster = 0.001
        global  MaxIterationADMM = 30   # max iterations of ADMM in the X update
        global  ResidualTarget = 0.01   # % in percentage
        println(gethostname(), " is started")
        using JuMP, DataFrames, Gurobi
        NTHREADS = 1
        global  SysDir = "/home/users/y/m/ymou/PSPADMMV20190208"
        SimDir = string(SysDir,"/sourceV2")
        cd(SimDir)
        include("./Constants.jl")
        include("./ReadSystemTables.jl")
        include("./MISC.jl")
        include("./ADMMOutput.jl")
        SlaveXSolverChoice = GurobiSolver(
            OutputFlag=0,
            Threads=NTHREADS
            )
        SlaveXMIPSolverChoice = GurobiSolver(
            OutputFlag=0, MIPGap=MIPGapSlave, Threads=NTHREADS,
                TimeLimit=36000)
    end

elseif HOST == "LENOVO"
    println("Date/time: ", now(), "\n\n")
    addprocs(6- nprocs())
    println("N.O. of proceses: ", nprocs())
    flush(STDOUT)
    @everywhere begin
        global  test_5D = 1
        global  SubHorizonChoice = "DAY"  # "DAY", "WEEK", "MONTH"
        global  MIPGapSlave = 0.001
        global  MIPGapMaster = 0.001
        global  MaxIterationADMM = 30  # max iterations of ADMM in the X update
        global  ResidualTarget = 0.1   # % in percentage
        println(gethostname(), " is started")
        using JuMP, DataFrames, Gurobi
        NTHREADS = 1
        global  SysDir = "C:/Users/Yuting/OneDrive - UCL/Demand Response Business Model/Julia/PSPADMMV20190208"
        SimDir = string(SysDir,"/sourceV2")
        cd(SimDir)
        include("./Constants.jl")
        include("./ReadSystemTables.jl")
        include("./time_cycling.jl")
        include("./ADMMOutput.jl")
        SlaveXSolverChoice = GurobiSolver(
            OutputFlag=0,
            Threads=NTHREADS)
        SlaveXMIPSolverChoice = GurobiSolver(
            OutputFlag=0, MIPGap=MIPGapSlave, Threads=NTHREADS,
                TimeLimit=36000)
    end

else
    println("Host is wrongly defined")
end

@eval @everywhere begin
    GasHourlyDerating = 1
    Generators, HeatRateCurves, FuelPrice, DlVl,
        SystemProfile, ConsumerProfile,	GasAvailability,
        HistoricalProfiles, ICCost,ICcostHourly,
        WindProfiles, SolarProfiles, u_consumer_initial = ReadSystemTables(InputDirectory)
    u_consumer_initial = convert(Array, u_consumer_initial)
    DlVl = convert(Array, DlVl)
    Dl = DlVl[:, 1]
    Vl = DlVl[:, 2]
    Vmax = 487
    Vmin = 0
    vb_LB = [Vmin + (Vmax - Vmin)*(n-1)/NPartitions for n = 1:NPartitions]    # lower bound of vb
    vb_UB = [Vmin + (Vmax - Vmin)*n/NPartitions for n = 1:NPartitions] # upper bound of vb
    #vb_LB = vb_UB   # set lower bound and upper bound to be the same
    VB = Array{Float64}(NOptions+1)   # valuation breakpoints parameters
    VB[1] = vb_LB[1]    # the first consumer is the first breakpoint
    for option = 2:NOptions+1
        VB[option] = vb_UB[sum(u_consumer_initial[:,1:option-1])]
    end
    Di = Array{Float64}(NOptions)
    for option = OptionsSet
        Di[option] = sum(u_consumer_initial[consumer,option]*Dl[consumer] for consumer = ConsumersSet)
    end
    ConsumerProfile = convert(Array, ConsumerProfile)
    SystemProfile = convert(Array, SystemProfile)
    # process data
    MinRunCapacityDerated = Array{Float64}(NGenerators, Thours)
    MaxRunCapacityDerated = Array{Float64}(NGenerators, Thours)
    MinRunCapacityDerated .= Generators[:MinRunCapacity]
    MaxRunCapacityDerated .= Generators[:MaxRunCapacity]
    # gas units are to be derated
    GasGeneratorsSet = view(GeneratorsSet, Generators[:FuelGenerator].=="NAT_GAS")
    if GasHourlyDerating == 1
        for t = TimeSet
            MinRunCapacityDerated[GasGeneratorsSet,t] =
                MinRunCapacityDerated[GasGeneratorsSet,t]*GasAvailability[:Ratio][t]*GasRatioScaleUp
            MaxRunCapacityDerated[GasGeneratorsSet,t] =
                MaxRunCapacityDerated[GasGeneratorsSet,t]*GasAvailability[:Ratio][t]*GasRatioScaleUp
        end
    else
        GasAvailabilityAverage = mean(GasAvailability[:Ratio][TimeSet])
        for t = TimeSet
            MinRunCapacityDerated[GasGeneratorsSet,t] =
                MinRunCapacityDerated[GasGeneratorsSet,t]*GasAvailabilityAverage*GasRatioScaleUp
            MaxRunCapacityDerated[GasGeneratorsSet,t] =
                MaxRunCapacityDerated[GasGeneratorsSet,t]*GasAvailabilityAverage*GasRatioScaleUp
        end
    end
    ICcostTotal = convert(Array, ICCost)
    ICcostTotal = sum(ICcostTotal[:,scenario]*ScenarioProb[scenario] for scenario = ScenariosSet)
    ICcostTotal = - ICcostTotal[1,1]
    ##################################
    ############ the constant term in Dual decompostion
    Profits = Profits - ICcostTotal # ic cost must be shifted
    #############################################################

    ICcostHourly = convert(Array, ICcostHourly)
    WindProfiles = convert(Array, WindProfiles)
    SolarProfiles = convert(Array, SolarProfiles)

    # in the 5d version, scale the production level of renewable
    if test_5D == 1
        WindEnergyRef = sum(HistoricalProfiles[:Wind][TimeSet])
        SolarEnergyRef = sum(HistoricalProfiles[:Solar][TimeSet])
        WindEnergySce = sum(WindProfiles[TimeSet,:],1)
        SolarEnergySce = sum(SolarProfiles[TimeSet,:],1)
        WindRatio = WindEnergySce/WindEnergyRef
        SolarRatio = SolarEnergySce/SolarEnergyRef
        for scenario = ScenariosSet
            WindProfiles[:,scenario] = WindProfiles[:,scenario]/WindRatio[scenario]
            SolarProfiles[:,scenario] = SolarProfiles[:,scenario]/SolarRatio[scenario]
        end
    end
    NSubproblemsX = length(DaysX)   # number of subproblems in the X updates in each scenario
    NSubproblemsX = NSubproblemsX*NScenarios
end

# Slave X is solved in each process
@everywhere begin
    include("./SlaveXScaled.jl")
end

include("./isADMMfinished.jl")
include("./SlaveZScaled.jl")
include("./PrimalRecovery.jl")

println("Loading moduels and data takes: ", round(toq(),2))
flush(STDOUT)

# creat slaveZ model using function DefineSlaveZ()
slaveZ, p, r, z_consumer, y_consumer, surplus, supplyR_Z, totalCost_Z = DefineSlaveZ(u_consumer_initial)
tic()

# could use different values of penalty for different variables

rho = zeros(Float64, NOptions+1)
rho_vec = zeros(Float64, NOptions+1, MaxIterationADMM+1)  # rho is updated in each iteration
rho_vec[:,1] = RHO

############################################################################
println("###########################################")
println("*** Solving  Using ADMM ***")
flush(STDOUT)
println("*** Z and U are initialzied as 0 ***")
# define solutions from Z updates
global totalCost_Z_sol = zeros(Float64, Thours, NScenarios)
global supplyR_Z_sol = zeros(Float64, NOptions, Thours, NScenarios)
# define U, to be used in the Z update and U update
global totalCost_U = zeros(Float64, Thours, NScenarios)
global supplyR_U = zeros(Float64, NOptions, Thours, NScenarios)

global dualResidual = 1
global dualResidual_rel = 100 # give a initial value as 100%
global primalResidual = 1
global primalResidual_rel = 100 # give a initial value as 100%
# collect the relative residuals in each ADMM iterations
global dualResidual_rel_vec = 100*ones(Float64, MaxIterationADMM)
global primalResidual_rel_vec = 100*ones(Float64, MaxIterationADMM)
global objval_slaveX_sol_vec = zeros(Float64, MaxIterationADMM)
global r_sol_vec = zeros(Float64, NOptions, MaxIterationADMM)
global p_sol_vec = zeros(Float64, NOptions, MaxIterationADMM)
global objvalue_slaveZPR = 1
global objvalue_slaveZPR_vec = zeros(Float64, MaxIterationADMM)

niterationsADMM = 1
println("ADMM starts at ", now())
ADMMDateTime = Dates.format(now(), "yyyy-mm-dd-HH-MM-SS")
flush(STDOUT)

# while ADMM is not finished, keep running it
while ~isADMMfinished(niterationsADMM, primalResidual_rel, dualResidual_rel, objvalue_slaveZPR)

    println("***** ADMM iterations: ", niterationsADMM)
    flush(STDOUT)
    if niterationsADMM > 1
        rho = RHO
    end
    ################################################################################
    ################################################################################
    println("***** X-Update ")
    #println("***** Preparing input for ADMM X")
    flush(STDOUT)
    inputSlaveXVec = []
    for subproblemIdx = 1:NSubproblemsX
        # calculate the month index and scenario index based on
        MonthIdx, ScenarioIdx= returnIDX(subproblemIdx,NSubhoriozns)
        totalCost_Z_sol_current = totalCost_Z_sol[StartHourX[MonthIdx]:EndHourX[MonthIdx], ScenarioIdx]
        supplyR_Z_sol_current = supplyR_Z_sol[:, StartHourX[MonthIdx]:EndHourX[MonthIdx], ScenarioIdx]
        totalCost_U_current = totalCost_U[StartHourX[MonthIdx]:EndHourX[MonthIdx], ScenarioIdx]
        supplyR_U_current = supplyR_U[:, StartHourX[MonthIdx]:EndHourX[MonthIdx], ScenarioIdx]
        push!(inputSlaveXVec, [subproblemIdx, MonthIdx, ScenarioIdx,
        u_consumer_initial, totalCost_Z_sol_current, supplyR_Z_sol_current,
        totalCost_U_current,  supplyR_U_current, niterationsADMM, rho])
    end

    #println("***** running ADMM X Update in Parallel")
    SlaveX_solutions = pmap(SolveSlaveX, inputSlaveXVec)
    # define the solution of X, collecting them together
    global totalCost_X_sol =  zeros(Float64, Thours, NScenarios)#Array{Float64}(Thours)
    global supplyR_X_sol = zeros(Float64, NOptions, Thours, NScenarios) #Array{Float64}(NOptions, Thours)
    global objval_slaveX_sol = 0.0
    # solution from the X updates, (could be in shared array)
    # after convergence of ADMM, we take the average
    # recover the solution form X SubProblems
    for subproblemIdx = 1:NSubproblemsX
        MonthIdx, ScenarioIdx= returnIDX(subproblemIdx,NSubhoriozns)

        totalCost_X_sol[StartHourX[MonthIdx]:EndHourX[MonthIdx],ScenarioIdx] =
            SlaveX_solutions[subproblemIdx][2]
        supplyR_X_sol[:,StartHourX[MonthIdx]:EndHourX[MonthIdx],ScenarioIdx] =
            SlaveX_solutions[subproblemIdx][3]
        objval_slaveX_sol = objval_slaveX_sol + ScenarioProb[ScenarioIdx]*SlaveX_solutions[subproblemIdx][4]
    end
    objval_slaveX_sol_vec[niterationsADMM] = objval_slaveX_sol
    ################################################################################
    ################################################################################
    println("***** Z-Update")
    rho = RHO
    # call SolveSlaveZ
        totalCost_Z_sol, supplyR_Z_sol, y_consumer_sol, primalResidual_rel, dualResidual =
    SolveSlaveZ(slaveZ, p, r, z_consumer, y_consumer, surplus, supplyR_Z, totalCost_Z,
        totalCost_X_sol, supplyR_X_sol, totalCost_Z_sol, supplyR_Z_sol, niterationsADMM)

    ################################################################################
    ################################################################################
    println("*****  U Update")
    include("./SolveSlaveU.jl")

    # store the residuals
    primalResidual_rel_vec[niterationsADMM] = primalResidual_rel
    dualResidual_rel_vec[niterationsADMM] = dualResidual_rel

    println("ADMM iteration until ",niterationsADMM, " runs to ", now())

    ################################################################################
    ################################################################################
    println("*****  primal recovery")

    objvalue_slaveZPR, r_sol, p_sol = SlaveZPR(u_consumer_initial, totalCost_X_sol, supplyR_X_sol)
    objvalue_slaveZPR_vec[niterationsADMM] = objvalue_slaveZPR
    r_sol_vec[:, niterationsADMM] = r_sol
    p_sol_vec[:, niterationsADMM] = p_sol
    println("primalResidual_rel_vec %: ", primalResidual_rel_vec)
    println("dualResidual_rel_vec %: ", dualResidual_rel_vec)
    println("objval_slaveX_sol_vec: ", objval_slaveX_sol_vec)
    println("r_sol_vec: ", r_sol_vec')
    println("p_sol_vec: ", p_sol_vec')
    println("objvalue_slaveZPR: ", objvalue_slaveZPR)
    println("r_sol: ", r_sol)
    println("p_sol: ", p_sol)
    flush(STDOUT)

    # output the results from this ADMM iterations
    if objvalue_slaveZPR <= 1e-5
        println("Output data in this ADMM iteration")
        flush(STDOUT)
        ADMMOutput(NSubproblemsX, niterationsADMM, SlaveX_solutions, r_sol, p_sol,
            primalResidual_rel_vec, dualResidual_rel_vec, supplyR_U,  ADMMDateTime)
    end


    # output to files of the results from every iterations
    if test_5D == 1
        OutputDir = string(SysDir,"/output/ADMM5D", ADMMDateTime)
    else
        OutputDir = string(SysDir,"/output/ADMMFullYear", ADMMDateTime )
    end
    mkpath(OutputDir)

    ## print these related to  the bounds
    BoundsIdx = find(objvalue_slaveZPR_vec .<= 1e-5)
    objval_slaveX_sol_vec_bounds = objval_slaveX_sol_vec[BoundsIdx]
    primalResidual_rel_vec_bounds = primalResidual_rel_vec[BoundsIdx]
    dualResidual_rel_vec_bounds = dualResidual_rel_vec[BoundsIdx]
    r_sol_vec_bounds = r_sol_vec[:,BoundsIdx]
    p_sol_vec_bounds = p_sol_vec[:,BoundsIdx]

    BoundsDF = DataFrame(iteration = BoundsIdx,
                        UB =objval_slaveX_sol_vec_bounds,
                        primalResidual_rel = primalResidual_rel_vec_bounds,
                        dualResidual_rel = dualResidual_rel_vec_bounds
                        )
    writetable(string(OutputDir,"/Bounds.csv"),  BoundsDF; header=true)
    PriceMenuDF =  DataFrame([r_sol_vec_bounds; p_sol_vec_bounds])
    writetable(string(OutputDir,"/PriceMenu.csv"),  PriceMenuDF; header=true)

    PriceMenuAllDF =  DataFrame([r_sol_vec; p_sol_vec])
    writetable(string(OutputDir,"/PriceMenuAll.csv"),  PriceMenuAllDF; header=true)

    # update iteration number of ADMM
    niterationsADMM = niterationsADMM + 1
end

println("objvalue_slaveZPR_vec", objvalue_slaveZPR_vec)

println("It", " takes ", round(toq()/3600,2), " hours in total")
flush(STDOUT)




#println(objval_slaveX_sol_vec[BoundsIdx])
#println(primalResidual_rel_vec[BoundsIdx])
#println(dualResidual_rel_vec[BoundsIdx])



println("*** ALL Done ! ***")
println("Date/time: ", now(), "\n\n")
flush(STDOUT)
