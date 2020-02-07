## function to read system data from CSV tables
#using DataFrames, CSV, Missings
function ReadSystemTables(systemdir::AbstractString; DisplayWarnings::Bool=false)

	# read static technical data using standard CSV.read function
	#Float64m = Union{Float64, Missings.Missing}
	# read technical specificatins of conventional generators
	Generators = readtable(string(systemdir, "/GeneratorsConv2050.csv"))
	HeatRateCurves = readtable(string(systemdir, "/HeatRateCurves.csv"))
	FuelPrice = readtable(string(systemdir, "/FuelPrices.csv"))

	#=
	Generators = CSV.read(string(systemdir, "/GeneratorsConv2050.csv"),
		types=[String, Float64, Float64, Float64, Int64,
		Int64, Float64, Float64, Float64, String, String, String, String],
		nullable=false)
	HeatRateCurves = CSV.read(string(systemdir, "/HeatRateCurves.csv"),
		types=[String, String, Float64, Float64], nullable=false)
	FuelPrice = CSV.read(string(systemdir, "/FuelPrices.csv"),
		types=[String, String, Float64], nullable=false)
		# initial values of consumers' choice
		=#
	if test_5D == 1
		# demand and valuation of 1000 consumers, i.e., V_l, D_l
		DlVl = readtable(string(systemdir, "/DlVl5D.csv"), header = false)
		SystemProfile = readtable(string(systemdir, "/SystemProfile5D.csv"), header = false)
		ConsumerProfile = readtable(string(systemdir, "/ConsumerProfile5D.csv"), header = false)
		ICcostHourly = readtable(string(systemdir, "/ICcostHourly5D.csv"), header = false)
		ICCost = readtable(string(systemdir, "/ICcost5D.csv"), header = false)
		u_consumer_sol = readtable(string(systemdir, "/u_consumer_sol5DTest.csv"), header = false)

#=
		DlVl = CSV.read(string(systemdir, "/DlVl5D.csv"), types=[Float64, Float64];
			header=[], datarow = 1, nullable=false)
		DiVi = CSV.read(string(systemdir, "/DiVi5D.csv"), types=[Float64, Float64];
			header=[], datarow = 1, nullable=false)
		SystemProfile =  CSV.read(string(systemdir, "/SystemProfile5D.csv"), types=[Float64];
			header=[], datarow = 1, nullable=false)
		ConsumerProfile = CSV.read(string(systemdir, "/ConsumerProfile5D.csv"), types=fill(Float64,Thours);
			header=[], datarow = 1, nullable=false)
		# read hourly IC cost
		ICcostHourly = CSV.read(string(systemdir, "/ICcostHourly5D.csv"), types=[Float64];
			header=[], datarow = 1, nullable=false)
		ICCost = CSV.read(string(systemdir, "/ICcost5D.csv"),types=[Float64];
			header=[], datarow = 1, nullable=false)
		u_consumer_sol = CSV.read(string(systemdir, "/u_consumer_sol5D.csv"), types=fill(Int,NOptions);
			header=[], datarow = 1, nullable=false)
=#
	else

		DlVl = readtable(string(systemdir, "/DlVl.csv"), header = false)
		SystemProfile = readtable(string(systemdir, "/SystemProfile.csv"), header = false)
		ConsumerProfile = readtable(string(systemdir, "/ConsumerProfile.csv"), header = false)
		ICcostHourly = readtable(string(systemdir, "/ICcostHourly.csv"), header = false)
		ICCost = readtable(string(systemdir, "/ICcost.csv"), header = false)

#=
		DlVl = CSV.read(string(systemdir, "/DlVl.csv"), types=[Float64, Float64];
			header=[], datarow = 1, nullable=false)

		DiVi = CSV.read(string(systemdir, "/DiVi.csv"), types=[Float64, Float64];
			header=[], datarow = 1, nullable=false)
		SystemProfile =  CSV.read(string(systemdir, "/SystemProfile.csv"), types=[Float64];
			header=[], datarow = 1, nullable=false)
		ConsumerProfile = CSV.read(string(systemdir, "/ConsumerProfile.csv"), types=fill(Float64,Thours);
			header=[], datarow = 1, nullable=false)
		# read hourly IC cost
		ICcostHourly = CSV.read(string(systemdir, "/ICcostHourly.csv"), types=[Float64];
			header=[], datarow = 1, nullable=false)
		ICCost = CSV.read(string(systemdir, "/ICcost.csv"),types=[Float64];
			header=[], datarow = 1, nullable=false)
=#
	u_consumer_sol = readtable(string(systemdir, "/u_consumer_solFullYearTest.csv"), header = false)


	end
	# ratio of available capacity of gas units
	GasAvailability =  readtable(string(systemdir, "/GasAvailableCapacityRatio2050Hourly.csv"))
	# Historicl profiles of TotalLoad, ICDemand Wind, Solar and Import
	HistoricalProfiles = readtable(string(systemdir, "/historicalProfile2050Hourly.csv"))
	# read renewable profiles
	WindProfiles = readtable(string(systemdir, "/WindRTHourlySce.csv"), header = false)
	SolarProfiles = readtable(string(systemdir, "/SolarRTHourlySce.csv"), header = false)

	#=
	# ratio of available capacity of gas units
	GasAvailability =  CSV.read(string(systemdir, "/GasAvailableCapacityRatio2050Hourly.csv"),
		types=[Float64], nullable=false)
	# Historicl profiles of TotalLoad, ICDemand Wind, Solar and Import
	HistoricalProfiles = CSV.read(string(systemdir, "/historicalProfile2050Hourly.csv"),
		types=[Float64, Float64, Float64, Float64, Float64], nullable=false)
	=#
	println("Tables are read and returned", now())
	# return all read data frames
	return Generators, HeatRateCurves, FuelPrice, DlVl,
		SystemProfile, ConsumerProfile,	GasAvailability,
		HistoricalProfiles, ICCost,ICcostHourly,
		WindProfiles, SolarProfiles, u_consumer_sol
end
