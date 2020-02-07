# Update U

totalCost_U = totalCost_U + totalCost_X_sol - totalCost_Z_sol
supplyR_U = supplyR_U + supplyR_X_sol - supplyR_Z_sol

dualDtemp = sum(totalCost_U.^2) +  sum(supplyR_U.^2)
dualResidual_rel = dualResidual/(dualDtemp^0.5)*100 # in percentage

#println("**** dualResidual_rel: ", round(dualResidual_rel,4), "%")
#flush(STDOUT)
