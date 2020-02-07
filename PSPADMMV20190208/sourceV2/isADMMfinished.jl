## chech convergence of ADMM
# if the number of interations have reached the max number of interations
# or both the primalResidual_rel and dualResidual_rel reach the target, we can
# stop ADMM
function isADMMfinished(niterationsADMM_input, primalResidual_rel_input, dualResidual_rel_input, objvalue_slaveZPR)
    if primalResidual_rel_input <= ResidualTarget && dualResidual_rel_input <= ResidualTarget && objvalue_slaveZPR <= 1e-6
        println("Residual targets reached")
        global primalResidual_rel = 100    # reset to 100%
        global dualResidual_rel = 100
        return true
    elseif niterationsADMM_input > MaxIterationADMM
        println("MaxIterationADMM reached")
        return true
    else
        return false
    end

end
