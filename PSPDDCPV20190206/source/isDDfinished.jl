## chech convergence of DD
# if the number of interations have reached the max number of interations
# or the gap have reached threshold
# stop DD
function isDDfinished(niterationsDD_input, gap_input)
    if abs(gap_input) <= GapTarget
        println("gap targets reached")
        global gap = 100    # reset to 100%
        return true
    elseif niterationsDD_input > MaxIterationDD
        println("MaxIterationDD reached")
        return true
    else
        return false
    end

end
