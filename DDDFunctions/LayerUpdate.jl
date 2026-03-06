#    Function Layer_update
#
#------------------------------------------------------------------------
#     Description:  Updates Layers by todays event and shifting to one timestepahead. 
#
#     Author: Thomas Skaugen
#     Revised: 17.12.2019
#--------------------------------------------------------------------------

function LayerUpdate!(ddist::Vector{Float64}, outx::Float64, Layers::Matrix{Float64}, layerUH::Matrix{Float64}, nodaysvector::Vector{Int}, NoL::Int)
  for j in 1:NoL
    multiplier = ddist[j] * outx
    Layers[j,1] = multiplier * layerUH[j,1]
    for h in 2:nodaysvector[j]
      Layers[j,h-1] = Layers[j,h] + multiplier * layerUH[j,h]
    end
  end
end                
