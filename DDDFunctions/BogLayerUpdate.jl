#    Function Layer_update
#
#------------------------------------------------------------------------
#     Description:  Updates Layers by todays event and shifting to one timestepahead. 
#
#     Author: Thomas Skaugen
#     Revised: 16.12.2019
#--------------------------------------------------------------------------

function BogLayerUpdate!(outbog, BogLayers, UHBog, nodaysvector)
    qlayer = outbog * UHBog    #finds response  in mm!!!!  for bog, a vector     
    for t in 2:nodaysvector
        BogLayers[t-1] = BogLayers[t] + qlayer[t] # shifts the level of the matrix one timestep ahead
    end
    if nodaysvector == 1
        BogLayers .= 0.
#        BogLayers .= qlayer
    end
end            


