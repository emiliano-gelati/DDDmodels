function SnowpackTemp(Temp)
#This function calculates a weighted (linearly decreasing) mean of the air temperature as an estimate of the skin tmeperature
# of the snowpack. It is assumed that the temperature at ground is zero such that the mean temperature is estT/2
# Temp: : 1 dim array float
# The temperature vector is organzed such that  the first element is todays temp
    len = length(Temp)
    wgt = collect(len:-1:1) / (len*(len+1)/2)
    snittT = sum(wgt.*Temp)          # This is the snowpack temperature. Tss is snittT *2
    if snittT > 0.0 
        snittT = 0.0 # Snowpack temperature cannot be more than zero
    end
return snittT
end
