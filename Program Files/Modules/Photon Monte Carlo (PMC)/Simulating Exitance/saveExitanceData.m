function temp = saveExitanceData(temp, eData, inPMC)
    if temp.iter ~= 1
        return
    elseif inPMC.saveExitance == "on"
        temp.solarExitance = eData.solar;
    end
end