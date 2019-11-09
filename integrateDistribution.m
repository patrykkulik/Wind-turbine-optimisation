function intResult = integrateDistribution(radiusPoints,inputDistribution, yFuncHandle)
    for i = 1:length(inputDistribution)
        intResult(i) = simpsonInt(1,i, [radiusPoints(:,1),inputDistribution], yFuncHandle);
    end
end