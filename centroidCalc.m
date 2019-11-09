function [xRoid, yRoid]  = centroidCalc(aeroFoilArea, aeroFoilPoints)

     contIntFuncHandle = @centIntFunc;
     xRoid = (1/aeroFoilArea) * -simpsonInt(1, length(aeroFoilPoints),aeroFoilPoints', contIntFuncHandle);
     yRoid = (1/aeroFoilArea) * simpsonInt(1, length(aeroFoilPoints),[aeroFoilPoints(2,:); aeroFoilPoints(1,:)]', contIntFuncHandle);
     
     function func = centIntFunc(points)
         func = points(:,1) .* points(:,2);
     end
end