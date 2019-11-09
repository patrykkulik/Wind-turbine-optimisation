%simpson rule
%INPUTS
%a - the index of the starting point
%b - the index of the final point
%points - an array of the function data (first colum is x, second is additional data etc.) 
%func - a function handle to the function to be simpton ruled! The function
%should be written to accept the points array
%
%OUTPUTS
%result - the result of the integration

function result = simpsonInt(a ,b ,points, func)
    result = sum( ((points(a+1:b,1) - points(a:b-1,1)) ./ 6) .*(func(points(a:b-1,:)) + 4 * func((points(a:b-1,:) + points(a+1:b,:))/2) + func(points(a+1:b,:))));
end
