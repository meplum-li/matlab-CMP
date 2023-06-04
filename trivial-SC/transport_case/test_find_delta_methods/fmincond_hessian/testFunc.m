function [outputArg1,outputArg2,out3] = testFunc(x,y,z)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

outputArg1 = x;
if nargout>1
outputArg2 = y;
end
 if nargout>2
     out3=z;
 end
end