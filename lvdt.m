function [di] = lvdt(x1,y1,x2,y2)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
di=sqrt((y2-y1).^2+(x2-x1).^2)
end

