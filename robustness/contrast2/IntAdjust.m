function [ A ] = IntAdjust( A )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[m,n] = size(A);

for i = 1 : m
    for j = 1 : n
        if A(i,j) < -0.5
            A(i,j) = -1;
        else if A(i,j) < 0.5
                A(i,j) = 0;
            else
                A(i,j) = 1;
            end
        end
    end
end

end

