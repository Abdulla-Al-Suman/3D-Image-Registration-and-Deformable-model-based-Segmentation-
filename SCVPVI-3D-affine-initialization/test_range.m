function [R0,I] = test_range(T11,T22)

global lev 
for lev = 1:3
    I0 = T11{lev};
    R0 = T22{lev};
[I] = register(I0,R0);
end

