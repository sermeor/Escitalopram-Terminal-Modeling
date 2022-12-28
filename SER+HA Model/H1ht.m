%% Factor of uptake 2 velocity (Vu2) depending on eht (a)
%and basal eht (c)
% UNITS IN uM. 
function f=H1ht(a, c)
    if a < c
          f = 0;
    elseif  a < (c + 0.02)
         f = (50)*(a-c);   
    else
        f = 1;
    end
end