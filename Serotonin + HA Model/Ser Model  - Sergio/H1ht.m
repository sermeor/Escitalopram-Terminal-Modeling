
%% Factor of uptake 2 velocity (Vu2) depending on eht (a)
% UNITS IN uM. 

function f=H1ht(a)
    if a < 0.0605
          f = 0;
    elseif  a < 0.0805
         f = (50)*(a-0.0605);   
    else
        f = 1;
    end
end