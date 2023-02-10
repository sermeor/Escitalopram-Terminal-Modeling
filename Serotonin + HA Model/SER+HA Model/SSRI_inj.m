%% Function that simulates the injection of escitalopram into the peritoneum.
%t_start is the time for the first dose, and t_repeat is the repetition
%time in hours. inj_time is time that the injection lasts. Output units in ug/h. 
function f = SSRI_inj(t, t_start, t_repeat, q)
inj_time = 1/3600;
n=length(t);
for i=1:n
    if t(i)>(t_start)
        n_stim = floor((t(i)-t_start)/t_repeat);
        if (t(i)>(t_start + n_stim * t_repeat)) && (t(i) <= (t_start + inj_time + n_stim * t_repeat))
            f(i) =  q/inj_time;
        else
            f(i) = 0;
        end
    else 
        f(i) = 0;
    end
end
end