%% Function of serotonin neuron firing.
%Commented functions are different firing paradigms. 
% UNITS of f() in events/h, time variables in seconds. 

function f = fireht(t, i_factor) %Repeated stim. 
n=length(t);
r = 8; 
repeat_time = 60*10;
t_start = 5;
t_flip = 7; 
t_end = 15;
basal = 1;
b = 2; 
max_f = 1.5; 
stim_boolean = 1;

if stim_boolean == 0
    for i=1:n
        f(i) = basal.*i_factor(i);
    end
    f(f>max_f) = max_f;
   return
end

for i=1:n
    time = t(i)*3600;
    n_stim = floor(time/repeat_time);
    if time <= (t_start + n_stim * repeat_time)
        f(i) = basal.*i_factor(i);
    elseif time > (t_start + n_stim * repeat_time) && time <= (t_flip + n_stim * repeat_time)

        f(i) = i_factor(i).*(basal + r.*(1 - exp(-b.*((time - n_stim * repeat_time)-t_start))));  
        
    elseif time < (t_end + n_stim * repeat_time)

        f(i) = i_factor(i).*(basal + r.*(exp(-b.*((time - n_stim * repeat_time)-t_flip)) - exp(-b.*(time - n_stim * repeat_time))));  

    else 
        f(i) = i_factor(i).*(basal);
    end
end
f(f>max_f) = max_f;
end








