%% Histamine firing function. 
% Units in events/h.
% 
function f = fireha(t, i_factor) %Repeated stim. 
n=length(t);
r = 150;  
repeat_time = 60*10;
t_start = 5;
t_flip = 7;
t_end = 15;
basal = 1;
b = 2; 
max_f = 10; 
stim_boolean = 0;

if stim_boolean == 0
    for i=1:n
        f(i) = basal.*i_factor;
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













% 
% function f = fireha(t) % NEW! SM 17.08.2022 FUNCTION DEFINED IN SECONDS AND ASSUMING t is hours.
% n=length(t);
% r = 12.5; 
% rest = 0;
% t_start = 5 + rest;
% t_flip = 7 + rest;
% basal = (1)*1;
% b = 2; 
% for i=1:n
%     time = t(i)*3600;
%     if time < t_start
%         f(i) = basal;
%     elseif time > t_start && time < t_flip
%         f(i) = basal + r.*(1 - exp(-b.*(time-t_start)));     
%     else 
%         f(i) = basal + r.*(exp(-b.*(time-t_flip)) - exp(-b.*time));
%     end
% end
