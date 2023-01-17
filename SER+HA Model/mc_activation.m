%% Function of mast cell travel and activation in the brain.
%mc_switch is the presence or lack of neuroinflammation,
% mc_start is the time of start of neuroinflammation, t is time. 
% 
function f = mc_activation(t, mc_switch, mc_start_time)
c = 1; %Amplitude of MC activation. 
b = 20; %Strength of increase. 
t1 = mc_start_time + log(999)/b; %Start time is 0.01% of sigmoid function.
n=length(t);
for i=1:n
    if mc_switch == 0
        f(i) = 0;
    else
        if t < mc_start_time
            f(i) = 0.001;
        else
            f(i) = c/(1 + exp(-b * (t(i) - t1)));
    end
end
end


% function f = mc_activation(t, mc_switch, mc_start_time) %Linear
% b = 0.25; %Strength of increase. 
% n=length(t);
% max_f = 1;
% for i=1:n
%     if (mc_switch == 0) || (t(i) < mc_start_time)
%         f(i) = 0;
%     else
%         f(i) = b .* (t(i) - mc_start_time);
%     end
% end
% f(f > max_f) = max_f;
% end