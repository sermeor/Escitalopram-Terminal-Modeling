%% Function of serotonin neuron firing.
%Commented functions are different firing paradigms. 
% UNITS IN events/h. 

function f=fireht(t) % For basal studies. 

f = 1;

end

% function f = fire(t) % NEW! SM 17.08.2022 FUNCTION DEFINED IN SECONDS AND ASSUMING t is hours.
% n=length(t);
% r = 1.5; %15;
% rest = 0;
% t_start = 5 + rest;
% t_flip = 6 + rest;
% basal = 1;
% b = 2;
% for i=1:n
%     time = t(i)*3600;
%     if time < t_start
%         f(i) = basal;
%     elseif time > t_start && time < t_flip
%         f(i) = basal + r.*(1 - exp(-b.*(time-t_start)));     
%     else f(i) = basal + r.*(exp(-b.*(time-t_flip)) - exp(-b.*time));
%     end
% end

% function f=fire(t); % 09/08/2022 SM : function of firing in hours. 
% n=length(t);
% basal = 1;
% max_amp = 15;
% t_start = 5/3600;
% t_flip = 6/3600;
% t_finish = 7/3600;
% slope_up = (max_amp-basal)/(t_flip - t_start); 
% slope_down = (basal-max_amp)/(t_finish - t_flip);
% for i=1:n
%     if t(i) < t_start
%         f(i) = basal;
%     elseif t(i) > t_start && t(i) < t_flip
%         f(i) = basal + slope_up * (t(i)-t_start); 
%     elseif t(i) > t_flip && t(i) < t_finish 
%         f(i) = basal + slope_up * (t_flip-t_start) + slope_down * (t(i) - t_flip);
%     else 
%         f(i) = basal;
%     end
% end
%  
% end







