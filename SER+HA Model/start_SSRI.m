%% Function to swtich on PK constant at the start of the dose.
% t is the time array or number and start_time is the time at which the dose is
% injected. f is a boolean (1 is active, 0 is inactive). 
function f = start_SSRI(t, start_time)
n=length(t);
for i=1:n
    if t(i)<start_time
        f(i) = 0;
    else
        f(i) = 1;
    end
end
end