%% Function of transport from vha reserve to readily vha. 
% vha is the current value of vesicular histamine
% vha_basal is the steady state value. 
function ratio = vha_trafficking(vha, vha_basal)
max_r = 100; 
min_r = 0;
s = 15; %Strength
diff = vha_basal - vha; %Difference respect to basal. 
ratio = s*diff;
ratio(ratio>max_r) = max_r; 
ratio(ratio<min_r) = min_r; 
end