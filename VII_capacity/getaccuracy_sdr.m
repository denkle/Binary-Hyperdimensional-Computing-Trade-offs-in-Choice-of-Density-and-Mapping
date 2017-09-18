function acc = getaccuracy_sdr (dp_hit,dp_rej,D)
%

%% Calculate accuracy of the encoding according to the equation for p_corr
%Use the heuristic which says that variance is quite similar to the mean
%value
 var_hit=1.1*dp_hit;
 var_rej=1.1*dp_rej;
 
 fun = @(u)  (1/(sqrt(2*pi*var_hit  ))).*exp(-((u-(dp_hit-dp_rej)).^2)/(2*(var_hit))).*((normcdf((u/sqrt(var_rej)))).^(D-1) );   
 acc= integral(fun,-5*sqrt(var_rej), dp_hit+5*sqrt(var_hit))   ;
 


end
