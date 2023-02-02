function [EI_p] = EI_acq(initial_Dom,initial_Domy,new_x,theta,sigma,alpha,num_samp,inv_R,min_obj)  

    % EI evaulation
    [pred_mean,pred_sig] = GPR(initial_Dom,initial_Domy,new_x,theta,sigma,alpha,num_samp,inv_R);

    % Based on the GPR model, EI process is progressed.
    u = (min_obj - pred_mean);     
    ss = sqrt(pred_sig);             
        
    if ss > 0
       EI_p =  (u*normcdf(u/ss) + ss*normpdf(u/ss));   % EI calculation process
    else
       EI_p = 0;
    end

end