function EI_p = EI_acq(initial_Dom,initial_Domy,new_x,theta,sigma,alpha,num_samp,inv_R,min_obj,EI_acq_mode,ratio_or_weight)  

    % EI evaulation
    [pred_mean,pred_sig] = GPR(initial_Dom,initial_Domy,new_x,theta,sigma,alpha,num_samp,inv_R);

    % Based on the GPR model, EI process is progressed.
    u = (min_obj - pred_mean);     
    ss = sqrt(pred_sig);             
        
    switch EI_acq_mode
        case 'normal'
            if ss > 0
                EI_p =  ((u-ratio_or_weight)*normcdf(u/ss) + ss*normpdf(u/ss));   % EI calculation process
            else
                EI_p = 0;
            end
        case 'weighted'
            % ref. On the Design of Optimization Strategies Based on...
            % Global Response Surface Approximation Models (2005)
            % Journal of Global Optimization
            if ss > 0
                EI_p =  (ratio_or_weight*u*normcdf(u/ss) + ss*(1-ratio_or_weight)*normpdf(u/ss));   % EI calculation process
            else
                EI_p = 0;
            end
        otherwise
    end
end