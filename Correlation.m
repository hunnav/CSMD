function r_r = Correlation(x_sample,new_x,theta)

    % Gaussian Process Regression
    Ath = diag(theta);  % make diagonal matrix
    if size(new_x,2) > size(new_x,1)   % if something wrong, make correct
        new_x = new_x';
    end
    
    Mxy = new_x'*Ath*new_x;
    Mxx = diag(Mxy);
    
    Kxy = x_sample'*Ath*x_sample;
    Kxx = diag(Kxy);
    
    Nxy = new_x'*Ath*x_sample;
    
    corr = Mxx + Kxx - 2*Nxy';
    r_r = exp(-corr);          % make correlation matrix r (between unknown point x and known point)
end