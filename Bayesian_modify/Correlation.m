function r_r = Correlation(x_sample,new_x)
    nSample = size(x_sample,2);   % # of the Samples 
    dim = size(x_sample,1);       % dim of inputs    
    % Gaussian Process Regression
    r_r = zeros(nSample,1,dim);
    new_x = new_x';

    for i = 1:dim
        Mxx = new_x(i)^2;
        Kxx = (x_sample(i,:).^ 2);
        Nxy = new_x(i)'*x_sample(i,:);
        
        r_r(:,:,i) = Mxx + Kxx' - 2*Nxy';
    end

end