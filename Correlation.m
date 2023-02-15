function r = Correlation(Domain,x)
    nSample = size(Domain,2);   % # of the Samples 
    dim = size(Domain,1);       % dim of inputs    
    % Gaussian Process Regression
    r = zeros(nSample,1,dim);
    x = x';

    for i = 1:dim
        Mxx = x(i)^2;
        Kxx = (Domain(i,:).^ 2);
        Nxy = x(i)'*Domain(i,:);
        
        r(:,:,i) = Mxx + Kxx' - 2*Nxy';
    end

end