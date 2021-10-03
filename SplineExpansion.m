%% Spline Expansion function
%%% returns the values of a spline of order n at given times t.
function x = SplineExpansion(d,t,n) 
    N = length(d);% number of coefficients
    L = length(t);

    % calculate beta of order zero.
    beta_0 = zeros(1,L);
    for l = 1:L
        if abs(t(l))<0.5
            beta_0(l) = 1;  
        elseif abs(t(l))==0.5
            beta_0(l) = 0.5;
        end
    end
    
    % calculate beta of order n by convoluting n times.
    beta_n = beta_0;
    for i = 1:n              
        beta_n = conv(beta_0,beta_n,'same'); %keep the same size as beta_0
        beta_n = beta_n * (t(2)-t(1));% normalize amplitude of beta_n
    end
    
    % calculating the spline expansion
    x = zeros(1,L);
    for k = 1:N
        %delta between time samples ((t(end)-t(1))/length(t));
        step = t(2)-t(1);
        %number of time samples in a time interval of length 1
        unit_step = ceil(1/step); 
        %shift beta_n by k steps
        beta_n_shift_k = [zeros(1,unit_step*k) beta_n(1:end-unit_step*k)];
        %sum the shifted B-splines according to the splie expansion formula
        x = x + d(k)* beta_n_shift_k;
    end
    
   end
    

