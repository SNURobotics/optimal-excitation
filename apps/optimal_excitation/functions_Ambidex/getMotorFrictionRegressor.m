function regressor = getMotorFrictionRegressor(motorvel)
    regressor = zeros(7, 14);
    
    cutoff = 4;
    
    for i = 1:7
        if motorvel(i) > cutoff
            regressor(i,i) = 1;
        elseif motorvel(i) > -cutoff
            regressor(i,i) = motorvel(i)/4;
        else
            regressor(i,i) = -1;
        end
        
        regressor(i,i+7) = motorvel(i);
    end
   
end