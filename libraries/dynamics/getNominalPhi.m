%% Generate Random Nominal Value for Single Link Inertia Phi
% 2018 Bryan Dongik Lee

%% Inputs
% [Name]        [Description]                                       [Size]
%  Phi           link inertia matrix for linear form dynamics        10*1

%% Outputs
% [Name]        [Description]                                       [Size]
%  nominal_Phi   physically consistant random nominal value of Phi   10*1

%% Implementation
function nominal_Phi = getNominalPhi(Phi)
    S = convertInertiaGToS(convertInertiaPhiToG(Phi));
    m = S(4,4);
    sigma_bar = S(1:3,1:3)/m;
    p = S(1:3,4)/m;
    sigma_bar_c = sigma_bar - p*p';
    
    % variance parameter for mass
    variance = 0.3;
    nominal_m = m *(1 + variance*(rand(1)-0.5));
    
    % variance paramter for inertia     
    num_sample = 50;
    x = zeros(3, num_sample);
    for i=1:num_sample
        x(:,i) = p + sqrtm(sigma_bar_c) * randn(3,1);
    end   
    
    nominal_p = mean(x,2);
    nominal_sigma_bar_c = zeros(3,3);
    for i=1:num_sample
        nominal_sigma_bar_c = nominal_sigma_bar_c + (x(:,i)-nominal_p)*(x(:,i)-nominal_p)';
    end
    nominal_sigma_bar_c = nominal_sigma_bar_c/num_sample;
    
    nominal_S = nominal_m * [nominal_sigma_bar_c+nominal_p*nominal_p' nominal_p; nominal_p' 1];
    nominal_Phi = convertInertiaGToPhi(convertInertiaSToG(nominal_S));
    assert(min(eig(nominal_S))>0)
end