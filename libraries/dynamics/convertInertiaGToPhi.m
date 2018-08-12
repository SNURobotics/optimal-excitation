%% Convert Single Inertia Matrix G to Phi
% 2018 Bryan Dongik Lee

%% Inputs
% [Name]  [Description]                                 [Size]
%  G       link inertia matrix for standard dynamics     6*6

%% Outputs
% [Name]  [Description]                                 [Size]
%  Phi     link inertia matrix for linear form dynamics  10*1

%% Implementation
function Phi = convertInertiaGToPhi(G)
    % // Exception
    if(size(G,1) ~= 6 || size(G,2) ~= 6)
        error('convert error: wrong input size for convertInertiaGToPhi().')
        Phi = nan(10,1);
        return
    end
    %     // parsing
    m = G(4,4);
    h = skew(G(1:3,4:6));
    I = [G(1,1); G(2,2); G(3,3); G(1,2); G(2,3); G(1,3)];

    %     // substitute
    Phi = [m; h; I];
end