%% Convert Single Inertia Matrix G to S
% 2018 Bryan Dongik Lee

%% Inputs
% [Name]  [Description]                                      [Size]
%  G       link inertia matrix for standard dynamics          6*6

%% Outputs
% [Name]  [Description]                                      [Size]
%  S      link inertia matrix in positive semi-definite form  4*4

%% Implementation
function S = convertInertiaGToS(G)
    % Exception
    if(size(G,1) ~= [6,6])
        error('convert error: wrong input size')
        S = nan(4,4);
        return 
    end

    % parsing
    m      = G(4,4);
    h      = skew(G(1:3,4:6));
    I      = G(1:3,1:3);
    Sigma  = 0.5 * trace(I) * eye(3) - I;

    % substitute
    S =[Sigma,  h;
        h',     m];
end