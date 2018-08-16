%% Convert Single Inertia Matrix S to G
% 2018 Bryan Dongik Lee

%% Inputs
% [Name]  [Description]                                       [Size]
%  S       link inertia matrix in positive semi-definite form  4*4

%% Outputs
% [Name]  [Description]                                       [Size]
%  G       link inertia matrix for standard dynamics           6*6

%% Implementation
function G = convertInertiaSToG(S)
    % Exception
    if(size(S,1) ~= [4,4])
        error('convert error: wrong input size')
        G = nan(6,6);
        return 
    end

    % // parsing
    m_Eye = S(4,4) * eye(3);
    h_bracket = skew(S(1:3,4));
    Sigma = S(1:3,1:3);
    I_moment = trace(Sigma) * eye(3) - Sigma;

    % // substitute
    G = [I_moment, h_bracket;
         h_bracket',  m_Eye];
end