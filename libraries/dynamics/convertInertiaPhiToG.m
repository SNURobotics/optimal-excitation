%% Convert Single Inertia Matrix Phi to G
% 2018 Bryan Dongik Lee

%% Inputs
% [Name]  [Description]                                 [Size]
%  Phi     link inertia matrix for linear form dynamics  10*1

%% Outputs
% [Name]  [Description]                                 [Size]
%  G       link inertia matrix for standard dynamics     6*6

%% Implementation
function G = convertInertiaPhiToG(Phi)
    % error check
    if(size(Phi,1) ~= 10)
        error(['ERROR: Wrong input size for PhiToG(): Phi.rows() = ' num2str(size(Phi,1))])
        G = nan(6,6);
        return
    end

    % parsing
    G = zeros(6,6);
    m_Eye = Phi(1) * eye(3);
    h_bracket = skew(Phi(2:4));
    I_moment = [Phi(5), Phi(8), Phi(10);
                Phi(8), Phi(6), Phi(9);
                Phi(10), Phi(9), Phi(7)];

    % substitute
    G(1:3,1:3) = I_moment;
    G(1:3,4:6) = h_bracket;
    G(4:6,1:3) = h_bracket';
    G(4:6,4:6) = m_Eye;
end