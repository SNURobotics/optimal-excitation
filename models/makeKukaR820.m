function robot = makeKukaR820()
       
    % degrees of freedom
    robot.dof = 7;
    
    % screws A_i, i-th screw described in i-th frame
    robot.A = [0 0 0 0 0 0 0
               0 0 0 0 0 0 0
               1 1 1 1 1 1 1
               0 0 0 0 0 0 0
               0 0 0 0 0 0 0
               0 0 0 0 0 0 0];
     
    % link frames M_{i,i-1}
    robot.M(:,:,1) = [1  0  0  0
                      0  1  0  0
                      0  0  1  0
                      0  0  0  1];
     
    robot.M(:,:,2) = [1  0  0  0
                      0  0  -1 0
                      0  1  0  0
                      0  0  0  1];
     
    robot.M(:,:,3) = [1  0  0  0
                      0  0  1  0.42
                      0  -1 0  0
                      0  0  0  1];
     
    robot.M(:,:,4) = [1  0  0  0
                      0  0  1  0
                      0  -1 0  0
                      0  0  0  1]; 

    robot.M(:,:,5) = [1  0  0  0
                      0  0  -1 -0.4
                      0  1  0  0
                      0  0  0  1];

    robot.M(:,:,6) = [1  0  0  0
                      0  0  -1 0
                      0  1  0  0
                      0  0  0  1];
     
    robot.M(:,:,7) = [1  0  0  0
                      0  0  1  0
                      0  -1 0  0
                      0  0  0  1];     
    
    for i = 1:robot.dof
        robot.M(:,:,i) = inverse_SE3(robot.M(:,:,i));
    end
    
    % matrix B for base parameters, Phi_B = B * Phi
    robot.B = [];
     
    % inertia matrix Phi for linear dynamics tau = Y * Phi
    robot.Phi = [];
    
%     robot.G = zeros(6,6,robot.dof);
%     for i = 1:robot.dof
%         robot.G(:,:,i) = convertInertiaPhiToG(robot.Phi(10*(i-1)+1:10*i));
%     end
%     
%     robot.pd_metric_Phi = zeros(10*robot.dof, 10*robot.dof);
%     for i = 1:robot.dof
%         robot.pd_metric_Phi(10*(i-1)+1:10*i, 10*(i-1)+1:10*i) = getPDMetricInertiaPhi(robot.Phi(10*(i-1)+1:10*i));
%     end
    
    % temporary..
%     robot.B_metric_inv_Phi_Bt = robot.B * pinv(robot.pd_metric_Phi) * robot.B';
%     robot.B_metric_inv_Phi_Bt = (robot.B_metric_inv_Phi_Bt+robot.B_metric_inv_Phi_Bt')/2;

    
%     robot.nominal_Phi = zeros(10*robot.dof,1);
%     for i=1:robot.dof
%         robot.nominal_Phi(10*(i-1)+1:10*i) = getNominalPhi(robot.Phi(10*(i-1)+1:10*i));
%     end
end