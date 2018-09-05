function robot = makeAtlasV5()

    root = parseXML('atlas_v5_core.xml');

    stack = CStack();
    stack.push(root.Children(6));

    num_link = 1;
    robot_tree{num_link, 1}.name = 'root';
    robot_tree{num_link, 1}.children = 2;
    robot_tree{num_link, 1}.name = 'root';
    
    
    islink = false;

    while(~stack.isempty)

        current_node = stack.pop;

        for i = 1:size(current_node.Children,2)
            if strcmp(current_node.Children(i).Name, 'joint')
                islink = true;
                num_link = num_link+1;
%                 disp([num2str(num_link) ' links found']);

                robot_tree{num_link,1}.name = 'none';
                robot_tree{num_link,1}.pos = [0 0 0];
                robot_tree{num_link,1}.quat = [1 0 0 0];

                robot_tree{num_link,1}.joint_q = [0 0 0];
                robot_tree{num_link,1}.joint_w = [0 0 0];
                robot_tree{num_link,1}.joint_range = [0 0];

                robot_tree{num_link,1}.inertial_pos = [0 0 0];
                robot_tree{num_link,1}.inertial_quat = [1 0 0 0];
                robot_tree{num_link,1}.inertial_mass = 0;
                robot_tree{num_link,1}.inertial_diag = [0 0 0];

                robot_tree{num_link,1}.children_names = ["itself"];
                robot_tree{num_link,1}.children = [];

                % link frame
                for j = 1:size(current_node.Attributes,2)
                    if strcmp(current_node.Attributes(j).Name, 'name')
                        robot_tree{num_link,1}.name = current_node.Attributes(j).Value;
                    elseif strcmp(current_node.Attributes(j).Name, 'pos')
                        robot_tree{num_link,1}.pos = str2num(current_node.Attributes(j).Value);
                    elseif strcmp(current_node.Attributes(j).Name, 'quat')
                        robot_tree{num_link,1}.quat = str2num(current_node.Attributes(j).Value);
                    end
                end

                % joint
                for j = 1:size(current_node.Children(i).Attributes,2)
                    if strcmp(current_node.Children(i).Attributes(j).Name, 'pos')
                        robot_tree{num_link,1}.joint_q = str2num(current_node.Children(i).Attributes(j).Value);
                    elseif strcmp(current_node.Children(i).Attributes(j).Name, 'axis')
                        robot_tree{num_link,1}.joint_w = str2num(current_node.Children(i).Attributes(j).Value);
                    elseif strcmp(current_node.Children(i).Attributes(j).Name, 'range')
                        robot_tree{num_link,1}.joint_range = str2num(current_node.Children(i).Attributes(j).Value);
                    end
                end        
            end
        end

        % inertia
        if islink
            for i = 1:size(current_node.Children,2)
                if strcmp(current_node.Children(i).Name, 'inertial')
                    for j = 1:size(current_node.Children(i).Attributes,2)
                        if strcmp(current_node.Children(i).Attributes(j).Name, 'pos')
                            robot_tree{num_link,1}.inertial_pos = str2num(current_node.Children(i).Attributes(j).Value);
                        elseif strcmp(current_node.Children(i).Attributes(j).Name, 'quat')
                            robot_tree{num_link,1}.inertial_quat = str2num(current_node.Children(i).Attributes(j).Value);
                        elseif strcmp(current_node.Children(i).Attributes(j).Name, 'mass')
                            robot_tree{num_link,1}.inertial_mass = str2num(current_node.Children(i).Attributes(j).Value);
                        elseif strcmp(current_node.Children(i).Attributes(j).Name, 'diaginertia')
                            robot_tree{num_link,1}.inertial_diag = str2num(current_node.Children(i).Attributes(j).Value);
                        end
                    end 
                end
            end
        end

        for i = 1:size(current_node.Children,2)
            if strcmp(current_node.Children(i).Name, 'body') && size(current_node.Children(i).Children, 1) ~= 0
                stack.push(current_node.Children(i));

                if islink
                    for j = 1:size(current_node.Children(i).Children,2)
                        if strcmp(current_node.Children(i).Children(j).Name, 'joint')
                            robot_tree{num_link,1}.children_names = [robot_tree{num_link,1}.children_names;
                                                                     current_node.Children(i).Attributes(1).Value()];
                        end
                    end                   
                end
            end
        end

        islink = false;
    end

    % set children nodes
    for i = 2:size(robot_tree)
        for j = 2:size(robot_tree{i}.children_names)
            for k = 2:size(robot_tree)
                if strcmp(robot_tree{i}.children_names(j), robot_tree{k}.name)
                    robot_tree{i}.children = [robot_tree{i}.children; k];
                end
            end
        end
    end

    % set parent nodes
    for i = 2:size(robot_tree)
        for j = 1:size(robot_tree{i}.children)
            robot_tree{robot_tree{i}.children(j)}.parent = i;
        end
    end

    % set robot properties
    for i = 2:size(robot_tree)
        robot_tree{i}.A = [robot_tree{i}.joint_w'; cross(robot_tree{i}.joint_q, robot_tree{i}.joint_w)'];

        robot_tree{i}.M = eye(4);
        robot_tree{i}.M(1:3,1:3) = quat2rotm(robot_tree{i}.quat);
        robot_tree{i}.M(1:3,4) = robot_tree{i}.pos';
        robot_tree{i}.M = inverse_SE3(robot_tree{i}.M);
        
        robot_tree{i}.G_i = zeros(6,6);
        robot_tree{i}.G_i(1,1) = robot_tree{i}.inertial_diag(1);
        robot_tree{i}.G_i(2,2) = robot_tree{i}.inertial_diag(2);
        robot_tree{i}.G_i(3,3) = robot_tree{i}.inertial_diag(3);
        robot_tree{i}.G_i(4,4) = robot_tree{i}.inertial_mass;
        robot_tree{i}.G_i(5,5) = robot_tree{i}.inertial_mass;
        robot_tree{i}.G_i(6,6) = robot_tree{i}.inertial_mass;
        
        T = eye(4);
        T(1:3,1:3) = quat2rotm(robot_tree{i}.inertial_quat);
        T(1:3,4) = robot_tree{i}.inertial_pos';
        T = inverse_SE3(T);
        Adjoint = large_Ad(T);
        
        robot_tree{i}.G = Adjoint'*robot_tree{i}.G_i*Adjoint;
    end
    
    robot_tree{2}.parent = 1;
    robot_tree{1}.A = zeros(6,1);
    robot_tree{1}.M = eye(4);
    
    robot.dof = size(robot_tree,1) - 1;
    robot.A = zeros(6,robot.dof);
    robot.M = zeros(4,4,robot.dof);
    robot.G = zeros(6,6,robot.dof);
        
    for i = 1:robot.dof
        robot.A(:,i) = robot_tree{i+1}.A;
        robot.M(:,:,i) = robot_tree{i+1}.M;
        robot.G(:,:,i) = robot_tree{i+1}.G;
        robot.name{i} = robot_tree{i+1}.name;
        robot.tree{i}.parent = robot_tree{i+1}.parent - 1;
        robot.tree{i}.children = robot_tree{i+1}.children - 1;
    end        
end