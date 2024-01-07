%                                                           ^             %
%                                                          /|\            %
%                                                           |             %
%      HI JACK/DYLAN. PRESS RUN -----------------------------             %
%                                                                         %  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







% This script is used to test the GenerateLattice() funciton for AERO3261
% A2
% Author: Thomas Ryan 

% Start with a clear workspace
clear;
clc;
close all

% Initialise geometry structure 
Geom.b = 5.84;            
Geom.cr = 1.58;
Geom.lambda = 0.5;   
Geom.ct = Geom.lambda*Geom.cr;         
Geom.sweep = 0;         
Geom.twist = 0;         
Geom.dihed = 0;
Geom.camber = 1; % 0 == Flat Plat, 1 == NACA 65(2)-415, 2 == NACA 23012

% Other Variables
alpha = 5;      % AoA (deg)
npx = 11;       % Number of x points
npy = 21;       % Number of y points
toPlot = 1;
L = (npx-1)*(npy-1);    % Number of panels
U = 10;
rho = 1.225;

% Generate Lattice
[Surface, Lattice, Collocation, Normals] = GenerateLattice(Geom, alpha, npx, npy, toPlot);

% Convert to cell functions
% Not the Lattice_cell contains the wake panel indexes in the final column
[Lattice_Cell, Collocation_Cell, Normals_Cell] = Mat2Cell(Lattice, Collocation, Normals);


% Generate ICM matrix 
[ICM_A, ICM_B] = GenerateICM(Normals_Cell, Collocation_Cell, Lattice_Cell, npx, npy);

% Generate Free Stream matrix 
RHS = GenFreeStreamMat(U, alpha, Normals_Cell, L);

% Calculate Gamma vector
Gamma = ICM_A\RHS;
Gamma_plot = reshape(Gamma,[npy-1,npx-1]);
Gamma_diff = [Gamma_plot(:,1), diff(Gamma_plot,1,2)];

% Area
[S, Area_mat] = PlanFormArea(Lattice_Cell, npx, npy);

% Chords at collocaiton points for sectional lift coefficient 
Chord_vec = ChordsAtColl(Surface, npy);

% Calculate lift 
[DeltaLift, Cl_sum, Cl_section, Cl_CL] = CalculateLift(rho, U, Gamma_diff, Lattice, Chord_vec, Area_mat, S, npx, npy);

% Calculate Drag
[DeltaDrag, Cd_sum] = CalculateDrag(ICM_B, Gamma, Gamma_diff, Lattice, S, rho, U, npx, npy);


% Plotting stuff
x = linspace(0,Geom.cr,npx-1);
x_plot = repmat(x,npy-1,1);
figure(3)
hold on
for i = 1:(npy-1)
    plot(x, Gamma_plot(i,:))
end
title('Gamma vs x')

disp(strcat('Total Lift Coefficient     :  ',string(Cl_sum)))
disp(strcat('Total Induced Drag Ceoff   : ',string(Cd_sum)))
disp(strcat('Max Net Circulation        : ',string(max(max(Gamma_diff)))))


% Create Gamma Mesh 
figure(4)
mesh(Collocation.X, Collocation.Y, Gamma_diff)
xlabel('x')
ylabel('y')
zlabel('\Delta\Gamma')
title('Exp Gamma Dist')

figure(5) 
mesh(Collocation.X,Collocation.Y,DeltaLift);
title('Exp Lift Dist')
xlabel('x')
ylabel('y')
zlabel('Lift [N]')

figure(6) 
mesh(Collocation.X, Collocation.Y,DeltaDrag)
title('Exp Induced Drag Dist')
xlabel('x')
ylabel('y')
zlabel('Drag [N]')

function [Surface, Lattice, Collocation, Normals] = GenerateLattice(Geom, AoA, npx, npy, toPlot)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AERO3260 Aerodynamics 1 A2
% Author: 470 390 467
%
% function [Surface, Lattice, Collocation, Normals] = GenerateLattice(Geom, AoA, npx, npy, toPlot)
% Info:
%   This function generates a non-linear mesh lattice for given wing input 
%   geometries. Camber: NACA23012
%
% Inputs:
%   Geom: Structure containing wing geometry. Required shown below:
%           Geom.b;             Semi span                     [m]
%           Geom.cr;            Root Chord                    [m]
%           Geom.ct;            Tip Chord*                    [m]
%           Geom.lambda;        Taper Ratio                    -
%           Geom.sweep;         Sweep angle                   [deg]
%           Geom.twist;         Wing tip twist angle          [deg]
%           Geom.dihed;         Dihedral angle                [deg]
%
%   AoA: Angle of attack [deg]
% 
%   npx: number of discretised points along x-axis
%
%   npy: number of discretised points along y-axis
% 
%   toPlot: Boolean value triggering mesh plots (True = plot)
% 
% Outputs: 
%   Surface: Stucture containing coordinates of surface under three fields:
%               - surface.X     Matrix of x coordinates [npy, npx]    
%               - surface.Y     Matrix of y coordinates [npy, npx]   
%               - surface.Z     Matrix of z coordinates [npy, npx]  
% 
%   Lattice: Stucture containing coordinates of lattice mesh under same
%            fields as above. Size of field: [npy, npx+1]  
%
%   Collocation: Structure containing coordinates of collocation points
%                with same fields as above. Size of field: [npy-1, npx-1]   
%
%   Normals: Structure containing unit vectors normal to the plane at
%            collocation points in the same order as the collocation 
%            structure. Contains three fields:
%                 - Normals.u     x component of unit vector [npy-1, npx-1]
%                 - Normals.v     y component of unit vector [npy-1, npx-1]
%                 - Normals.w     z component of unit vector [npy-1, npx-1]
%
% correct?: Updated to linear 10/11/2020 correct
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    %% 0 Initialise 
    % Exract wing geometry 
    b = Geom.b;                 % Semi span                     [m]
    cr = Geom.cr;               % Root Chord                    [m]
    ct = Geom.ct;               % Tip Chord                     [m]
    lambda = Geom.lambda;       % Taper Ratio                    -
    sweep = Geom.sweep;         % Sweep angle                   [deg]
    twist = Geom.twist;         % Wing tip twist angle          [deg]
    dihed = Geom.dihed;         % Dihedral Angle                [deg]
    camber = Geom.camber;       % Camber identifier flag        -

    %% 1 Discretise 
    % Create discretisation meshes using non linear methods Cosine in x and 
    % half-sine in y 
    y = linspace(0,b,npy);   % LINEAR
    x = linspace(0,cr,npx);  % LINEAR
    [Y_0, X_0] = meshgrid(y, x);
    Z_0 = zeros(npx,npy);   
    
    % Quarter Chord Panel Positions
    % We will calculat the quarter chord positions now and apply all the
    % same transformations as the surface panels. This will make it easier
    % to calculat the coordinates under multiple geometric transformations
    % Vortex lattice x coordinates
    X_lat = 0.25*(X_0(2:end,:)-X_0(1:end-1,:))+X_0(1:end-1,:);
    X_lat = [X_lat; X_lat(end,:) + (X_0(end,:)-X_0(end-1,:))];
    
    %% 2 Camber 
    % This function will use a NACA 5-series camber profile for ease 

 

    switch camber 
        
        % Flat Plate
        case 0
            z_c = 0*X_lat(:,1)';  
            z_c_lat = z_c;
            
        %NACA 65(2)-415
        case 1
            [P, z_c] = NACA65415(npx);
            x_lat = X_lat(:,1)./cr;
            z_c_lat = polyval(P,x_lat)';
            
        % NACA23012
        case 2
            % Preallocate memory for performance       
            x = x./cr;                        % x/c ratio vector
            x_lat = X_lat(:,1);
            p = 0.15;                       
            m = 0.2025;
            k1 = 15.957;
            z_c = zeros(1,length(x));
            z_c_lat = z_c;
            for i = 1:length(x)
                pos = x(i);
                if pos > m
                    z_c(i) = k1/6*m^3*(1-x(i));
                    z_c_lat(i) = k1/6*m^3*(1-x_lat(i));
                else 
                    z_c(i) = k1/6*(x(i)^3-3*m*x(i)^2+m^2*(3-m)*x(i));
                    z_c_lat(i) = k1/6*(x_lat(i)^3-3*m*x_lat(i)^2+m^2*(3-m)*x_lat(i));
                end
            end
            
        otherwise
            error('Incorrect Camber Identifier: 0 == Flat Plat, 1 == NACA 65(2)-415, 2 == NACA 23012')
    end
    
    % Calulate Z values across entire half span
    Z_1 = cr*repmat(z_c',1,npy);
    Z_1_lat = cr*repmat(z_c_lat',1,npy);
    

    

    %% 3 Twist 
    % Apply rotation matrix to coordinates to simulate linear twist
    % +ve angle twists tip nose down
    shift = 0.25;                               % Point of y-axis rotation
    twist = linspace(0,deg2rad(twist));   %LINEAR
    
    % Loop through every y positon
    for i = 1:npy
        
        % Rotation matrix 
        R = [cos(twist(i)), -sin(twist(i)); 
             sin(twist(i)), cos(twist(i))];
         
        % Surface
        xloc = X_0(:,i)';    
        zloc = Z_1(:,i)';
        xShift = shift*(xloc(end)-xloc(1));
        xloc = xloc - xShift;
        % Rotate
        temp = R*[xloc;zloc];
        xt = temp(1,:);
        zt = temp(2,:);

        % Lattice
        xloc_lat = X_lat(:,i)';    
        zloc_lat = Z_1_lat(:,i)';
        xShift_lat = shift*(xloc_lat(end)-xloc_lat(1));
        xloc_lat = xloc_lat - xShift_lat;
        % Rotate
        temp_lat = R*[xloc_lat;zloc_lat];
        xt_lat = temp_lat(1,:);
        zt_lat = temp_lat(2,:);
        
        
        % Re-contruct meshgrid from values
        X_3(:,i) = xt' + xShift;                % Surface X
        Z_3(:,i) = zt' + Z_1(:,i);              % Surface Z
        X_3_lat(:,i) = xt_lat' + xShift_lat;    % Lattice X
        Z_3_lat(:,i) = zt_lat' + Z_1_lat(:,i);  % Lattice Z
    end
    
    
    
    %% 4 Taper
    % Add taper ratio to the wingspan, keeping the trailing edge straight
    chord = cr + (cr*y)/b*(lambda-1);     % Calculate chrod at each y position
    cTaper = repmat(chord./cr, npx, 1);         % Expand to matrix
    
    % Taper X and Z coordinates 
    % NOTE because twist and camber are linear functions of Z,X we can just
    % multiply the taper ratio through each to estimate wing tapering
    % Surface
    %X_4 = (X_3 - cr).*cTaper + cr;  %w.r.t trailing edge
    X_4 = (X_3 - 0.25*cr).*cTaper + 0.25*cr; % w.r.t 1/4 chord 
    Z_4 = Z_3.*cTaper;              %w.r.t bottom surface 
    
    % NOTE lattice mesh cannot be scaled in the same way. We must find the
    % difference between the surface and lattice points, add the taper
    % scaling factor and then add them back onto the surface points to get
    % the tapered lattice points
    X_diff = X_3_lat-X_3;
    Z_diff = Z_3_lat-Z_3;
    % Lattice
    X_4_lat = X_4 + cTaper.*X_diff;
    Z_4_lat = Z_4 + cTaper.*Z_diff;
    
    
    %% 5 Sweep and Dihedral
    % Apply sweep estimation by adding x scalar depending on y position
    dXSweep = tand(sweep)*y;
    X_5 = X_4 + dXSweep;
    X_5_lat = X_4_lat + dXSweep;
    
    % Apply dihedral estimation by adding z scalar depending on y position
    dZDihed = tand(dihed)*y;
    Z_5 = Z_4 + dZDihed;
    Z_5_lat = Z_4_lat + dZDihed;
    
    % NOTE in order to keep these estimations accurate dihedral and sweep
    % angle must be kept small  < 10 deg
    
    %% CREATE OUPUT 
    % Create surface strucrure, transpose to correlate to panel numbering
    Surface.X = X_5';
    Surface.Y = Y_0';
    Surface.Z = Z_5';
    
    
    %% Add wake panel
    xwake = 100;
    X_6_lat = [X_5_lat; X_5_lat(end,:) + xwake];
    Y_6_lat = [Y_0; y];
    Z_6_lat = [Z_5_lat; Z_5_lat(end,:) + xwake*tand(AoA)];  % Include angle of attack
    
    % Create Lattice output
    % transpose to correlate to panel numbering
    Lattice.X = X_6_lat';
    Lattice.Y = Y_6_lat';
    Lattice.Z = Z_6_lat';
    
    %% Collocation Points and Normals
    % NOTE Loop goes through in this order, collocation points will be
    % assigned in this order
    %   
    %        TIP  
    %   |---|---|---|
    %   | 3 | 6 | 9 |
    %   |---|---|---|
    % j | 2 | 5 | 8 |
    %   |---|---|---|
    %   | 1 | 4 | 7 |
    %   |---|---|---|
    %         i
    %        ROOT 
    % 
    % This will be store into the matrix in this order:
    % 
    %   [1, 4, 7]
    %   [2, 5, 8]
    %   [3, 6, 9]
    %
    %
    %
    % Preallocate memroy for speed
    X_col = zeros((npy-1),(npx-1));
    Y_col = zeros((npy-1),(npx-1));
    Z_col = zeros((npy-1),(npx-1));
    u_norm = zeros((npy-1),(npx-1));
    v_norm = zeros((npy-1),(npx-1));
    w_norm = zeros((npy-1),(npx-1));
    x = 1;
    y = 2;
    z = 3;
    for i = 1:(npx-1)
        for j = 1:(npy-1)
            
            % X point
            p_11_x = X_6_lat(i,j); 
            p_21_x = X_6_lat(i+1,j); 
            p_12_x = X_6_lat(i,j+1); 
            p_22_x = X_6_lat(i+1,j+1);
            c_x = mean([p_11_x,p_21_x,p_12_x,p_22_x]);
            
            % Z Point
            p_11_z = Z_6_lat(i,j); 
            p_21_z = Z_6_lat(i+1,j); 
            p_12_z = Z_6_lat(i,j+1); 
            p_22_z = Z_6_lat(i+1,j+1);
            c_z = mean([p_11_z,p_21_z,p_12_z,p_22_z]);
            
            % Y Point
            p_11_y = Y_6_lat(i,j);
            p_12_y = Y_6_lat(i,j+1);
            c_y = mean([p_11_y,p_12_y]);
            
            % Assign to matrix
            X_col(j,i) = c_x;
            Y_col(j,i) = c_y;
            Z_col(j,i) = c_z;
            
            % Create veoctors to calculate normals
            A = [X_6_lat(i+1,j+1)-X_6_lat(i,j), Y_6_lat(i+1,j+1)-Y_6_lat(i,j), Z_6_lat(i+1,j+1)-Z_6_lat(i,j)];
            B = [X_6_lat(i,j+1)-X_6_lat(i+1,j), Y_6_lat(i,j+1)-Y_6_lat(i+1,j), Z_6_lat(i,j+1)-Z_6_lat(i+1,j)];

            % Calculate cross product
            AxB_x = A(y)*B(z)-A(z)*B(y);
            AxB_y = A(z)*B(x)-A(x)*B(z);
            AxB_z = A(x)*B(y)-A(y)*B(x);
            
            % Calculate Magnitude
            AxB_abs = sqrt(AxB_x^2+AxB_y^2+AxB_z^2);
            
            % Assign to matrix 
            u_norm(j,i) = AxB_x./AxB_abs;
            v_norm(j,i) = AxB_y./AxB_abs;
            w_norm(j,i) = AxB_z./AxB_abs;
        end
    end
    
    % Create output 
    Collocation.X = X_col;
    Collocation.Y = Y_col;
    Collocation.Z = Z_col;
    Normals.u = u_norm;
    Normals.v = v_norm;
    Normals.w = w_norm;
    
   
    %%  Make plots 
    if toPlot 
        mesh(X_5,Y_0,Z_5) 
        alpha(0.5);
        axis equal
        hold on 
        plot3(X_col,Y_col,Z_col,'rx')
        %mesh(X_6_lat,Y_6_lat,Z_6_lat,'LineStyle','--') Include wake panel
        plot3(X_5_lat,Y_0,Z_5_lat,'ob')
        alpha(.5)
        quiver3(Collocation.X,Collocation.Y,Collocation.Z,Normals.u,Normals.v,Normals.w,0.4)
        view(2)  
        xlabel('X [m]')
        ylabel('Y [m]')
        zlabel('Z [m]')
    end    
end
    
function [u,v,w] = vortexLine(x,y,z,x1,y1,z1,x2,y2,z2,gamma,eps)
%% Subroutine for computing the induced velocity at some point P (x,y,z) due
% to a vortex line segment between P1 (x1, y1, z1) and P2 (x2, y2, z2) 
%
% Inputs
% 
% Output 
%
% Author: 470390467

    % Compute cross product
    r1r2c_x = (y-y1)*(z-z2)-(z-z1)*(y-y2);
    r1r2c_y = -(x-x1)*(z-z2)+(z-z1)*(x-x2);
    r1r2c_z = (x-x1)*(y-y2)-(y-y1)*(x-x2);

    % Compute abs of vector product
    r1r2abs = (r1r2c_x^2 + r1r2c_y^2 + r1r2c_z^2);

    % Compute distances
    r1 = sqrt((x-x1)^2+(y-y1)^2+(z-z1)^2);
    r2 = sqrt((x-x2)^2+(y-y2)^2+(z-z2)^2);

    % Check for singular values
    if (r1 < eps) || (r2 < eps) || (r1r2abs < eps)
        u = 0;
        v = 0;
        w = 0;
    else 
        % Compute dopt products
        r0r1d = (x2-x1)*(x-x1)+(y2-y1)*(y-y1)+(z2-z1)*(z-z1);
        r0r2d = (x2-x1)*(x-x2)+(y2-y1)*(y-y2)+(z2-z1)*(z-z2);

        % Compute scalar constant
        K = gamma/4/pi/r1r2abs *(r0r1d/r1-r0r2d/r2);

        % Compute velocity components 
        u = K*r1r2c_x;
        v = K*r1r2c_y;
        w = K*r1r2c_z;

    end
end

function [u,v,w] = vortexParrallel(P, points, gamma,eps)
    %% Compute induced velocty for two parrallel lines used in the drag
    % calculation
    % placed between P1?P2, P2?P3, P3?P4, P4?P1.

    % Define some index variables
    x = 1;
    y = 2; 
    z = 3;

    % Unpack input matrix 
    P1 = points(1,:);
    P2 = points(2,:);
    P3 = points(3,:);
    P4 = points(4,:);

    [u1,v1,w1] = vortexLine(P(x),P(y),P(z),P1(x),P1(y),P1(z),P2(x),P2(y),P2(z),gamma,eps);
    [u3,v3,w3] = vortexLine(P(x),P(y),P(z),P3(x),P3(y),P3(z),P4(x),P4(y),P4(z),gamma,eps);


    % Sum the induced velocities
    u = u1 + u3;
    v = v1 + v3;
    w = w1 + w3;

end


function [u,v,w] = vortexRing(P, ring, gamma,eps)
%% Compute the induced velocity due to a rectilinear vortex ring by use of line segments
% placed between P1?P2, P2?P3, P3?P4, P4?P1.

    % Define some index variables
    x = 1;
    y = 2; 
    z = 3;

    % Unpack input matrix 
    P1 = ring(1,:);
    P2 = ring(2,:);
    P3 = ring(3,:);
    P4 = ring(4,:);

    [u1,v1,w1] = vortexLine(P(x),P(y),P(z),P1(x),P1(y),P1(z),P2(x),P2(y),P2(z),gamma,eps);
    [u2,v2,w2] = vortexLine(P(x),P(y),P(z),P2(x),P2(y),P2(z),P3(x),P3(y),P3(z),gamma,eps);
    [u3,v3,w3] = vortexLine(P(x),P(y),P(z),P3(x),P3(y),P3(z),P4(x),P4(y),P4(z),gamma,eps);
    [u4,v4,w4] = vortexLine(P(x),P(y),P(z),P4(x),P4(y),P4(z),P1(x),P1(y),P1(z),gamma,eps);

    % Sum the induced velocities
    u = u1 + u2 + u3 + u4;
    v = v1 + v2 + v3 + v4;
    w = w1 + w2 + w3 + w4;
end



function [Lattice_Cell, Collocation_Cell, Normals_Cell] = Mat2Cell(Lattice, Collocation, Normals)
%% AERO3260 Aerodnynamics 1 A2
% Author: 470390467
% 
% [Lattice_Cell, Collocation_Cell, Normals_Cell] = Mat2Cell(Lattice, Collocation, Normals)
%
% Info:
%   This function takes the Lattice, Collocation and Normal structures and
%   converts them into cell arrays that can be accessed in order of panels
%   definition. The panel definition is shown below:
%                               TIP  
%                          |---|---|---|
%                          | 3 | 6 | 9 |
%                          |---|---|---|
%                        j | 2 | 5 | 8 |
%                          |---|---|---|
%                          | 1 | 4 | 7 |
%                          |---|---|---|
%                                i
%                               ROOT 
%   The output decription is described below but the aim is to access panel
%   information by a single parameter n = 1:L, L is number of panels. Input
%   matricies are in the form:
%                           [1, 4, 7]
%                           [2, 5, 8]
%                           [3, 6, 9]
%
% Inputs: 
%   Lattice: Structure containing X,Y and Z matrixes of vortex panel corner
%   coordinates
%
%   Collocation: Structure containing X,Y and Z matrixes of collocation
%   point coordinates
%
%   Normals: Structure containing u,v and w matrixes of normal unit vectors
%   components 
%
%
% Outputs: 
%   Latttice_cell: [npy-1, npx] cell array containing 4x3 matricies
%                  containing points of vortex lattice corners. NOTE the
%                  wake panel lattice coordinates are included in the final
%                  column of this array.
%   Create Lattice Cell Vector with corner designation: 
% 
%                  (1)---(2)   
%                 j |     |
%                  (4)---(3)     
%                      i 
%    Each panel at i,j indexes a 4x3 matrix of corner points in the form:
%             
%               [x_1, y_1, z_1;
%                x_2, y_2, z_2;
%                x_3, y_3, z_3;
%                x_4, y_4, z_4]
%
%    Collcation_cell: [npy-1, npx-1] cell array containing 1x3 vector of
%                    collocation points coordinates [x,y,z]
%   Normals_cell: [npy-1, npx-1] cell array containing 1x3 vector of normal
%                 unit vector components [u,v,w]
%
% correct?: yes, think so. Careful with matrix indexes. 3/11/2020

    % Find npx and npy 
    [r, c] = size(Normals.u);

    % Preallocate memory for performance
    Lattice_Cell = cell(r,c);
    Collocation_Cell = cell(r,c);
    Normals_Cell = cell(r,c);

    % Loop through every panel 
    for j = 1:(c+1)
        for i = 1:r 

            % Create Lattice Cell Vector
            p1 = [Lattice.X(i+1,j), Lattice.Y(i+1,j), Lattice.Z(i+1,j)];
            p2 = [Lattice.X(i+1,j+1), Lattice.Y(i+1,j+1), Lattice.Z(i+1,j+1)];
            p3 = [Lattice.X(i,j+1), Lattice.Y(i,j+1), Lattice.Z(i,j+1)];
            p4 = [Lattice.X(i,j), Lattice.Y(i,j), Lattice.Z(i,j)];
            Lattice_Cell{i,j} = [p1; p2; p3; p4];

            if ~(j > c)
                % Create Collocation Cell Vector
                Collocation_Cell{i,j} = [Collocation.X(i,j),Collocation.Y(i,j), Collocation.Z(i,j)];

                % Create Normals Cell Vector
                Normals_Cell{i,j} = [Normals.u(i,j),Normals.v(i,j), Normals.w(i,j)];
            end
        end
    end
end
    
function [ICM_A, ICM_B] = GenerateICM(Normals_Cell, Collocation_Cell, Lattice_Cell, npx, npy)
%% AERO3261 Aerodynamics 1 A2
% Author: 470 390 467
%
% Info:
%   This function generates the Influence Ceofficient Matrix which is used
%   to calculate the circulation. 
%
% Inputs:
%   Normals_Cell
%   Collocation_Cell
%   Lattice_Cell
%   npx: number of x points
%   npy: number of y points
% 
% Outputs:
%   ICM: LxL Influence Coefficient Matrix 
%
% Correct?: Correct! Fucking finally

eps = 1e-9;   % Vortex element width - maybe make this a fn of panel area? eps = min(min(area))*10^-3
gamma = 1;  % Circulation (set to 1 to normalise w.r.t gamma)
M = npx-1;  % Chord-wise panels
N = npy-1;  % Span-wise panels
L = N*M;    % Total number of panels

% Preallocate memory for performance
ICM_A = zeros(L);         % Square LxL matrix 
ICM_B = zeros(L);

    % Loop through every panel^2
    for p = 1:L % Panel m
        for i = 1:L % Vortex that is influencing panel n 

            % Normal vector
            n = Normals_Cell{p};

            % Collocation Point
            coll = Collocation_Cell{p};

            % Panel coordinates 
            ring = Lattice_Cell(i);

            % Convert ring to matrix 
            ring = cell2mat(ring);

            % Induced velocity 
            [u_r, v_r, w_r] = vortexRing(coll, ring, gamma, eps);
            [u_l, v_l, w_l] = vortexRing(coll, LeftSide(ring), -gamma, eps);
            [ub_r, vb_r, wb_r] = vortexParrallel(coll, ring, gamma, eps);
            [ub_l, vb_l, wb_l] = vortexParrallel(coll, LeftSide(ring), -gamma, eps);

            % Influence coefficient of vortex ring i at panel p
            A_curr = dot([u_r+u_l,v_l-v_r,w_r+w_l],n);
            B_curr = dot([ub_r+ub_l,vb_l-vb_r,wb_r+wb_l],n);

            % Add to ICM matrix 
            ICM_A(p,i) = A_curr;
            ICM_B(p,i) = B_curr;

            % Check if influence panel is on the trailing edge
            if i > (L-N)

                % If so add the wake panel term aswell
                ring_wake = Lattice_Cell(i+N);
                ring_wake = cell2mat(ring_wake);
                [u_wr, v_wr, w_wr] = vortexRing(coll, ring_wake, gamma, eps);
                [u_wl, v_wl, w_wl] = vortexRing(coll, LeftSide(ring_wake), -gamma, eps);
                [ub_wr, vb_wr, wb_wr] = vortexParrallel(coll, ring_wake, gamma, eps);
                [ub_wl, vb_wl, wb_wl] = vortexParrallel(coll, LeftSide(ring_wake), -gamma, eps);

                % Wake panel influence coefficient
                A_wake = dot([u_wr+u_wl,v_wl-v_wr,w_wr+w_wl], n);
                B_wake = dot([ub_wr+ub_wl,vb_wl-vb_wr,wb_wr+wb_wl], n);

                % Add the wake influence coefficient to the trailing edge panel
                ICM_A(p,i) = ICM_A(p,i) + A_wake;
                ICM_B(p,i) = ICM_B(p,i) + B_wake;
            end

        end
    end
end
    
function RHS = GenFreeStreamMat(U, alpha, Normals_Cell, L)
%% AERO3260 Aerodynamics 1 A3
% Author: 470 390 467 
% 
% Info
%   This function calcultes the free stream matrix on the RHS of the
%   impermiability boundary condition equaiton in tutorial 10:
%   RHS = -[u_inf,u_inf,w_inf] * n
%
% Inputs;
%   U: Free Stream magnitude [m/s]
%   alpha: Angle of Attack [deg]
%   Normals_Cell: Cell Array of normal unit vectors for each panel
%   L: Total number of Panels
%
% Outputs:
%   RHS: Free Stream Influence Matrix 
% 
% Correct? YES 

    % Convert alpha to radians
    alpha = deg2rad(alpha);

    % Calculate free stream velocity components 
    u_inf = U*cos(alpha);
    v_inf = 0;
    w_inf = U*sin(alpha);
    U_inf = -[u_inf, v_inf, w_inf];
    U_inf = repmat(U_inf, L, 1);

    % Get colum vector of normals 
    Normals_col = cell2mat(Normals_Cell(:));

    % Evaluate RHS free stream matrix vector
    RHS = dot(U_inf',Normals_col')';
end

function [S, Areas_mat] = PlanFormArea(Lattice_Cell, npx, npy)
%% Calculates total planform area S as well as a matrix of all panle areas
    i = 1;
    Areas_mat = zeros(npy-1,npx-1);
    for x = 1:(npx-1)
        for y = 1:(npy-1)
            ring = Lattice_Cell{i}; 
            area_i = PanelArea(ring);
            Areas_mat(y,x) = area_i;
            i = i + 1;
        end
    end
    S = sum(sum(Areas_mat));
end

function Chord_vec = ChordsAtColl(Surface, npy)
%% This function calcualtes the chord at the collocation points along y axis
% by averaging the chord at j+1 and j. This can be done assuming a linear
% resolution mesh.
% Ouptut: (npy-1)x1 column vec of chords at collocation points

    Chord_vec = zeros(npy-1,1); % Column vector
    X = Surface.X;

    % Loop through every Collocation point
    for j = 1:(npy-1)

        % Chord at j+1
        chord_plus = X(j+1,end)-X(j+1,1);

        % Chord at j
        chord = X(j,end)-X(j,1);

        % Average between the two
        Chord_vec(j) = (chord+chord_plus)/2;
    end
end

function [DeltaLift, CL_sum, Cl_section, Cl_CL] = CalculateLift(rho, U, Gamma_diff, Lattice, Chord_vec, Area_mat, S, npx, npy)
%% AERO3260 Aerodynamics 1 A2
% Author: 470 390 467
% 
% Info:
%   This function calculates the lift distribution and sum from the vortex 
%   lattice ciruculation values calculated before.
%
% Inputs:
%   rho:    density
%   U:      Free stream velocity magnitude
%   Gamma_diff:  Matrix of gamma values at each 0.25x_ij location,
%   including the difference due to shared vortex line
%   Lattice: vortex lattice structure output from GenerateLattice() func
%   L: number of panels
%
% Outputs: 
%   Lift_mat: (nnpx-1)x(npy-1) matrix of Delta_lift values
%   Lift_sum: Overall lift of the wing
%
% Correct? Correct - Sectional lift coefficient is correct I THINK 


    % DeltaY values 
    L = (npx-1)*(npy-1);
    yLin = Lattice.Y(:,1);
    Y_diff = diff(yLin);
    
    % Reshape gamma
    Gamma_diff_col = Gamma_diff(:);
    
    % Pre-Allocate memory for speed
    DeltaLift_col = zeros(L,1);
    
    CL_sum = 0;
    i = 1; % Panel Counter
    for x = 1:(npx-1)
        for y = 1:(npy-1)
            
            % DeltaLift Col
            DeltaLift_col(i) = rho*U*Gamma_diff_col(i)*Y_diff(y);
            
            % Lift coefficient 
            delta = 2*DeltaLift_col(i)/rho/U.^2/S;  
            CL_sum = CL_sum + delta;
            
            % Update current panel
            if i <= L
                i = i + 1;
            else
                error('Panel counter exceeding max panels L')
            end
        end
    end
    % Reshape Delta_Lift matrix
    DeltaLift = reshape(DeltaLift_col,[npy-1,npx-1]);
    
    % Preallocate memory
    Cl_section = zeros(npy-1,1);
    
    
    % Find section lift coefficient 
    for y = 1:(npy-1)
        
        % Sectional lift coefficient 
        Cl_section(y) = sum(DeltaLift(y,:)*2/rho/U.^2/Chord_vec(y)); 
    end
    
    % 3D Lift coefficient
    CL_section = CL_sum*sum(Area_mat,2);
    
    % Cl/CL vector
    Cl_CL = Cl_section./CL_section;
end

function [DeltaDrag, Cd_sum] = CalculateDrag(ICM_B, Gamma, Gamma_diff, Lattice, S, rho, U, npx, npy)
%% AERO3260 Aerodynamics 1 A2
% Author: 470 390 467

    % Calculate the induced down wash vector 
    W_vec = ICM_B*Gamma;
    
    % Define some parameters
    i = 1;                  % Pamel counter
    Cd_sum = 0;             % Induced Drag Coefficient 
    yLin = Lattice.Y(:,1);  
    Y_diff = diff(yLin);    % Panel Chord lengths
    
    % Preallocate memory 
    L = (npx-1)*(npy-1);
    DeltaDrag = zeros(L,1);
    
    % Loop through every panel
    for x = 1:(npx-1)
        for y = 1:(npy-1)

            % Drag increment 
            DeltaDrag_col(i) = -rho*W_vec(i)*Gamma_diff(i)*Y_diff(y);
            
            % Drag Coefficient 
            delta = 2*DeltaDrag_col(i)/rho/U.^2/S; 
            Cd_sum = Cd_sum + delta;
            
            % Update Panel Counter
            if i <= L
                i = i + 1;
            else
                error('Panel counter exceeding max panels L')
            end
        end
    end
    
    % Reshape Delta_Lift matrix
    DeltaDrag = reshape(DeltaDrag_col,[npy-1,npx-1]);
end 

function [P, z_c] = NACA65415(npx)
%% Calculates mean camber line for the NACA 65(2)-415
% P is polynomial coefficients for trend fitting line (6th order)
% z_c is the camber line z coordinates normalised by root chrod
% Author: 470 390 467

    % DATA from Airfoil Tools 
    %             {x/c      , y/c}
    Camber = [0.0023, 0.011, 0.0163, 0.0196, 0.0214, 0.022, 0.0215, 0.0196 ,0.0161, 0.0105 ,0];
    interpData = [0.0020    0.4780;
                  0.0110    0.2170;
                  0.0980    1.0000;
                  0.1420    1.2610;
                  0.1800    1.5220;
                  0.2210    1.6520;
                  0.2970    1.9130;
                  0.3980    2.0430;
                  0.4970    2.0430;
                  0.5970    2.1740;
                  0.6450    2.0430;
                  0.6970    1.7830;
                  0.7460    1.7830;
                  0.7980    1.6520;
                  0.8430    1.2610;
                  0.8990    1.1300;
                  0.9410    0.6090;
                  0.9980         0];
    x = linspace(0,1,npx);
    P = polyfit(x,Camber,6);
    z_c = polyval(P,x);
%     plot(interpData(:,1),polyval(P,interpData(:,1)))
%     hold on
%     plot(x,z_c,'o')
%     axis equal
%     grid on
%   


end

function [ring_LS] = LeftSide(ring)
%% AERO3260 Aerodynamics 1
% Author: 470 390 467
%
% Info
%   This function takes RHS ring coordinate matrix 'ring' and negates the y
%   components, mirroring it to the left hand side of the aircraft
% ring is a 4x3 matrix containing coordinates of vortex lattice panel
    
    y_coor = ring(:,2);
    ring_LS = ring;
    ring_LS(:,2) = -y_coor;
end

function area = PanelArea(ring)
%% AERO3260 Aerodynamics 1 A3
% Author: 470 390 467
% 
% Info:
%   This function takes a set of coordinates describing a trapezium in
%   3D and calculates the surface area between all four points
%
% Inputs: 
%   ring: 4x3 matrix of ring corner coordinates in 3D
% 
% Outputs:
%   area: Area of current panel [m^2]
%
% Correct? Should be - NOT TESTED

    P1 = ring(1,:);
    P2 = ring(2,:);
    P3 = ring(3,:);
    P4 = ring(4,:);

    % Treat panel as trapezium and calculate area 
    a_vec = P1-P2;
    b_vec = P4-P3;
    a = sqrt(a_vec(1).^2+a_vec(2).^2+a_vec(3).^2);
    b = sqrt(b_vec(1).^2+b_vec(2).^2+b_vec(3).^2);
    h = abs(P4(2) - P1(2));
    area = h/2*(a+b);
end





% Woo Hoo 1000 lines of code! 







