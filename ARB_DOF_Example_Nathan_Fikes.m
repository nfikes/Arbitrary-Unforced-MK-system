%% Example MK system Arbitrary DOF
% File Name: Coupled EOM decoupling

% Name: Nathan Fikes 
% Virginia Polytechnic Institute and State University
% Class: ME 3524
% Program: Example for solving coupled EOM, massive one

clc;
clear;
close all;

%% System Parameters

% Mass of each mass.
m = [2.4    ;
     9.3    ;
     10.4   ;
     12.3   ;
     10     ;
     27.4   ;
     20.7   ;
     2.6    ;
     3.5    ;
     2.4    ];

% The stiffness of the spring to the left of each mass.
k = [534    ;
     296    ;
     134    ;
     850    ;
     250    ;
     123    ;
     194    ;
     823    ;
     240    ;
     194    ];

% The origin of each mass at static equilibrium.
massorig = [10;
            20; 
            30; 
            40; 
            50;
            60;
            70;
            80;
            90;
            100];

%ARB DOF IC (Both the Position and Velocity).
x0 = [0; 1; -5; 1; 0; 2; -6; -3; 4; 5];
xdot0 = [13; -5; 12; 2; 9; 1; 4; -5; 7; -5];


%{
m = [10      ;
     100     ];

k = [80     ;
     60     ];
 
massorig = [10;
            20];

%ARB DOF IC
x0 = [-1; 2];
xdot0 = [0; 0];
%}
 
DOF = length(k);

simulation_speed = 60;

tspan = [0:0.001:20];

% Extra Element in K to help with Stiffness Matrix Calculation
k = [k; 0];

%% Guard Statement
if length(m) ~= (length(k) - 1)
    disp("ERROR: Mass and Spring vector sizes do not agree.")
    disp("Please check the dimensions.")
    return
end


%% Solving the Arbitrary DOF system (UNFORCED)
%Mass Matrix
M = [];
for i = 1:DOF
    zero11_line = [zeros(1, i-1) 1 zeros(1, length(m)-i)]';
    mass_indexed = zero11_line*m(i);
    M = [M mass_indexed];
end
% Example of M matrix:
%   M = [m1,  0,  0,  0;
%         0, m2,  0,  0;
%         0,  0, m3,  0;
%         0,  0,  0, m4;



%Symmetric Stiffness Matrix
K_untrimmed = [];
for i = 1:DOF
    zero113_line = [zeros(1, i-1) 1 0 0 zeros(1, length(m)-i)]';
    zero123_line = [zeros(1, i-1) 0 1 0 zeros(1, length(m)-i)]';
    zero133_line = [zeros(1, i-1) 0 0 1 zeros(1, length(m)-i)]';
    
    stiffness_indexed = zero113_line * -k(i) + zero123_line * (k(i) + k(i + 1)) + zero133_line * -k(i + 1);
    
    K_untrimmed = [K_untrimmed stiffness_indexed];
end
K = K_untrimmed(2:(length(k)), 1:(length(k)-1));
%Example of K matrix:
%   K = [k1 + k2,     -k2,       0,   0;
%            -k2, k2 + k3,     -k3,   0;
%              0,     -k3, k3 + k4, -k4;
%              0,       0,     -k4,  k4;]



%Intermediate Mass Matrix for Mass Normalization
M_HInv = M^(-1/2);  

%Packaged system dynamics in a mass normalized stiffness matrix
K_Te = M_HInv * K * M_HInv;

%Eigensolution on mass normalized stiffness matrix to allow the 
%representation of it as a diagonal matrix.
[V_Pr, D] = eig(K_Te);

%Orthonormalized Eigenvector matrix.
P = [];
for i = 1:DOF
    normalized_vector = V_Pr(:,i) / sqrt( V_Pr(:,i)' * V_Pr(:,i) );
    P = [P normalized_vector];
end

%Modal Matrix.
S = M_HInv * P; %A portal to get from Modal Coordinates to Physical Coordinates
invS = inv(S);  %A portal to get from Physical Coordinates to Modal Coordinates

%Modal IC conversion from physical coordinates.
r0 = invS * x0;
rdot0 = invS * xdot0;

% Unforced Homogeneous solution: The constants using MODAL IC
A = r0;

B = [];
for i = 1:DOF
B = [B rdot0(i) ./ sqrt(D(i,i))];
end

%Modal Coordinates Unforced:
RT_untransposed = [];
for i = 1:DOF
    rt = @(t) A(i)*cos(sqrt(D(i,i))*t) + B(i)*sin(sqrt(D(i,i))*t);
    RT_untransposed = [RT_untransposed rt(tspan)'];
end
RT = RT_untransposed';

% Transform the modal coordinate values into physical coordinates
% using the modal matrix.
xt = S * RT;

masize = 0.5;

figure_ARB_DOF = figure("NumberTitle","off", "Name","Arbitrary Degree of Freedom System", "Position",[0 0 1000 1000], "Color", "#DCBAF9");
movegui(figure_ARB_DOF, "southeast");

colors = [[0.75 0 0]' [0.75 0.75 0]' [0.75 0 0.75]' [0 0.75 0.75]' [0 0 0.75]' [0.75 0.5 0.75]' [0.5 0.75 0.75]' [0.5 0.5 0.75]' [0.25 0.75 0.75]' [0.25 0.25 0.75]']';

pause(10)
for i = 1:simulation_speed:length(tspan)
    clf()
    subplot(3,1,1);
    title('All Individual Positions Visualized as displacements in 1 dimension')
    hold on
    grid on
    for j = 1:DOF
        plot(tspan, xt(j,:), 'Color', colors(j,:));
    end
    xline(tspan(i), 'r-.')

    subplot(3,1,2);
    title('All Individual Positions Visualized as a physical system')
    ylim([-(max(massorig)+20)/10, (max(massorig)+20)/10]);
    xlim([min(massorig)-10,max(massorig)+20]);
    hold on
    grid on

    %Plotting the Masses
    for j = 1:DOF
        x_vector = xt(j,:);
        rectangle('Position',[massorig(j) + x_vector(i) - masize, -masize , 2*masize, 2*masize],'FaceColor',colors(j,:),'EdgeColor','k','LineWidth',1)
    end

    mod_massorig = [0; massorig];
    %Plotting the Lines
    for j = 0:DOF-1
        if j ==0
            x_vector_1 = xt(1,:);
        else
            x_vector_1 = xt(j,:);
        end
        x_vector_2 = xt((j+1),:);
        if j == 0
            line([0 massorig(j+1) + x_vector_1(i) - masize], [0 0],'Color',colors(j+1,:),'LineStyle','--')
        else
        line([massorig(j) + x_vector_1(i) + masize massorig(j+1)  + x_vector_2(i) - masize], [0 0],'Color',colors((j+1),:),'LineStyle','--')
        end
    end

    subplot(3,1,3);
    title('Trucated Modal Responses (r1, r2...)')
    hold on
    grid on
    %Trucated Modal Response
    for j = 1:DOF
        plot(tspan, RT(j,:));
    end
    xline(tspan(i), 'r-.')

    pause(1/60)
end
