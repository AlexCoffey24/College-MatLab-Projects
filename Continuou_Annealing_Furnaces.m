%Variables:
T_melt = 1356;%K
T_anneal = T_melt*0.6%K
T_coil = T_melt*0.4 %K
T_min = .98*T_anneal; %K
T_max = 1.02*T_anneal; %K

T_f1 = 800;%K
T_f2 = 800;%K
T_f3 = 900;%K
T_f4 = T_anneal;%K
V_b = .5/8;%m/s
%3.000 to 3.5 (as of now) = 8 sec

%Constants:
p1 = 100;%slices
p2 = 250;%slices
p3 = 250;%slices
p4 = 100;%slices
p5 = 500;%slices
T_room = 298;%K
F_top = 0.2;
F_wall = 0.4;
H_m = 0.0005;%m
W_m = 0.4;%m
em_m = 0.25;
h_f = 1;%m
w_f = 1;%m
L_f1 = 0.5;%m
L_f2 = 1.25;%m
L_f3 = 1.25;%m
L_f4 = 0.5;%m
L_5 = 2.5; %m outside of furnace
em_f = 0.75;
C_p = 380;%J/Kg-K
rho = 8800;%Kg/m^3
sig = 5.68e-8;
in1 = zeros(1,p1); %input(iterations,position)
in2 = zeros(1,p2); %input(iterations,position)
in3 = zeros(1,p3); %input(iterations,position)
in4 = zeros(1,p4); %input(iterations,position)
in5 = zeros(1,p5);

T_i1 = 298;
dx = (L_f1)/(p1);%m
dt = dx/V_b;%s
t1 = L_f1/V_b;%s
t2 = L_f2/V_b;%s
t3 = L_f3/V_b;%s
t4 = L_f4/V_b;%s
t5 = L_5/V_b;%s

A_f1 = 4*L_f1+2;%m^2
A_f2 = 4*L_f2+2;%m^2
A_f3 = 4*L_f3+2;%m^2
A_f4 = 4*L_f4+2;%m^2
A_p = dx*W_m*2;%m^2    why was this A_p again can we make it dA???
R1 = ((1-em_f)/(A_f1*em_f))+(1/A_f1)+((1-em_m)/(A_p*em_m));
R2 = ((1-em_f)/(A_f2*em_f))+(1/A_f2)+((1-em_m)/(A_p*em_m)); 
R3 = ((1-em_f)/(A_f3*em_f))+(1/A_f3)+((1-em_m)/(A_p*em_m));
R4 = ((1-em_f)/(A_f4*em_f))+(1/A_f4)+((1-em_m)/(A_p*em_m));
P_p = 2*dx + 2*W_m;
L = A_p/P_p;



clf
figure(1)
hold on



%Furnace 1
in1(1,1) = T_room;

for m = 2:p1
    
    B = in1(1,m-1);%before position
    
    %h values
        DT = T_f1 - B;%K
        %top face
        if (5e-5 < L^3 * DT < 0.30)
            h_t = 1.32 * (DT/L)^(1/4);
        elseif (0.30 < L^3 * DT < 470)
            h_t = 1.52 * (DT)^(1/3);
        else
           disp("coeff not in range");
        end
        %bottom face
        if 0.005 < (L^3 * DT) < 470
            h_b = 0.59 * (DT/L)^(1/4);
        else
          disp("coeff not in range");
        end

    in1(1,m) = B + (dt/(dx*W_m*H_m*rho*C_p))*((((T_f1^4)-(B^4))*sig/R1)+((T_f1-B)*A_p*(1/(1/h_t + 1/h_b))));

    
    
end

x1 = linspace(0,L_f1,p1);
plot(x1, in1(1,:),'b')
xlim([0,6])
ylim([0,1200])
title('Distance vs Temp FDM_C_u')
xlabel('Distance (m)')
ylabel('Temp (K)')



%Furnace 2
in2(1,1) = in1(1,p1);

for m = 2:p2
    
    D = in2(1,m-1);%before position
    %h values
        DT = T_f2 - D;%K
        %top face
        if (5e-5 < L^3 * DT < 0.30)
            h_t = 1.32 * (DT/L)^(1/4);
        elseif (0.30 < L^3 * DT < 470)
            h_t = 1.52 * (DT)^(1/3);
        else
           disp("coeff not in range");
        end
        %bottom face
        if 0.005 < (L^3 * DT) < 470
            h_b = 0.59 * (DT/L)^(1/4);
        else
          disp("coeff not in range");
        end
    in2(1,m) = D + (dt/(dx*W_m*H_m*rho*C_p))*(((T_f2^4-D^4)*sig/R2)+((T_f2-D)*A_p*(1/(1/h_t + 1/h_b))));

end

x2 = linspace(L_f1,L_f2+L_f1,p2);
plot(x2, in2(1,:),'c')



%Furnace 3
in3(1,1) = in2(1,p2);

for m = 2:p3
    
    E = in3(1,m-1);%before position
    DT = T_f3 - E;%K
        %top face
        if (5e-5 < L^3 * DT < 0.30)
            h_t = 1.32 * (DT/L)^(1/4);
        elseif (0.30 < L^3 * DT < 470)
            h_t = 1.52 * (DT)^(1/3);
        else
           disp("coeff not in range");
        end
        %bottom face
        if 0.005 < (L^3 * DT) < 470
            h_b = 0.59 * (DT/L)^(1/4);
        else
          disp("coeff not in range");
        end
    in3(1,m) = E + (dt/(dx*W_m*H_m*rho*C_p))*(((T_f3^4-E^4)*sig/R3)+((T_f3-E)*A_p*(1/(1/h_t + 1/h_b))));

end

x3 = linspace(L_f2+L_f1,L_f3+L_f2+L_f1,p3);
plot(x3, in3(1,:),'m')



%Furnace 4
in4(1,1) = in3(1,p3);

for m = 2:p4
    
    G = in4(1,m-1);%before position
    DT = 1*(T_f4 - G);%K change on furnace temp
        %top face
        if (5e-5 < L^3 * DT < 0.30)
            h_t = 1.32 * (DT/L)^(1/4);
        elseif (0.30 < L^3 * DT < 470)
            h_t = 1.52 * (DT)^(1/3);
        else
           disp("coeff not in range");
        end
        %bottom face
        if 0.005 < (L^3 * DT) < 470
            h_b = 0.59 * (DT/L)^(1/4);
        else
          disp("coeff not in range");
        end
    in4(1,m) = G + (dt/(dx*W_m*H_m*rho*C_p))*(((T_f4^4-G^4)*sig/R4)+((T_f4-G)*A_p*(1/(1/h_t + 1/h_b))));

end

x4 = linspace(L_f3+L_f2+L_f1,L_f4+L_f3+L_f2+L_f1,p4);
plot(x4, in4(1,:),'-r')

%Outside of furnace
in5(1,1) = in4(1,p4);
for m = 2:p5
    
       J = in5(1,m-1);%before position
       
        DT = J - T_room; %K
        %top face
        if (5e-5 < L^3 * DT < 0.30)
            h_t = 1.32 * (DT/L)^(1/4);
        elseif (0.30 < L^3 * DT < 470)
            h_t = 1.52 * (DT)^(1/3);
        else
           disp("coeff not in range");
        end
        %bottom face
        if 0.005 < (L^3 * DT) < 470
            h_b = 0.59 * (DT/L)^(1/4);
        else
          disp("coeff not in range");
        end
       
    in5(1,m) = J + (dt/(dx*W_m*H_m*rho*C_p))*((sig*A_p*em_m*-1*((J^4)-(T_room^4))+(-1*(J - T_room)*A_p*(1/(1/(4*h_t) + 1/(4*h_b))))));
    
end
disp(in5(1,p5))
x5 = linspace(L_f4+L_f3+L_f2+L_f1,L_5+L_f4+L_f3+L_f2+L_f1,p5);
plot(x5, in5(1,:),'-g') 


