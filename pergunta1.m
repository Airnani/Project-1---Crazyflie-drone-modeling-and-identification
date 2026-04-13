% Unidade Curricular: UAV's
% Project 1 - Crazyflie drone modeling and identification

close;clear all; clc;

% Dados
m = 0.045;                      
g = 9.81;                       
L = 0.046;                      
l = L * sin(pi/4);              
CQ = 7.75e-11;                  
CT = 3.72e-8;                   
Ctau = CQ/CT;                   
J = diag([2.4e-5, 2.4e-5, 3.2e-5]); 
D = diag([7.24e-10, 7.24e-10, 7.24e-10]);

% MATRIZ DE MAPEAMENTO

M = [ 1,      1,      1,      1;
     -l,     -l,      l,      l;
     -l,      l,      l,     -l;
     Ctau, -Ctau,  Ctau, -Ctau];

% ENTRADA:Thrust de Equilíbrio (Hover)
T_individual = ((m * g) / 4)*1.25 ; 
u_motores = [T_individual; T_individual; T_individual; T_individual];

u_controlo = M * u_motores;
T = u_controlo(1);      
tau = u_controlo(2:4);  

% CONFIGURAÇÃO DA SIMULAÇÃO
tspan = [0 1.5];        
x0 = zeros(12, 1);      

% Resolver as equações diferenciais
[t, x] = ode45(@(t, x) drone_nonlinear_dynamics(t, x, T, tau, m, g, J), tspan, x0);

% GRÁFICO DA ALTURA (Z)
figure;
plot(t, x(:,3), 'LineWidth', 2, 'Color', [0 0.4470 0.7410]);
grid on;
xlabel('Tempo (s)');
ylabel('Altura Z (m)');
title('Simulação 1.5: Resposta do Drone em Altura ');

% GRÁFICO DA VELOCIDADE VERTICAL (Vz)
figure;
plot(t, x(:,6),'LineWidth', 2,'Color', [0.8500 0.3250 0.0980]); 
grid on;
xlabel('Tempo (s)');
ylabel('Velocidade Vz (m/s)');
title('Simulação 1.5: Variação da Velocidade Vertical');
subtitle(['Aceleração Constante: ', num2str((T/m)-g), ' m/s^2']);

subtitle(['Thrust Total: ', num2str(T), ' N | Peso: ', num2str(m*g), ' N']);

% FUNÇÃO DO MODELO NÃO LINEAR
function dxdt = drone_nonlinear_dynamics(~, x, T, tau, m, g, J)
    
    v = x(4:6);
    omega = x(10:12);
    az = (T/m) - g;
    
    % Derivadas do Vetor de Estado
    d_pos = v;                   
    d_vel = [0; 0; az];          
    d_ang = omega;               
    d_omega = inv(J) * (tau - cross(omega, J*omega));
    
    dxdt = [d_pos; d_vel; d_ang; d_omega];

end

%% Matrizes

% Matrizes auxiliares 
O3 = zeros(3);
O3x4 = zeros(3,4);
I3 = eye(3);

% Bloco da gravidade para Hover
A_grav_OP1 = -g * [0, -1, 0; 
                   1,  0, 0; 
                   0,  0, 0];

% Montagem da Matriz A (12x12)
A_OP1 = [ O3,  I3,  O3,          O3 ;
          O3, -D,   A_grav_OP1,  O3 ;
          O3,  O3,  O3,          I3 ;
          O3,  O3,  O3,          O3 ];

% Parâmetros específicos do OP2 (Exemplos)
Vx = 1.0;       
theta = 0.1;    
phi = 0; psi = 0; 

% R(lambda) simplificado e T(lambda) simplificado
R_lam = I3; 
T_lam = I3; 

% Bloco Posição/Velocidade (L1,C3)
A_pos_vel_OP2 = [0, -Vx*sin(theta), 0; 
                 0,  0,             Vx*cos(theta); 
                 0, -Vx*cos(theta), 0];

% Bloco Gravidade (L2,C3)
A_grav_OP2 = -g * [0, -cos(theta), 0; 
                   cos(theta), 0,  0; 
                   0,  sin(theta), 0];

% Montagem da Matriz A_OP2
A_OP2 = [ O3,  R_lam, A_pos_vel_OP2, O3 ;
          O3, -D,     A_grav_OP2,    O3 ;
          O3,  O3,    O3,            T_lam ;
          O3,  O3,    O3,            O3 ];
% Bloco de Thrust (aceleração linear)
B_thrust = [0, 0, 0, 0; 
            0, 0, 0, 0; 
            1/m, 0, 0, 0];

B_thrust2 = [0, 0, 0, 0; 
            0, 0, 0, 0; 
            -1/m, 0, 0, 0];


% Bloco de Momentos (aceleração angular)
Jx = 2.4e-5;
Jy = 2.4e-5;
Jz = 3.2e-5;
B_moments = [0, 1/Jx, 0,    0; 
             0, 0,    1/Jy, 0; 
             0, 0,    0,    1/Jz];

% Montagem da Matriz B (12x4)
B_OP1 = [ O3x4 ; 
      B_thrust ; 
      O3x4 ; 
      B_moments ];

B_OP2 = [ O3x4 ; 
      B_thrust2 ; 
      O3x4 ; 
      B_moments ];
% Valores próprios para Hover
ev_OP1 = eig(A_OP1);

% Valores próprios para Voo Horizontal
ev_OP2 = eig(A_OP2);

% Display
disp('Eigenvalues OP1:'); disp(ev_OP1);
disp('Eigenvalues OP2:'); disp(ev_OP2);

%% CONTROLABILIDADE E OBSERVABILIDADE 
C = eye(12); 

% 1 Matriz de Controlabilidade
Co = ctrb(A_OP1, B_OP1); 
rank_Co = rank(Co);

% 2 Matriz de Observabilidade
Ob = obsv(A_OP1, C);
rank_Ob = rank(Ob);

fprintf('\n--- Análise de Sistemas (1.8) ---\n');
fprintf('Número de Estados (n): 12\n');

% Verificação de Controlabilidade
if rank_Co == 12
    fprintf('Rank Co: %d -> O sistema é TOTALMENTE CONTROLÁVEL.\n', rank_Co);
else
    fprintf('Rank Co: %d -> O sistema NÃO é totalmente controlável (Cuidado!)\n', rank_Co);
end

% Verificação de Observabilidade
if rank_Ob == 12
    fprintf('Rank Ob: %d -> O sistema é TOTALMENTE OBSERVÁVEL.\n', rank_Ob);
else
    fprintf('Rank Ob: %d -> O sistema NÃO é totalmente observável.\n', rank_Ob);
end


% Matriz de Controlabilidade 
Co = ctrb(A_OP2, B_OP2); 
rank_Co2 = rank(Co);

% Matriz de Observabilidade
Ob = obsv(A_OP2, C);
rank_Ob2 = rank(Ob);

fprintf('\n--- Análise de Sistemas (1.8) ---\n');
fprintf('Número de Estados (n): 12\n');

% Verificação de Controlabilidade
if rank_Co2 == 12
    fprintf('Rank Co: %d -> O sistema é TOTALMENTE CONTROLÁVEL2.\n', rank_Co);
else
    fprintf('Rank Co: %d -> O sistema NÃO é totalmente controlável2 (Cuidado!)\n', rank_Co);
end

% Verificação de Observabilidade
if rank_Ob2 == 12
    fprintf('Rank Ob: %d -> O sistema é TOTALMENTE OBSERVÁVEL2.\n', rank_Ob);
else
    fprintf('Rank Ob: %d -> O sistema NÃO é totalmente observável2.\n', rank_Ob);
end

%% EXTRAÇÃO DAS 6 FUNÇÕES DE TRANSFERÊNCIA (1.10) 

% Definir C (identidade para extrair estados) e D
C = eye(12);
D_mat = zeros(12, 4);

% Criar os sistemas 
sys_OP1 = ss(A_OP1, B_OP1, C, D_mat);
G_OP1 = tf(sys_OP1);

% 1 MALHAS INTERNAS (Atitude)
G_nx_phi   = minreal(G_OP1(7, 2)); % Roll
G_ny_theta = minreal(G_OP1(8, 3)); % Pitch
G_nz_psi   = minreal(G_OP1(9, 4)); % Yaw

% 2 MALHAS EXTERNAS (Posição e Altitude)
G_T_pz = minreal(G_OP1(3, 1));

% Translação: 
G_phi_py = minreal(G_OP1(2, 2) / G_OP1(7, 2)); % py/phi
G_theta_px = minreal(G_OP1(1, 3) / G_OP1(8, 3)); % px/theta

%% EXIBIÇÃO FORMATADA
fprintf('\n FUNÇÕES DE TRANSFERÊNCIA\n');

fprintf('\n MALHAS INTERNAS\n');
fprintf('G_nx_phi (phi/nx):'); display(G_nx_phi);
fprintf('G_ny_theta (theta/ny):'); display(G_ny_theta);
fprintf('G_nz_psi (psi/nz):'); display(G_nz_psi);

fprintf('\nMALHAS EXTERNAS\n');
fprintf('G_phi_py (py/phi):'); display(G_phi_py);
fprintf('G_theta_px (px/theta):'); display(G_theta_px);
fprintf('G_T_pz (pz/T):'); display(G_T_pz);

%% ANÁLISE DE ESTABILIDADE EM MALHA FECHADA (Root Locus)

% 2 Criar Figuras para os Root Locus
figure('Name', 'Análise de Estabilidade - Root Locus');

subplot(2,1,1);
rlocus(G_T_pz);
title('Root Locus: Malha de Altitude G_{T,p_z}(s)');
grid on;

subplot(2,1,2);
rlocus(G_phi_py);
title('Root Locus: Malha Lateral G_{\phi,p_y}(s)');
grid on;

%% CÁLCULO DOS POLOS EM MALHA FECHADA 
K = 1;

% Malha fechada com Feedback Unitário: T = (K*G) / (1 + K*G)
sys_CL_alt = feedback(K * G_T_pz, 1);
sys_CL_lat = feedback(K * G_phi_py, 1);

% Extração dos Polos
polos_alt_MA = pole(G_T_pz);       % Malha Aberta
polos_alt_MF = pole(sys_CL_alt);   % Malha Fechada

polos_lat_MA = pole(G_phi_py);     % Malha Aberta
polos_lat_MF = pole(sys_CL_lat);   % Malha Fechada

%% EXIBIÇÃO DOS RESULTADOS

fprintf('DISCUSSÃO DE ESTABILIDADE (Root Locus)\n');

fprintf('\n[ MODO ALTITUDE - G_T_pz ]\n');
fprintf('Polos MA: [%s]\n', num2str(polos_alt_MA'));
fprintf('Polos MF (K=%d): [%s]\n', K, num2str(polos_alt_MF'));

fprintf('\n[ MODO LATERAL - G_phi_py ]\n');
fprintf('Polos MA: [%s]\n', num2str(polos_lat_MA'));
fprintf('Polos MF (K=%d): [%s]\n', K, num2str(polos_lat_MF'));

% Conclusão
fprintf('\n Conclusão da Estabilidade \n');
if any(real(polos_alt_MF) == 0)
    fprintf('Altitude: Marginalmente Estável (Oscilatório puro).\n');
end
if any(real(polos_lat_MF) == 0)
    fprintf('Lateral: Marginalmente Estável (Oscilatório puro).\n');
end
fprintf('Nota: Os polos na origem (1/s^2) exigem um controlador PD/PID para estabilidade estável.\n');