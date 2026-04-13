% Unidade Curricular: UAV's
% Project 1 - Crazyflie drone modeling and identification

close;clear all; clc;

% Ler ficheiro
csvfilename = '2024-04-04_log07.csv';
array = dlmread(csvfilename,',',1,0);

tempo = array(:,1)'*1e-3;
pos = array(:,2:4)';
vel = array(:,5:7)';
lbd = array(:,8:10)'*pi/180;
om = array(:,11:13)'*pi/180;
pos_ref = array(:,14:16)';
yaw_ref = array(:,17)';
motores = array(:,18:21)';

% Converter informação
t = tempo - tempo(1);
x = [pos;vel;lbd;om];
x_ref = [pos_ref;0*vel;lbd*0;om*0];
x_ref(9,:) = yaw_ref;
uint16_max = 2^16;
u = motores/uint16_max;

% plot dados
initPlots;
vehicle3d_ref_show_data(t,x,u,x_ref);

% preparar dados para ID
u_id = lbd(1,:)';
y_id = pos(2,:)';

%% 3 IDENTIFICAÇÃO AUTOMÁTICA E STEP RESPONSE

% 1 Configuração de Dados e Amostragem
Ts = mean(diff(t)); 

data_px = iddata(pos(1,:)', lbd(2,:)', Ts);
data_py = iddata(pos(2,:)', lbd(1,:)', Ts);
data_pz = iddata(pos(3,:)', u(1,:)', Ts);

% 2 Funções de Transferência 
np = 2; 
nz = 0;

disp('A calcular funções de transferência ...');
sys_px = tfest(data_px, np, nz);
sys_py = tfest(data_py, np, nz);
sys_pz = tfest(data_pz, np, nz);

%% 3. PLOT DAS STEP RESPONSES COM EIXOS ADEQUADOS

figure('Name', 'Identificação Experimental - Step Responses', 'Color', [1 1 1]);

% Eixo X (Px)
subplot(3,1,1);
t_final_x = 2500; 
step(sys_px, t_final_x);
grid on; grid minor;
title('Figura 17: Step response da função G_{\theta,p_x}(s)');
ylabel('Amplitude (p_x)');
xlim([0 t_final_x]); 

% Eixo Y (Py)
subplot(3,1,2);
t_final_y = 300; 
step(sys_py, t_final_y);
grid on; grid minor;
title('Figura 18: Step response da função G_{\phi,p_y}(s)');
ylabel('Amplitude (p_y)');
xlim([0 t_final_y]);

% Eixo Z (Pz)
subplot(3,1,3);
t_final_z = 14; 
step(sys_pz, t_final_z);
grid on; grid minor;
title('Figura 19: Step response da função G_{T,p_z}(s)');
ylabel('Amplitude (p_z)');
xlabel('Time (seconds)');
xlim([0 t_final_z]);

%% 4 RESULTADOS
fprintf('\n Funções de Transferência Calculadas\n');
fprintf('Modo Px:'); display(tf(sys_px));
fprintf('Modo Py:'); display(tf(sys_py));
fprintf('Modo Pz:'); display(tf(sys_pz));

%% 5 COMPARAÇÃO DE OUTPUTS MODELO COM EXPERIMENTAL
figure('Name', 'Comparação: Medido vs Simulado', 'Color', [1 1 1]);

% Comparação Eixo X (G_theta,px)
subplot(3,1,1);
compare(data_px, sys_px); 
grid on;
title('Figura 20: Comparação de G_{\theta,p_x}(s)');
ylabel('p_x [m]');

% Comparação Eixo Y (G_phi,py)
subplot(3,1,2);
compare(data_py, sys_py);
grid on;
title('Figura 21: Comparação de G_{\phi,p_y}(s)');
ylabel('p_y [m]');

% Comparação Eixo Z (G_T,pz)
subplot(3,1,3);
compare(data_pz, sys_pz);
grid on;
title('Figura 22: Comparação de G_{T,p_z}(s)');
ylabel('p_z [m]');
xlabel('Time (seconds)');

%% 5 ANÁLISE DE QUALIDADE DO AJUSTE

[~, fit_px] = compare(data_px, sys_px);
[~, fit_py] = compare(data_py, sys_py);
[~, fit_pz] = compare(data_pz, sys_pz);

fprintf('\n Qualidade do Ajuste\n');
fprintf('Fit para Px: %.2f%%\n', fit_px);
fprintf('Fit para Py: %.2f%%\n', fit_py);
fprintf('Fit para Pz: %.2f%%\n', fit_pz);

%% 6 IDENTIFICAÇÃO DE ORDEM SUPERIOR (eixo X)
% Configuração: 6 polos e 5 zeros
np_high = 6; 
nz_high = 5;

fprintf('\n--- Estimando Modelo de 6ª Ordem para Eixo X ---\n');
sys_px_6p5z = tfest(data_px, np_high, nz_high);

% COMPARAÇÃO DIRETA DE PERFORMANCE
figure('Name', 'Comparação: 2 Polos vs 6 Polos (Eixo X)', 'Color', [1 1 1]);

[~, fit_2p] = compare(data_px, sys_px);
[~, fit_6p] = compare(data_px, sys_px_6p5z);

% Plot de comparação
compare(data_px, sys_px, sys_px_6p5z);
grid on;
legend('Dados Experimentais', ...
       ['Original (2 Polos) - Fit: ' num2str(fit_2p, '%.2f') '%'], ...
       ['Superior (6 Polos, 5 Zeros) - Fit: ' num2str(fit_6p, '%.2f') '%'], ...
       'Location', 'best');

title('Figura 23: Impacto da Ordem do Modelo na Fidelidade dos Dados (Eixo X)');
ylabel('Posição p_x [m]');

% RESULTADOS
fprintf('\nResultados Comparativos\n');
fprintf('Ajuste Modelo 2ª Ordem: %.2f%%\n', fit_2p);
fprintf('Ajuste Modelo 6ª Ordem: %.2f%%\n', fit_6p);
fprintf('Melhoria na Identificação: %.2f%%\n', fit_6p - fit_2p);

fprintf('\nNova Função de Transferência (Px - 6ª Ordem):\n');
display(tf(sys_px_6p5z));
