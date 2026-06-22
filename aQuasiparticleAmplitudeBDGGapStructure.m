% Parameters
t = 1.0; % Hopping parameter
mu = 0.0; % Chemical potential
Delta = 1.0; % Superconducting gap parameter
Nk = 100; % Number of k points

% Define k-space grid
kx = linspace(-pi, pi, Nk);

% Initialize arrays for eigenvalues and colors
eigenvalues_1d = zeros(Nk, 2);
colors = zeros(Nk, 3); % RGB colors

for i = 1:Nk
    % Hamiltonian for 1D
    e_k = -t * cos(kx(i)) - mu;
    H_BdG = [e_k, Delta; Delta, -e_k];
    
    % Diagonalize the Hamiltonian
    [e_vecs, e_vals] = eig(H_BdG);
    eigenvalues = diag(e_vals)';
    eigenvalues_1d(i, :) = eigenvalues;
    
    % Amplitudes for color mixing
    electron_amplitude = abs(e_vecs(1, :)).^2;
    hole_amplitude = abs(e_vecs(2, :)).^2;
    
    % Normalizing amplitudes for color scaling
    total_amplitude = max(electron_amplitude + hole_amplitude);
    colors(i, 1) = electron_amplitude(1) / total_amplitude; % Red
    colors(i, 3) = hole_amplitude(2) / total_amplitude; % Blue
end

% Plotting
figure;
hold on;
for i = 1:Nk
    scatter(kx(i), eigenvalues_1d(i, 1), [], [colors(i, 1), 0, colors(i, 3)], 'filled');
    scatter(kx(i), eigenvalues_1d(i, 2), [], [colors(i, 1), 0, colors(i, 3)], 'filled');
end
title('Quasiparticle Energy Bands with Color Mixing');
xlabel('$k_x$', 'Interpreter', 'latex');
ylabel('Energy', 'Interpreter', 'latex');
line(xlim, [0 0], 'Color', 'black', 'LineStyle', '--'); % Fermi level for reference
hold off;

%% Zeeman splitting

% Parameters
t = 1.0; % Hopping parameter
mu = 0.0; % Chemical potential
Delta = 1.0; % Superconducting gap parameter
V_Z = 0.5; % Maximum Zeeman energy
Nk = 200; % Number of k points

% Define k-space grid
kx = linspace(-pi, pi, Nk+1);kx(end)=[];
ky=linspace(-pi,pi,Nk+1);ky(end)=[];

% Initialize arrays for eigenvalues and eigenvectors
eigenvalues_4by4 = zeros(Nk,Nk, 4);
eigenvectors_4by4 = zeros(4, 4, Nk,Nk);

for i = 1:Nk
    for j=1:Nk
    % Momentum-dependent Zeeman splitting
    V_Z_kx = V_Z * sin(kx(i))*sin(ky(j));
    
    % Hamiltonian components for electron and hole
    e_k = -t * (cos(kx(i))+cos(ky(j))) - mu;
    h_k = -e_k; % Assuming symmetry for simplicity, adjust as needed
    
    % 4x4 BdG Hamiltonian with momentum-dependent Zeeman splitting
    H_BdG = [e_k+V_Z_kx, 0, 0, Delta;
             0, e_k-V_Z_kx, -Delta, 0;
             0, -Delta, h_k-V_Z_kx, 0;
             Delta, 0, 0, h_k+V_Z_kx];
    
    % Diagonalize the Hamiltonian
    [e_vecs, e_vals] = eig(H_BdG);
    eigenvalues = diag(e_vals)';
    eigenvalues_4by4(i,j, :) = eigenvalues;
    eigenvectors_4by4(:, :, i,j) = e_vecs;
    end
end

% Plotting with color mixing
figure;
ax1=gca;
hold on;
color=zeros(Nk,Nk,4,3);
for i = 1:Nk
    for j=1:Nk
        for band = 1:4
            % Extracting the amplitudes for color mixing
            amplitude = abs(eigenvectors_4by4(:, band, i,j)).^2;
            color(i,j,band,:) = [amplitude(1) + amplitude(2), 0, amplitude(3) + amplitude(4)];
            color = color / max(max(max(max(color)))); % Normalize color intensity
            % scatter(kx(i), eigenvalues_4by4(i, band), 36, color, 'filled');
        end
    end
end

title('Quasiparticle Energy Bands with Color Mixing');
xlabel('$k_x$', 'Interpreter', 'latex');
ylabel('Energy', 'Interpreter', 'latex');
hold off;
%% 

%%%%
figure
ax1=gca;
plot([EEE1(101:200,101);EEE1(1,101:200)';flip(diag(EEE1(101:200,101:200)));diag(EEE1(101:-1:2,101:200));EEE1(1:101,1);EEE1(101,200:-1:101)'],'r','LineWidth',0.8)
hold on
plot([EEE2(101:200,101);EEE2(1,101:200)';flip(diag(EEE2(101:200,101:200)));diag(EEE2(101:-1:2,101:200));EEE2(1:101,1);EEE2(101,200:-1:101)'],'b','LineWidth',0.8)
plot([EEE3(101:200,101);EEE3(1,101:200)';flip(diag(EEE3(101:200,101:200)));diag(EEE3(101:-1:2,101:200));EEE3(1:101,1);EEE3(101,200:-1:101)'],'r','LineWidth',0.8)
plot([EEE4(101:200,101);EEE4(1,101:200)';flip(diag(EEE4(101:200,101:200)));diag(EEE4(101:-1:2,101:200));EEE4(1:101,1);EEE4(101,200:-1:101)'],'b','LineWidth',0.8)
xticks([0 N/2 N 3*N/2 2*N 5*N/2 3*N])
xticklabels({'\Gamma', 'X', 'M','\Gamma','M', 'Y', '\Gamma'})
xlim([0 3*N])
ax1.XGrid='on'
ax1.GridAlpha=1;
ax1.FontName='Times New Roman';
ax1.FontSize=12;
