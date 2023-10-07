%% Clear Everything
close all;clear all;clc;
warning('off','all')

%% Create & Plot System
n=2;
x(:,1)=ones(1,n);
% x(:,1)=[0.5; 0.8];

t(1)=0;
dt=1e-4;

% tf=5;tfcc=20;
tf=5;
nA=5;
c=5;
epsil=0.05;
M=tf/dt;
t=dt:dt:tf;

for i=1:M
x(1,i+1) = x(1,i)+cos (x(1,i+1)
x(2,i+1) =cos x(1,I+1)

% figure;
% plot(t,x(1,1:M),'k');hold on;
% plot(t,x(2,1:M),'r');
% legend('x_1', 'x_2')

%% Find the Koopman operator matrix using a dictionary of observables


Sai_of_u=@(u) [u u^2 u^3 ];
Sai_nabla=@(u) [1 2*u 3*u^2];

X1=x(:,1:M);
X2=(x(:,2:M+1)-X1)/dt;

% ----EDMD-------------
% Observable Set 1
% obs_size = 4;
% Sai=@(xx) [xx(1) xx(2) (-xx(1)+xx(2)^2) (-2*xx(2)+xx(1)*xx(2))];
% Sai_u=@(xx) [1 0 1 xx(2)];
% Sai_z=@(xx) [0 1 2*xx(2) xx(1)];

% Observable Set 2
% obs_size = 4;
% Sai=@(xx) [xx(1) xx(2) xx(2)^2 xx(1)*xx(2)];
% Sai_u=@(xx) [1 0 1 xx(2)];
% Sai_z=@(xx) [0 1 2*xx(2) xx(1)];

% Observable Set 3
% obs_size = 6;
% Sai=@(xx) [xx(1) xx(2) xx(2)^2 xx(1)*xx(2) (-4*xx(2)^2 + 2*xx(1)*xx(2)^2) (3*xx(1)*xx(2) + xx(1)^2 * xx(2) + xx(2)^2)];
% Sai_u=@(xx) [1 0 1 xx(2) 2*xx(2)^2 (3*xx(2) + 2*xx(1)*xx(2))];
% Sai_z=@(xx) [0 1 2*xx(2) xx(1) (-8*xx(2) + 4*xx(1)*xx(2)) (3*xx(1) + xx(1)^2 + 3*xx(2)^2)];

% Observable Set 4
% obs_size = 6;
% Sai=@(xx) [xx(1) xx(2) (xx(2)^2) (xx(1)*xx(2)) (- 3*xx(2)^2 + 2*xx(1)*xx(2)^2) (- 5*xx(1)*xx(2) + xx(1)^2 * xx(2) + xx(2)^3)];
% Sai_u=@(xx) [1 0 1 xx(2) (2*xx(2)^2) (-5*xx(2) + 2*xx(1)*xx(2))];
% Sai_z=@(xx) [0 1 2*xx(2) xx(1) (-6*xx(2) + 4*xx(1)*xx(2)) (-5*xx(1) + xx(1)^2 + 3*xx(2)^2)];

% Ideal Observalbe Set
% obs_size = 8;
% Sai=@(xx) [xx(1) xx(2) xx(2)^2 xx(1)*xx(2) xx(2)^3 ...
%     xx(1)^2*xx(2) xx(1)*xx(2)^2 xx(1)*xx(2)^3];
% Sai_u=@(xx) [1 0 0 xx(2) 0 2*xx(1)*xx(2) xx(2)^2 xx(2)^3];
% Sai_z=@(xx) [0 1 2*xx(2) xx(1) 3*xx(2)^2 xx(1)^2  2*xx(1)*xx(2) 3*xx(1)*xx(2)^2];

% Observalbe Set 6
% obs_size = 7;
% Sai=@(xx) [xx(1) xx(2) xx(2)^2 xx(1)*xx(2) xx(2)^3 ...
%     xx(1)^2*xx(2) xx(1)*xx(2)^2];
% Sai_u=@(xx) [1 0 0 xx(2) 0 2*xx(1)*xx(2) xx(2)^2];
% Sai_z=@(xx) [0 1 2*xx(2) xx(1) 3*xx(2)^2 xx(1)^2  2*xx(1)*xx(2)];

% Observalbe Set 7
obs_size = 8;
Sai=@(xx) [xx(1) xx(2) xx(2)^2 xx(1)*xx(2) xx(2)^3 ...
    xx(1)^2*xx(2) xx(1)*xx(2)^2 xx(1)^3*xx(2)];
Sai_u=@(xx) [1 0 0 xx(2) 0 2*xx(1)*xx(2) xx(2)^2 3*xx(1)^2*xx(2)];
Sai_z=@(xx) [0 1 2*xx(2) xx(1) 3*xx(2)^2 xx(1)^2  2*xx(1)*xx(2) xx(1)^3];

% Observalbe Set 8
% obs_size = 8;
% Sai=@(xx) [xx(1) xx(2) xx(2)^2 xx(1)*xx(2) xx(2)^3 ...
%     xx(1)^2*xx(2) xx(1)*xx(2)^2 xx(1)^3];
% Sai_u=@(xx) [1 0 0 xx(2) 0 2*xx(1)*xx(2) xx(2)^2 3*xx(1)^2];
% Sai_z=@(xx) [0 1 2*xx(2) xx(1) 3*xx(2)^2 xx(1)^2  2*xx(1)*xx(2) 0];

q1=tic;
% SaiX = zeros(n, M);
for i=1:M
        SaiX(i,:)=Sai(X1(:,i));
%         M1(i,:)=X2(2:3,i)'*Sai_z(X1(:,i));
%         M2(i,:)=X2(1,i)'*Sai_u(X1(:,i));
        Gamma(i,:)=X2(:,i)'*[Sai_u(X1(:,i));Sai_z(X1(:,i))];
end
Kedmd=pinv(SaiX)*Gamma;
elapsed_EDMD=toc(q1)

%% Plot Poles and Eigenvalues

[Vedmd,Lambda_edmd]=eig(Kedmd);
for i=1:size(Lambda_edmd,1)
    poles_edmd(i)=Lambda_edmd(i,i);
end

% figure;
% plot(real(poles_edmd),imag(poles_edmd),'b*','linewidth',1);hold on;
% legend('EDMD');
% grid

%% EDMD Reconstruction

B=[1 zeros(1,obs_size-1);0 1 zeros(1,obs_size-2)];
Phi_edmd=@(x) Sai(x)*Vedmd;

% Koopman Modes
V=(Vedmd\B')';
%Phi0=pinv(V)*X1(:,1);
zedmd(:,1)=X1(:,1);

for k=1:M
    zedmd(:,k+1)=V*(expm(Lambda_edmd.*dt*(k+1)))*Phi_edmd(zedmd(:,1))';
end

figure('Position', [100, 100, 1100, 400]);
subplot(1, 2, 1)
plot(t,x(1,1:M),'k:','linewidth',2.5);hold on;
plot(t,zedmd(1,1:M),'b--','linewidth',1.5 );hold on;
plot(t,x(1,1:M)- zedmd(1,1:M), 'g', 'linewidth',1.5); hold on;
legend('Actual System','EDMD Approximation', 'Error');
xlabel('t')
ylabel('x_1')

subplot(1, 2, 2)
plot(t,x(2,1:M),'k:','linewidth',2.5);hold on;
plot(t,zedmd(2,1:M),'b--','linewidth',1.5 );hold on;
plot(t,x(2,1:M)- zedmd(2,1:M), 'g','linewidth',1.5); hold on;
legend('Actual System','EDMD Approximation', 'Error');
xlabel('t')
ylabel('x_2')

%% ROA Approximation
A = Kedmd;
A = Lambda_edmd;
Q = eye(obs_size);
P = lyap(A', Q);
% Lyapunov Function
x_sym = sym('x', [n 1], 'real');
phi_x = Phi_edmd(x_sym).'
V_x = ctranspose(phi_x) * P * phi_x

rad = 3;
[x1_mesh x2_mesh] = meshgrid(-1*rad:0.1:rad);
x1_data = linspace(-1*rad, 1*rad, 61);
x2_data = linspace(-1*rad, 1*rad, 61);

Vd1 = matlabFunction(diff(V_x, x_sym(1)))
Vd2 = matlabFunction(diff(V_x, x_sym(2)))

[x1d, x2d] = f_x_mat(x1_mesh, x2_mesh);
V_d = Vd1(x1_mesh, x2_mesh).*x1d + Vd2(x1_mesh, x2_mesh).*x2d;

c = max_c(x1_mesh, x2_mesh, V_x, V_d)
plot_ROA(@(t,xx) f_sys1(xx), 10)
hold on

fimplicit(V_x == c, '-r', 'LineWidth', 1.5)

function [x1d, x2d] = f_x_mat(x1, x2)
    x1d = -x1 + x2.^2;
    x2d = -2*x2 + x1.*x2;
end




