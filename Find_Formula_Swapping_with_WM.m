clear,clc
close all;

% Definition 
syms p r real
assume((0<=r)&(r<=1))
assume((0<=p)&(p<=1))
assume(r, 'real')
assume(p, 'real')

% basic state
e0 = [1 0; 0 sqrt(1-r)];
e1 = [0 sqrt(r); 0 0];
I = eye(2);
I2 = kron(I, I);
H = [1; 0]; % |0>
V = [0; 1]; % |1>

% Pauli operators
sigmax = [0 1; 1 0];
sigmaz = [1 0; 0 -1]; 

% W state
ket_000 = [1; 0; 0; 0; 0; 0; 0; 0];
ket_100 = [0; 0; 0; 0; 1; 0; 0; 0];
ket_010 = [0; 0; 1; 0; 0; 0; 0; 0];
ket_001 = [0; 1; 0; 0; 0; 0; 0; 0];
%Alice  |W>_{123}
W_123 = (1/2) * (ket_100 + ket_010 + sqrt(2) * ket_001);

%Bob  |W>_{456}
W_456 = (1/2) * (ket_100 + ket_010 + sqrt(2) * ket_001);

% density matrix
rho_123 = W_123 * W_123';
rho_456 = W_456 * W_456';

rt0=kron(rho_123,rho_456);

%% joint measurement operators
psi_P1 =1/2*(kron(H, kron(V, H)) + kron(H, kron(H, V)) + sqrt(2)*kron(V, kron(H, H)));
psi_P2 = 1/2*(kron(H, kron(V, H)) + kron(H, kron(H, V)) -sqrt(2)* kron(V, kron(H, H)));
psi_P3 = 1/2*(kron(V, kron(V, H)) + kron(V, kron(H, V)) + sqrt(2)*kron(H, kron(H, H)));
psi_P4 = 1/2*(kron(V, kron(V, H)) + kron(V, kron(H, V)) - sqrt(2)*kron(H, kron(H, H)));

rho_P1 = psi_P1*psi_P1';
rho_P2 = psi_P2*psi_P2';
rho_P3 = psi_P3*psi_P3';
rho_P4 = psi_P4*psi_P4';

% Alice applies measurements on the input state and her qubit of the shared W entangled state
Proj1 = kron(I2,kron(rho_P1 , I));
Proj2 = kron(I2,kron(rho_P2 , I));
Proj3 = kron(I2,kron(rho_P3 , I));
Proj4 = kron(I2,kron(rho_P4 , I));

%% ADC
% Apply ADC to rho_123 & rho_456
e0_3 = kron(kron(I, I), e0);
e1_3 = kron(kron(I, I), e1);
rho_123_prime = e0_3 * rho_123 * e0_3' + e1_3 * rho_123 * e1_3';

e3_0 = kron(kron(e0, e0), I);
e3_1 = kron(kron(e0, e1), I);
e3_2 = kron(kron(e1, e0), I);
e3_3 = kron(kron(e1, e1), I);
rho_456_prime = e3_0 * rho_456 * e3_0' + e3_1 * rho_456 * e3_1' + e3_2 * rho_456 * e3_2' + e3_3 * rho_456 * e3_3';

disp('Alice after ADC ｜density matrix:');
disp(simplify(rho_123_prime));
A=(1/2) * (ket_100 + ket_010 + sqrt(2*(1-r)) * ket_001);
display('Formula verify | Alice after ADC:');
assume((0<=r)&(r<=1))
disp(simplify(rho_123_prime-A*A'-(r/2)*ket_000*ket_000'));

disp('Bob after ADC｜density matrix:');
disp(simplify(rho_456_prime));
B=(1/2) * (sqrt(1-r)*ket_100 + sqrt(1-r)*ket_010 + sqrt(2) * ket_001);
display('Formula verify | Bob after ADC:');
assume((0<=r)&(r<=1))
disp(simplify(rho_456_prime-B*B'-(r/2)*ket_000*ket_000'));


%% PartialTrace and Swap
%%Initial shared state for comparing
R0=(PartialTrace(Proj1*rt0*Proj1',[3,4,5]));
R0=(R0/trace(R0))

%% Shared state after ADc
rt=kron(rho_123_prime,rho_456_prime);
R1e=(PartialTrace(Proj1*rt*Proj1',[3,4,5]));

display('Density matrix | damped state:')
assume((0<=r)&(r<=1))
R1e=simplify(R1e/(trace(R1e)));
disp(R1e);
C=(1/(2*sqrt(1+r))) * (ket_100 + ket_010 + sqrt(2) * ket_001);
assume((0<=r)&(r<=1))
display('Formula verify | damped state:')
disp(simplify(R1e-C*C'-r/(1+r)*ket_000*ket_000'))

%% Weak measurement  on shared state
% Apply WM on each qubit
M = [sqrt(1-p), 0; 0, 1];
M3=kron(kron(M,M),M);
display('Density matrix | final shared state:')
assume((0<=r)&(r<=1))
assume((0<=p)&(p<=1))
Rf=simplify(M3*R1e*M3');
P=simplify(trace(Rf));
Rf=Rf/P
D=(1/2) * (ket_100 + ket_010 + sqrt(2) * ket_001);
trace(Rf);
display('Formula verify | final shared state:')
disp(simplify(Rf-(1/P)*((1-p)^2/(1+r))*(D*D'+r*(1-p)*ket_000*ket_000')))


%% Fidelity and Success probility
assume((0<=r)&(r<=1))
assume((0<=p)&(p<=1))
% success probability
display('P | final shared state:')
P
% Fidelity
display('Fid | final shared state with WM:')
fidelity = simplify(trace(Rf*R0))
display('Fid | damped state:')
simplify(trace(R1e*R0))