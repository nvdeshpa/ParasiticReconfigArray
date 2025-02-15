%% MATLAB simulation code for generating results from the paper "Beamforming with hybrid reconfigurable parasitic antenna arrays"
%% Authors: Nitish Vikas Deshpande, Miguel R. Castellanos, Saeed R. Khosravirad, Jinfeng Du, Harish Viswanathan, and Robert W. Heath Jr.
clear all
close all
rng default % For reproducibility

%% Universal constants 
c=3*10^8;
%% Fixed simulation parameters
f_c=7*10^9;
lambda_c= c/f_c;
l=0.5*lambda_c;
Z_0=50;  %% reference impedance
sigma_c= 5.7*10^7;
a=lambda_c/500; %% Radius of dipole antenna
k=2*pi/lambda_c; %% wavenumber 
Rzr_reference=0.05; %% Real part of parasitic impedance is assumed to be small to avoid losses
dx_by_lambda_c=0.4; %% antenna spacing along x axis
dy_by_lambda_c=0.5; %% antenna spacing along y axis
%% Parameters for path loss
r= 250; %% Distance between transmitter and receiver
Kb=1.380649*10^(-23); %% Boltzmann's constant
TA=300; %% Antenna noise temperature in K
BW=20*10^6; %% Communication signal bandwidth
G=(lambda_c/(4*pi*r))^2; %% free space path loss

%% Parameters for power consumption: ref Hybrid MIMO Architectures for Millimeter Wave Communications: Phase Shifters or Switches?
P_rf= 240*10^(-3); %% Power consumed by RF chain
P_ps= 30*10^(-3); %% phase shifter


%% Generate impedance matrix from MATLAB closed-form expressions or use the impedance matrices generated from Feko simulations
impedance_mat_from_matlab=1;
impedance_mat_from_feko=0;

%% Number of Monte Carlo iterations used for random channel generation
num_realizations=5000;

%% Vary transmit power in dBm
vary_transmit_power=1;
if vary_transmit_power==1
    %P_t_arr_dBm= [-10:1:30];
    %P_t_arr_dBm= [10:0.1:20];
    P_t_arr_dBm= [10,11,12];
else
    P_t_arr_dBm= [30];
end
%% Power conversion from dBm to linear scale
P_t_arr=10.^((P_t_arr_dBm-30)./10);


%% Vary number of parasitic elements and active antennas
vary_number_of_parasitics=1;
%% number of parasitic elements per active antenna
if vary_number_of_parasitics==1
    Np_arr = [1,2,3,4,5,6];
else
    Np_arr = 2;
end

vary_number_of_active_antennas=0;
if vary_number_of_active_antennas==1
    Na_arr = [4:2:32];
else
    Na_arr = [4];
    %Na_arr = [6];
end


%% Store SNR, power consumed, spectral and energy efficiency results 
SNR_MAMP_arr=zeros(num_realizations, length(P_t_arr_dBm), length(Na_arr), length(Np_arr));
SE_MAMP_arr=zeros(num_realizations, length(P_t_arr_dBm), length(Na_arr), length(Np_arr));
SNR_only_active_arr=zeros(num_realizations, length(P_t_arr_dBm), length(Na_arr), length(Np_arr)); %% SNR with all parasitic elements removed (FA-ULA)
SE_only_active_arr=zeros(num_realizations, length(P_t_arr_dBm), length(Na_arr), length(Np_arr));
SNR_fully_active_arr=zeros(num_realizations, length(P_t_arr_dBm), length(Na_arr), length(Np_arr)); %% SNR with fully active elements (FA-UPA) 
SE_fully_active_arr=zeros(num_realizations, length(P_t_arr_dBm), length(Na_arr), length(Np_arr));

SNR_hybrid_ps_arr=zeros(num_realizations, length(P_t_arr_dBm), length(Na_arr), length(Np_arr)); %% SNR with sub-connected hybrid array with Phase shifters (HPS-UPA)
SE_hybrid_ps_arr=zeros(num_realizations, length(P_t_arr_dBm), length(Na_arr), length(Np_arr));

p_consumed_MAMP_arr=zeros(length(P_t_arr),length(Na_arr), length(Np_arr) );
p_consumed_only_active_arr=zeros(length(P_t_arr),length(Na_arr), length(Np_arr) );
p_consumed_fully_active_arr=zeros(length(P_t_arr),length(Na_arr), length(Np_arr) );
p_consumed_hybrid_PS_arr=zeros(length(P_t_arr),length(Na_arr), length(Np_arr) );


EE_MAMP_arr=zeros(num_realizations, length(P_t_arr_dBm), length(Na_arr), length(Np_arr));
EE_only_active_arr=zeros(num_realizations, length(P_t_arr_dBm), length(Na_arr), length(Np_arr));
EE_fully_active_arr=zeros(num_realizations, length(P_t_arr_dBm), length(Na_arr), length(Np_arr));
EE_hybrid_PS_arr=zeros(num_realizations, length(P_t_arr_dBm), length(Na_arr), length(Np_arr));


for np_i = 1:length(Np_arr)
    for na_i=1:length(Na_arr)
        Na=Na_arr(na_i);
        Np=Np_arr(np_i);
        fprintf("Number of active antennas is %f\n", Na);
        fprintf("Number of parasitic elements per active antenna is %f\n", Np);

        dx=dx_by_lambda_c*lambda_c;
        dy=dy_by_lambda_c*lambda_c;
        Z_s= self_impedance_dip(k, l, a);

        Nx=Np+1;
        Ny=Na;
        N=(Np+1)*Na;
        if mod(Np+1,2)==0
            index_active_antenna= [(Np+1)/2: Np+1:(Np+1)/2+  (Np+1)*(Na-1)];
        else
            index_active_antenna= [(Np+2)/2:  Np+1: (Np+2)/2+(Np+1)*(Na-1)] ;    
        end
        %% Generate antenna impedance matrices
        if impedance_mat_from_matlab==1
            fprintf("Using impedance matrices from closed-form expressions implemented in MATLAB")
            Z_ref0= zeros(Nx, Ny); %% Impedances of all antennas wrt the 1st antenna

            for nx=1:Nx
                for ny=1:Ny
                    if nx==1 && ny ==1 
                        %% Self impedance 
                        Z_ref0(nx, ny)=  self_impedance_dip(k, l, a);
                    else
                        d_upa= sqrt(dx^2.*((nx-1)^2)+dy^2.*((ny-1)^2));
                        Z_ref0(nx, ny)= MI_sbs_dip(d_upa/lambda_c,l/lambda_c);
                    end
                end
            end

            Z_tx=zeros(N,N);
            for n1=1:N
                for n2=1:N
                    n1x= 1+mod(n1-1, Np+1);
                    n1y=ceil(n1/(Np+1));
                    n2x=1+mod(n2-1, Np+1);
                    n2y=ceil(n2/(Np+1));
                    Z_tx(n1, n2)= Z_ref0(1+abs(n1x-n2x),1+abs(n1y-n2y));
                end
            end

            R_loss=  compute_loss_resistance(k, l, a, f_c, sigma_c);
            Z_tx=Z_tx+R_loss*eye(Na*(Np+1));

            [Z_P, Z_M, Z_A, parasitic_antenna_indices]= compute_Zp_Zo_from_Zmat(Z_tx, Na*(Np+1), index_active_antenna );
        elseif impedance_mat_from_feko==1
            fprintf("Using impedance matrices from Feko simulations");
            if Na==4
                if Np==1
                    Z_matrix=load('Feko_data/Z_matrix_Na4Np1.mat');
                elseif Np==2
                    Z_matrix=load('Feko_data/Z_matrix_Na4Np2.mat');
                elseif Np==3
                    Z_matrix=load('Feko_data/Z_matrix_Na4Np3.mat');
                elseif Np==4
                    Z_matrix=load('Feko_data/Z_matrix_Na4Np4.mat');
                end
                    
            elseif Na==6
                if Np==2
                    Z_matrix=load('Feko_data/Z_matrix_Na6Np2.mat');
                end
            end
            Z_tx=Z_matrix.Z_matrix;
            Z_A=Z_tx(1:Na, 1:Na);
            Z_P=Z_tx(Na+1:N, Na+1:N);
            Z_M=Z_tx(Na+1:N, 1:Na);
            parasitic_antenna_indices=setdiff([1:1:N], index_active_antenna);
        end


        sigma_sq= 4*Kb*TA*BW/real(Z_tx(1,1)); %% noise variance 

        %% Generate channel realizations
        for ti=1:num_realizations

            %% Multipath channel params

            L=4;
            channel_gain_var=1;
            theta_min=-pi;
            theta_max=pi;
            [h, hp, ha]= generate_multipath_channel(L, channel_gain_var, theta_min, theta_max, N, parasitic_antenna_indices, index_active_antenna, f_c, dx, dy,Na, Np,  c);

            %% decoupled closed-form solution
            Izr_decoupled=zeros(Np, Ny);
            for ny=1:Ny
                hp_j=hp((ny-1)*Np+1:1:ny*Np);
                zm_jj= Z_M((ny-1)*Np+1:1:ny*Np, ny);
                for np=1:Np
                   zeta=1/(real(Z_P(1,1))+Rzr_reference);
                   phi_opt= -angle(1-1/ha(ny)*(   zeta/2*hp_j.'*(zm_jj) )) + angle(1/ha(ny)*hp_j(np)*zm_jj(np)  );
                   Izr_decoupled(np, ny)=  -imag(Z_P(1,1))- cot(0.5*(phi_opt))/zeta;
                end
            end
            Z_R_star= diag(ones(Np*Na,1)*Rzr_reference+ 1j*Izr_decoupled(:));
            Z_eff_star=real(Z_A)-(Z_M'*(inv(Z_P+Z_R_star))')*real(Z_M) - real(Z_M.')*inv(Z_P+Z_R_star)*Z_M + (Z_M'*(inv(Z_P+Z_R_star))'*real(Z_P)*inv(Z_P+Z_R_star)*Z_M);    %% Effective impedance matrix
            h_eff_star=ha-Z_M.'*inv(Z_P+Z_R_star)*hp;   %% Effective channel

            Z_fully_active=[Z_A  Z_M.';  Z_M  Z_P];   %% Fully active impedance matrix
            %%  Solution for sub-connected hybrid phase-shifter

            f_PS_vec= exp(1j*angle(conj(h)));
            F_PS_subconnected=  construct_subconnected_F_PS_matrix(f_PS_vec, Na, Np);
            h_eff_ps=(F_PS_subconnected*h);

            Z_eff_ps=conj(F_PS_subconnected)*real(Z_fully_active)*F_PS_subconnected.';

            %% Vary transmit power

            for p_i=1:length(P_t_arr)
                P= P_t_arr(p_i);

                i_A_star=sqrt(P)*inv(Z_eff_star)*conj(h_eff_star)/norm(inv(sqrtm(Z_eff_star))*conj(h_eff_star));

                
                SNR_MAMP_arr(ti, p_i, na_i, np_i)= G/sigma_sq*(i_A_star.'*h_eff_star)^2;
                SE_MAMP_arr(ti, p_i, na_i, np_i)=log2(1+SNR_MAMP_arr(ti, p_i, na_i, np_i));
                SNR_only_active_arr(ti, p_i, na_i, np_i)=G/sigma_sq*P*ha'*inv(real(Z_A))*ha;
                SE_only_active_arr(ti, p_i, na_i, np_i)=log2(1+SNR_only_active_arr(ti, p_i, na_i, np_i));

                SNR_hybrid_ps_arr(ti, p_i, na_i, np_i)=  G/sigma_sq*P*h_eff_ps'*inv(Z_eff_ps)*h_eff_ps;
                SE_hybrid_ps_arr(ti, p_i, na_i, np_i)= log2(1+SNR_hybrid_ps_arr(ti, p_i, na_i, np_i));

                SNR_fully_active_arr(ti, p_i, na_i, np_i)=G/sigma_sq*P*h'*inv(real(Z_fully_active))*h;
                SE_fully_active_arr(ti, p_i, na_i, np_i)=log2(1+SNR_fully_active_arr(ti, p_i, na_i, np_i));
                p_consumed_MAMP_arr(p_i, na_i, np_i)=P+Na*P_rf;
                p_consumed_only_active_arr(p_i, na_i, np_i)=P+Na*P_rf;
                p_consumed_fully_active_arr(p_i, na_i, np_i)=P+(Na+Na*Np)*P_rf;
                p_consumed_hybrid_PS_arr(p_i, na_i, np_i)= P+ Na*P_rf + (Na+Na*Np)*P_ps  ;

                EE_MAMP_arr(ti, p_i, na_i, np_i)=SE_MAMP_arr(ti, p_i, na_i, np_i)/p_consumed_MAMP_arr(p_i, na_i, np_i);
                EE_only_active_arr(ti, p_i, na_i, np_i)=SE_only_active_arr(ti, p_i, na_i, np_i)/p_consumed_only_active_arr(p_i, na_i, np_i);
                EE_fully_active_arr(ti, p_i, na_i, np_i)=SE_fully_active_arr(ti, p_i, na_i, np_i)/p_consumed_fully_active_arr(p_i, na_i, np_i);
                EE_hybrid_PS_arr(ti, p_i, na_i, np_i)= SE_hybrid_ps_arr(ti, p_i, na_i, np_i)/p_consumed_hybrid_PS_arr(p_i, na_i, np_i);
            end
        end

    end
end

SE_MAMP_averaged= mean(SE_MAMP_arr, 1);
SE_only_active_averaged= mean(SE_only_active_arr, 1);
SE_fully_active_averaged= mean(SE_fully_active_arr, 1);
SE_hybrid_ps_arr_averaged=mean(SE_hybrid_ps_arr,1);

EE_MAMP_averaged= mean(EE_MAMP_arr, 1);
EE_only_active_averaged= mean(EE_only_active_arr, 1);
EE_fully_active_averaged= mean(EE_fully_active_arr, 1);
EE_hybrid_PS_averaged=mean(EE_hybrid_PS_arr,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generating plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if vary_transmit_power==1 && vary_number_of_parasitics==0 && vary_number_of_active_antennas ==0
    %% Generate plots with transmit power on x axis keeping remaining parameters fixed
    
    figure
    
    %% Plot spectral efficiency vs transmit power
    
    plot(P_t_arr_dBm, SE_MAMP_averaged, '-x', 'Markersize', 10, 'LineWidth', 3, 'DisplayName', 'HRP-UPA (Proposed)')
    hold on
    plot(P_t_arr_dBm, SE_hybrid_ps_arr_averaged, '-v', 'Markersize', 10, 'LineWidth', 3, 'DisplayName', 'HPS-UPA')
    plot(P_t_arr_dBm, SE_fully_active_averaged, '-o', 'Markersize', 10, 'LineWidth', 3, 'DisplayName', 'FA-UPA')
    plot(P_t_arr_dBm, SE_only_active_averaged, '-pent', 'Markersize', 10, 'LineWidth', 3, 'DisplayName', 'FA-ULA')
    legend
    xlabel('Transmit power P_{max} (dBm)')
    ylabel('Spectral efficiency (bps/Hz)')
    grid on
    title('Spectral efficiency vs transmit power for different architectures')
    
    
    %% Plot energy efficiency vs transmit power
    
    figure
    plot(P_t_arr_dBm, EE_MAMP_averaged, '-x', 'Markersize', 10, 'LineWidth', 3, 'DisplayName', 'HRP-UPA (Proposed)')
    hold on
    plot(P_t_arr_dBm, EE_hybrid_PS_averaged, '-v', 'Markersize', 10, 'LineWidth', 3, 'DisplayName', 'HPS-UPA')
    plot(P_t_arr_dBm, EE_fully_active_averaged, '-o', 'Markersize', 10, 'LineWidth', 3, 'DisplayName', 'FA-UPA')
    plot(P_t_arr_dBm, EE_only_active_averaged, '-pent', 'Markersize', 10, 'LineWidth', 3, 'DisplayName', 'FA-ULA')
    legend
    xlabel('Transmit power P_{max} (dBm)')
    ylabel('Energy efficiency (bps/Hz/W)')
    grid on
    title('Energy efficiency vs transmit power for different architectures')

    
end


if vary_number_of_parasitics==1 && vary_transmit_power==1
    figure
    plot([0, squeeze(Np_arr)],[SE_only_active_averaged(1,1,:,1); squeeze(SE_MAMP_averaged(:,1,:,:))] , '-ro','Markersize', 9,'LineWidth', 3, 'DisplayName', "P_{max}=" +num2str(P_t_arr_dBm(1))+" dBm")
    hold on
    plot([0, squeeze(Np_arr)],[SE_only_active_averaged(1,2,:,1); squeeze(SE_MAMP_averaged(:,2,:,:))] , '-bpent','Markersize', 9,'LineWidth', 3, 'DisplayName', "P_{max}=" +num2str(P_t_arr_dBm(2))+" dBm")
    plot([0, squeeze(Np_arr)],[SE_only_active_averaged(1,3,:,1); squeeze(SE_MAMP_averaged(:,3,:,:))] , '-gsq','Markersize', 9,'LineWidth', 3, 'DisplayName', "P_{max}=" +num2str(P_t_arr_dBm(3))+" dBm")
    legend
    xlabel('Number of parasitic antennas')
    ylabel('Spectral efficiency (bps/Hz)')
    grid on
    title('Spectral efficiency vs number of parasitic elements for different transmit powers')
end

if vary_number_of_active_antennas ==1
    figure
    plot(Na_arr, squeeze(SE_MAMP_averaged), '-x', 'Markersize', 10, 'LineWidth', 3, 'DisplayName', 'HRP-UPA (Proposed)')
    hold on
    plot(Na_arr, squeeze(SE_hybrid_ps_arr_averaged), '-v', 'Markersize', 10, 'LineWidth', 3, 'DisplayName', 'HPS-UPA')
    plot(Na_arr, squeeze(SE_fully_active_averaged), '-o', 'Markersize', 10, 'LineWidth', 3, 'DisplayName', 'FA-UPA')
    plot(Na_arr, squeeze(SE_only_active_averaged), '-pent', 'Markersize', 10, 'LineWidth', 3, 'DisplayName', 'FA-ULA')
    legend
    xlabel('Number of active antennas')
    ylabel('Spectral efficiency (bps/Hz)')
    grid on
    title('Spectral efficiency vs number of active antennas for different architectures')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Functions 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Z= self_impedance_dip(k0, l, a)
%% Self impedance of a dipole antenna (eqn 8.57 from Balanis Antenna theory textbook 4th edition)
C=0.5772; %% Euler constant 
epsilon=8.854*10^(-12); %% permittivity of free space
c=3*10^8; %% speed of light
R= (  C+ log(k0*l) - cosint(k0*l)   +0.5*sin(k0*l)*(sinint(2*k0*l)  - 2*sinint(k0*l)) +0.5*cos(k0*l)*(C+log(k0*l/2) ...
     + cosint(2*k0*l)  -2*cosint(k0*l))   )/(2*pi*epsilon*c);

X= ( 2*sinint(k0*l) + cos(k0*l)*(2*sinint(k0*l)- sinint(2*k0*l))- ...
    sin(k0*l)*(2*cosint(k0*l)  - cosint(2*k0*l)-   cosint(2*k0*a^2/l))  )/(4*pi*epsilon*c);

Z=R+1j*X;
end

function Z= MI_sbs_dip(d2lambda, l2lambda)
%% Mutual  impedance with the side by side configuration (eqn 8.68 from Balanis Antenna theory textbook 4th edition)
mu_0=4*pi*10^(-7);
epsilon_0=8.854*10^(-12); 
eta=sqrt(mu_0/ epsilon_0);

u0=2*pi*d2lambda;
u1=2*pi*(sqrt(d2lambda.^2+ l2lambda^2) + l2lambda);
u2=2*pi*(sqrt(d2lambda.^2+ l2lambda^2) - l2lambda);

R= eta./(4*pi).*(2*cosint(u0)-cosint(u1)-cosint(u2));
X= -eta./(4*pi).*(2*sinint(u0)- sinint(u1)-sinint(u2));
Z=R+1j*X;
end


function [Z_p, Z_0, Z_00, parasitic_antenna_indices ]=compute_Zp_Zo_from_Zmat(Z_mat, N, index_active_antenna )

parasitic_antenna_indices=setdiff([1:1:N], index_active_antenna);
Z_p= Z_mat(parasitic_antenna_indices, parasitic_antenna_indices);
Z_0= Z_mat(parasitic_antenna_indices, index_active_antenna);
Z_00=  Z_mat(index_active_antenna, index_active_antenna);
end

function R_loss=  compute_loss_resistance(k, l, rho, f, sigma_c)
%% Code for computing the resistive loss in the dipole array (eqn (7) in the paper "Superdirective Antenna Pairs for Energy-Efficient Terahertz Massive MIMO"
mu=4*pi*10^(-7);
R_loss= (k*l-sin(k*l))/(4*k*rho*(sin(k*l/2))^2)*sqrt(f*mu/(pi*sigma_c));

end


function [h, hp, ha]= generate_multipath_channel(L, channel_gain_var, theta_min, theta_max, N, parasitic_antenna_indices, index_active_antenna, f_c, dx, dy,Na, Np,  c)

theta_arr=unifrnd(theta_min, theta_max, [L,1]);
%theta_arr=1.977467242167018;
gamma_l= sqrt(channel_gain_var)/2*(randn(L,1)+1j*randn(L,1));
%gamma_l=0.916942507297543 - 1.129423430501824i;
h=zeros(N,1);
for l=1:L
    ay_theta=exp(-1j*2*pi*f_c*dy/c*cos(theta_arr(l)).*[0:1:Na-1]).';
    ax_theta=exp(1j*2*pi*f_c*dx/c*sin(theta_arr(l)).*(parasitic_antenna_indices(1:1:Np)-index_active_antenna(1))).';
    axy_theta=[ay_theta.', (kron(ay_theta, ax_theta)).'].';
    h=h+1/sqrt(L)*gamma_l(l)*axy_theta;
end
ha=h(1:Na);
hp=h(Na+1:1:end);

end


function F_PS_subconnected=  construct_subconnected_F_PS_matrix(f_ps, NA, NP)
%% Code to generate the beamforming matrix for the sub-connected phase shifter architecture
F_PS_subconnected=zeros(NA,  NA*(NP+1));
for ni=1:NA
    for np=1:NP+1
        if np==1
            F_PS_subconnected(ni, ni)= f_ps(ni);
        else 
            F_PS_subconnected(ni, NA+ NP*(ni-1)+np-1)= f_ps(NA+ NP*(ni-1)+np-1);
            
        end
    end
end

end
