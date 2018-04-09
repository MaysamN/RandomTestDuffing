clc;
clearvars;

zetaVec = [0.05 0.1 0.2 0.3 0.4 0.5 0.7 1 1.5];
r = (0:0.01:5);

h1 = figure;h2 = figure;h3 = figure;h4 = figure;h5 = figure;
for i=1:length(zetaVec)
    zeta = zetaVec(i);
    
    X_Force = 1./sqrt((1-r.^2).^2+(2*zeta*r).^2);
    X_Unblc = r.^2.*X_Force ;                  % Rotating Unbalance
    X_Base  = sqrt(1+(2*zeta*r).^2).*X_Force ; % Base Excitation
                                               % Support Reaction
    phi = atand(2*zeta*r./((1-r.^2)));
    phi = phi+(phi<0)*180; % to smooth the phase plot
    
    psi = atand(2*zeta*r.^3./((1-r.^2)+(2*zeta*r).^2));
    psi = psi+(psi<0)*180; 
    
    figure(h1);plot(r,X_Force,'k');hold on;grid on;
    figure(h2);plot(r,X_Unblc,'k');hold on;grid on;
    figure(h3);plot(r,X_Base,'k');hold on;grid on;
    figure(h4);plot(r,phi,'k');hold on;grid on;
    figure(h5);plot(r,psi,'k');hold on;grid on;
end
figure(h1);xlabel('\omega/\omega_n');ylabel('Xk/F0');ylim([0 3.5]);
figure(h2);xlabel('\omega/\omega_n');ylabel('MX/me or Z/Y');ylim([0 3.5]);
figure(h3);xlabel('\omega/\omega_n');ylabel('X/Y');ylim([0 3.5]);
figure(h4);xlabel('\omega/\omega_n');ylabel('\phi');
figure(h5);xlabel('\omega/\omega_n');ylabel('\psi');