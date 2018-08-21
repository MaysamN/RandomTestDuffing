function SDOF6_cubic_multisine_Version1
clc;
clear all;
fprintf('In The Name of ALLAH');fprintf('\n');
fprintf('SDOF_cubic_MultiSine');fprintf('\n');
MyMemory = memory;disp(MyMemory);
Myclock=clock;
Myclock=[num2str(Myclock(4)),':',num2str(Myclock(5))];
fprintf(['Started @ ',Myclock]);fprintf('\n');
%% inputs
k = 10;knl = 1;m = 1;c = 0.01;
Frms = sqrt(2)/2; Nsample =2^13;Fs = 30;%Fs:sampling frequency
dt = 1/Fs;
t = (0:dt:(Nsample-1)*dt);
% t vector is used in ode45 solver.
df = Fs/Nsample;f = (df:df:Fs);

Deltaf = Fs/Nsample;
fmax = 6; % this parameter is used in excitation signal production
fmin = 0;%Fs/Lines;% this parameter is used in excitation signal production
% should be set multiples of 1/Nsample and should not be set to zero.
fprintf(['m=',num2str(m)]);fprintf(['    k=',num2str(k)]);
fprintf(['    c=',num2str(c)]);fprintf(['    knl=',num2str(knl)]);
fprintf('\n');
fprintf(['RMS:',num2str(Frms)]);
fprintf(['     Sampling Freq:',num2str(Fs)]);
fprintf(['     Nsample:',num2str(Nsample)]);
fprintf(['    fmin:',num2str(fmin)]);
fprintf(['    Deltaf:',num2str(Deltaf)]);
fprintf(['    fmax:',num2str(fmax)]);fprintf('\n');
%% Excitation signal
F = MultisineSignal_Normalized_ifft(fmin/Fs,fmax/Fs,Deltaf/Fs,Nsample);F=F';

%% setting RMS 
RMS_F = norm(F)/sqrt(length(F));
F=Frms*(F/RMS_F); 
%% checking the created input signal (not necessary)
SFF = (fft(F).*conj(fft(F)))/Nsample/Fs;% double-sided power spectral density
disp(sqrt(sum(SFF)*df));%it should be equal to RMS
figure;
subplot(2,1,1);plot(t,F);xlabel('t(sec)');ylabel('F');xlim([0 t(length(t))]);
% subplot(2,1,2);plot(f,SFF);xlabel('Frequency(Hz)');ylabel('double-sided PSD([F]^2/Hz)');
subplot(2,1,2);plot(f(1:Nsample/2),2*SFF(1:Nsample/2));xlabel('Frequency(Hz)');ylabel('single-sided PSD([F]^2/Hz)');

%% solving by rung-kutta method
fprintf('End Time:');fprintf('%f \n',t(length(t)));
tstart=tic;
initial_condition = [0;0];
[T x] = nonlinear_dynamic_solution(k,knl,m,c,t,F,initial_condition);
MyRunTime = toc(tstart);

%% saving the result
Myclock=clock;
Myclock=[num2str(Myclock(4),'%02.0f'),num2str(Myclock(5),'%02.0f')];
filename = ['SDOF_cubic_MultiSine(',Myclock,'-',date,')'...
            ,'_m=',num2str(m),'_k=',num2str(k),...
             '_c=',num2str(c),'_knl=',num2str(knl)...
            ,'_Fs=',num2str(Fs),'(fmax=',num2str(fmax),...
            ')_RMS=',num2str(Frms),'_Nsample=',num2str(Nsample),'.mat'];
ud.InputType = 'MultiSine';
ud.k=k;ud.m=m;ud.c=c;ud.knl=knl;
ud.Fs = Fs;ud.Frms = Frms;
ud.fmin=fmin;ud.fmax=fmax;ud.Deltaf=Deltaf;
ud.x = x;ud.F = F;ud.t = t;


ud.Mydate=date;ud.Myclock=Myclock;ud.MyRunTime=MyRunTime;

save(filename,'ud');

%%%%%%%%% simple post process
win = hanning(4096);overlap = length(win)/2;
df2 = (Fs/length(win));
f2 =df2*(0:df2:round(length(win))/2);
myxlim = [0 fmax];

%%%% cpsd
% --> It is very important to use length(win) instead of [] in the fifth
% input of cpsd and tfe.
[PFx f1] = cpsd(F,x(:,1),win,overlap,[],Fs);
[PxF f1] = cpsd(x(:,1),F,win,overlap,[],Fs);
[Pxx f1] = cpsd(x(:,1),x(:,1),win,overlap,[],Fs);
[PFF f1] = cpsd(F,F,win,overlap,[],Fs);
H_cpsd1 = (PxF)./(PFF);
H_cpsd2 = (Pxx)./(PFx);
H_cpsd3 = sqrt((Pxx)./(PFF));
%handle3 = figure('name','cpsd: averaged');
figure('Name','cpsd');
plot(f1,db(H_cpsd1),f1,db(H_cpsd2),f1,db(H_cpsd3));grid on;hold on;
xlabel('Frequency(Hz)');ylabel('FRF(dB)');xlim(myxlim);
legend('H1 (cpsd)','H2 (cpsd)','H3 (cpsd)');

%%%% tfestimate 
[H_tfe1 f1] = tfestimate(F,x(:,1),win,overlap,[],Fs);
[H_tfe2 f1] = tfestimate(x(:,1),F,win,overlap,[],Fs);H_tfe2 = H_tfe2.^-1;
figure('Name','tfestimate');
plot(f1,db(H_tfe1),f1,db(H_tfe2)); hold on; grid on;
xlabel('Frequency(Hz)');ylabel('FRF(dB)');xlim(myxlim);

% linear FRF
Hlinear = 1./(-m*((2*pi*f1).^2)+c*1i*2*pi*f1+k); 
plot(f1,db(Hlinear),'r');hold on;grid on;
legend('H1 (tfestimate)','H2 (tfestimate)','Linear System');

%% PDF and CDF plots
PDF_CDF_plot(x(:,1),100);
PDF_CDF_plot(F,100);

end

function F = MultisineSignal_Normalized_ifft(f0,f1,df,Nsample)

tic;
fsample=1;Ts=1/fsample;         % sample frequency and sample period
% fsample is 1 in Normalized method
%Nsample: length of the signal
% time=[0:Nsample-1]*Ts;              % time axis
% Nsines=round(f1*Nsample)-round(f0*Nsample);          % number of sines
Nsines=length(1+round(f0*Nsample):round(df*Nsample):round(f1*Nsample));
% f=(0:Nsample-1)*fsample/Nsample;   % multisine frequencies
% f=(0:Nsample-1)*df;   % multisine frequencies

U2=zeros(Nsample,1);               % choose random phases
F=zeros(Nsample,1);               % choose random phases
% U2(1+round(f0*Nsample):round(f1*Nsample))=exp(1i*2*pi*rand(Nsines,1));
U2(1+round(f0*Nsample):round(df*Nsample):round(f1*Nsample))...
    =exp(1i*2*pi*rand(Nsines,1));
F=2*real(ifft(U2));
F=F/std(F);
% U2m=fft(F)/sqrt(Nsample);       % spectrum of the actual generate multisine
% wF= 2*pi*(f0:df:f1);
% phiF = (2*rand(1,length(wF))-1)*pi;
Fruntime = toc;
fprintf('F is created in %10.1f secs',Fruntime);fprintf('\n');
%%%%%%%%%%%%%%%%%%%%%%%%print and write%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filename = ['MultisineNormalized_Excitation_ifft','_Nsample=',num2str(Nsample)];
% save(filename,'F');

end
function [t x] = nonlinear_dynamic_solution(k,knl,m,c,tvector,F,initial_condition)
dt = tvector(2)-tvector(1);
N = length(tvector);

Refine = 1;
RelTol = 1e-5;
MaxStep = 0.01;
InitialStep = 0.01;

graphic = 'odeprint';... 'odeprint' 'odeplot' 'none'
F_calculation_method ='interpolation';...'interpolation' 'exact'
odetype = 'ode45';... 'ode45' 'ode23' 'ode23s'

switch graphic
    case 'odeplot'
        Outputsel = [1 2];
        options = odeset('Mass',[1 0;0 m],'Jacobian',[0 1 ;-k -c],'MaxStep',MaxStep,'Refine',Refine,...
            'InitialStep',InitialStep,...
            'OutputFcn',@odeplot,'OutputSel',Outputsel,'RelTol',RelTol,'stats','on');        
        figure;
    case 'odeprint'
%         Outputsel = [1 2];
        options = odeset('Mass',[1 0;0 m],'Jacobian',[0 1 ;-k -c],'MaxStep',MaxStep,'Refine',Refine,...
            'InitialStep',InitialStep,...
            'OutputFcn',@myodeprint,'RelTol',RelTol,'stats','on');...'OutputSel',Outputsel,
    case 'none'
        fprintf('RK starts');fprintf('\n');
        options = odeset('Mass',[1 0;0 m],'Jacobian',[0 1 ;-k -c],'MaxStep',MaxStep,'Refine',Refine,...
            'InitialStep',InitialStep,'RelTol',RelTol,'stats','on');
end
switch F_calculation_method
    case 'interpolation'
        switch odetype
            case 'ode45'                 
                [t,x] = ode45(@func_Rungfunc_interpolation,tvector,initial_condition ,options );                
            case 'ode23'
                [t,x] = ode23(@func_Rungfunc_interpolation,tvector,initial_condition ,options );
            case 'ode23s'
                [t,x] = ode23s(@func_Rungfunc_interpolation,tvector,initial_condition ,options );
        end
    case 'exact'
         odetype = 'ode45';
        [t,x] = ode45(@func_Rungfunc_exact,tvector,initial_condition ,options );
end

    function dy = func_Rungfunc_exact(t,y)
        F1 = FG*sin(omega*t);
        % adding nonlinear part as an extra excitation
        F1(2*node_spring-1) = F1(2*node_spring-1) - k_cubic * y(2*node_spring-1)^3;
        dy = [O I ;-KG -CG] * y + [OF;F1];        
    end
    function dy = func_Rungfunc_interpolation (t,y)
        tstep_counter = floor((t-tvector(1))/dt)+1;
        if tstep_counter==N
            %%extrapolation
            tstep_counter =tstep_counter-1;
        end
        t1 = tvector(tstep_counter);t2 = tvector(tstep_counter+1) ;
        F1 = F(tstep_counter);F2 = F(tstep_counter+1) ;
        Ft = F1 + (F2-F1)*(t-t1)/(t2-t1);
%         Ft = Ft * FG;
        Ft = Ft - knl * y(1)^3;% adding nonlinear part as an excitation force
        dy = [0 1 ;-k -c] * y + [0;Ft];
    end
    function status = myodeprint(t,y,flag,varargin)
        %ODEPRINT
        if nargin < 3 || isempty(flag) % odeprint(t,y) [v5 syntax] or odeprint(t,y,'')
            fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%015.5f',t);
        else
            switch(flag)
                case 'init'               % odeprint(tspan,y0,'init')
                    fprintf('t:')
                    fprintf('%015.5f',t(1))
                    
                case 'done'               % odeprint([],[],'done')
                    
                    fprintf('\n\n');
                    
            end
        end
        
        status = 0;
    end

end
function PDF_CDF_plot(signal,N1)
Nsample = length(signal);
%%PDF
figure;subplot(2,1,1);
N=floor(Nsample/N1);
A=min(signal); B=max(signal);
Delta=(B-A)/N;
t=A-Delta/2+[1:N]*Delta; f=hist(signal,t)/(Delta*Nsample);
bar(t,f);
title('Estimated PDF');
%%CDF
subplot(2,1,2);
p=hist(signal,t)/Nsample;
CDF=cumsum(p);
plot(t,CDF)
title('Estimated CDF')
end
