function POST6_WelchEstimates_H1H2H3Hv
% In the Name of ALLAH

% In this function we want to show the differences between the 4 famous FRF
% estimators. the user can use this function for linear and nonlinear
% multisine excited or sweptsine excited systems. 
% one should note that for H3 estimator we lose the phase informations.

clc;
clear all;
fprintf('In the Name of ALLAH');fprintf('\n');
%% inputs:
win = hanning(2*1024);overlap =length(win)/2;
%% Raw estimate
% win = rectwin(length(F));overlap =0;
%% Loading a data file 
uiopen('LOAD');
disp(ud);
Frms = ud.Frms;
m=ud.m;k=ud.k;c=ud.c;knl=ud.knl;
x=ud.x;F=ud.F;Fs=ud.Fs;
myudfields = fieldnames(ud);
% Checking the data file for 'fmax'
if sum(strcmp(myudfields,'fmax'))==1
    fmax=ud.fmax;
    fmin=ud.fmin;
    fprintf('fmax is read from file');fprintf('\n');
else
    fmax=6;% in case when the ud excludes fmax
    fmin=0;
    fprintf('fmax is set by user');fprintf('\n');
end
fmin=0;fmax=6;
myxlim = [fmin fmax];
clear ud;
% Reshaping
F = reshape(F,length(F),1);
x = reshape(x(:,1),length(x),1);

%% printing inputs
fprintf('\n');
fprintf(['Nsample:',num2str(length(F)),',WinLen:',num2str(length(win))...
        ,',Overlap:',num2str(overlap)]);
fprintf('\n');
%% cpsd method 
[PFx f1] = cpsd(F,x,win,overlap,length(win),Fs);
[PxF f1] = cpsd(x,F,win,overlap,length(win),Fs);
[Pxx f1] = cpsd(x,x,win,overlap,length(win),Fs);
[PFF f1] = cpsd(F,F,win,overlap,length(win),Fs);
% Hcpsd1 = (PFx)./(PFF);
Hcpsd1 = (PxF)./(PFF);
Hcpsd2 = (Pxx)./(PFx);
Hcpsd3 = sqrt((Pxx)./(PFF));
Hcpsdv = (Hcpsd1+Hcpsd2)/2;

% indexmax1 = length(find(f1<fmax));
indexmax1 = sum(f1<fmax);
indexmax1 = indexmax1-0;
indexmin = sum(f1<fmin)+1;
% indexmin =1;
f1 = f1(indexmin:indexmax1);
Hcpsd1 = Hcpsd1(indexmin:indexmax1);
Hcpsd2 = Hcpsd2(indexmin:indexmax1);
Hcpsd3 = Hcpsd3(indexmin:indexmax1);
Hcpsdv = Hcpsdv(indexmin:indexmax1);

df = Fs/length(F);f2 = (0:df:Fs-df);
Hlinear = 1./(-m*((2*pi*f2).^2)+c*1i*2*pi*f2+k); 

RMS_x = norm(x(:,1))/sqrt(length(x(:,1)));
H_hb = 1./(-m*((2*pi*f2).^2)+c*1i*2*pi*f2+k+(0.75*knl* (1.4286*RMS_x)^2));

indexmax2 = length(find(f2<fmax));
indexmin2 = sum(f2<fmin)+1;

f2 = f2(indexmin2:indexmax2);
Hlinear = Hlinear(indexmin2:indexmax2);
H_hb = H_hb(indexmin2:indexmax2);

Lg1 = [',',num2str(Frms),',',num2str(length(F)),...
            ',(',num2str(length(win)),')'];
       
Lg2 = [',',num2str(Frms),',',num2str(length(F))];
%% FRF magnitude
figure('Name','FRF Magnitude 1');
plot(f1,db(Hcpsd1),'.',f1,db(Hcpsd2),'o',f1,db(Hcpsd3),'s',f1,db(Hcpsdv),'d',f2,db(Hlinear),'--',f2,db(H_hb));
legend(['H1',Lg1],['H2',Lg1],['H3',Lg1],['Hv',Lg1],'Linear',['HB',Lg2]);
xlabel('Frequency');ylabel('FRF(db)');
% xlim(myxlim);

% figure('Name','FRF Magnitude 2');
% plot(f1,db(Hcpsd1),f1,db(Hcpsd2),'--',f1,db(Hcpsd3),'-.',f1,db(Hcpsdv),':',f2,db(Hlinear),'--',f2,db(H_hb));
% legend(['H1',Lg1],['H2',Lg1],['H3',Lg1],['Hv',Lg1],'Linear',['HB',Lg2]);
% xlabel('Frequency');ylabel('FRF(db)');
% xlim(myxlim);

%% phase : H1 and H2 and H3 angle are the same for welch and raw estimates
figure('Name','phase 1');
plot(f1,angle(Hcpsd1),'.',f1,angle(Hcpsd2),'o',f1,angle(Hcpsdv),'d',f2,angle(Hlinear),'--',f2,angle(H_hb));
legend(['H1',Lg1],['H2',Lg1],['Hv',Lg1],'Linear',['HB',Lg2]);
xlabel('Frequency(Hz)');ylabel('Phase(rad)');
title('Phase Plot');
% xlim(myxlim);

% figure('Name','phase 2');
% plot(f1,angle(Hcpsd1),f1,angle(Hcpsd2),'--',f1,angle(Hcpsdv),':',f2,angle(Hlinear),'--',f2,angle(H_hb));
% legend(['H1',Lg1],['H2',Lg1],['Hv',Lg1],'Linear',['HB',Lg2]);
% xlabel('Frequency(Hz)');ylabel('Phase(rad)');
% title('Phase Plot');

% xlim(myxlim);
% figure;plot(f1,angle(Hcpsd3));

%% Bode diagram

%% Nyquist diagram

figure('Name','Nyquist 1');
plot(real(Hcpsd1),imag(Hcpsd1),'.',real(Hcpsd2),imag(Hcpsd2),'o',...real(Hcpsd3),imag(Hcpsd3),'s',...
    real(Hcpsdv),imag(Hcpsdv),'d',real(Hlinear),imag(Hlinear),'--',real(H_hb),imag(H_hb));grid on;
xlabel('real');ylabel('Image');
title('Nyquist diagram');
legend(['H1',Lg1],['H2',Lg1],['Hv',Lg1],'Linear',['HB',Lg2]);

figure('Name','Nyquist 2');
plot(real(Hcpsd1),imag(Hcpsd1),real(Hcpsd2),imag(Hcpsd2),'--',...real(Hcpsd3),imag(Hcpsd3),'-.',...
    real(Hcpsdv),imag(Hcpsdv),':',real(Hlinear),imag(Hlinear),'b',real(H_hb),imag(H_hb),'r');grid on;
xlabel('real');ylabel('Image');
title('Nyquist diagram');
legend(['H1',Lg1],['H2',Lg1],['Hv',Lg1],'Linear',['HB',Lg2]);

% figure;plot(real(Hcpsd3),imag(Hcpsd3));

%%
Inputfig = figure('Name','Input Summary','Position',[200 200 380 200]);
dat1 = [m,k,c,knl,fmin,fmax]';
cnames1 = {'m','k','c','knl','fmin','fmax'};
t1 = uitable('Data',dat1,'RowName',cnames1,...'ColName',{'',''},...
    'Parent',Inputfig,'Position',[20 57 132 130]);
dat2 = [Frms,Fs,length(F),length(win),f1(1),f1(end),indexmin,indexmax1]';
cnames2 = {'Frms','Fs','samples','blocksize','fminPlot','fmaxPlot','indexmin','indexmax'};
t2 = uitable('Data',dat2,'RowName',cnames2,...
    'Parent',Inputfig,'Position',[160 20 176 166]);

end