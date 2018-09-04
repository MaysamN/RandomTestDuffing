function POST6_WelchEstimates_And_RawEstimate_Published
% In the Name of ALLAH
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

fmin=0;fmax=6;
myxlim = [fmin fmax];
clear ud;
F = reshape(F,length(F),1);
x = reshape(x(:,1),length(x),1);

%% printing inputs
fprintf('\n');
fprintf(['Nsample:',num2str(length(F)),',WinLen:',num2str(length(win))...
        ,',Overlap:',num2str(overlap)]);fprintf('\n');
%% Welch estimates - cpsd method 
[PFx f1] = cpsd(F,x,win,overlap,length(win),Fs);
[PxF f1] = cpsd(x,F,win,overlap,length(win),Fs);
[Pxx f1] = cpsd(x,x,win,overlap,length(win),Fs);
[PFF f1] = cpsd(F,F,win,overlap,length(win),Fs);
% Hcpsd1 = (PFx)./(PFF);
Hcpsd1 = (PxF)./(PFF);
Hcpsd2 = (Pxx)./(PFx);
Hcpsd3 = sqrt((Pxx)./(PFF));
Hcpsdv = (Hcpsd1+Hcpsd2)/2;

%% Raw estimate
win = rectwin(length(F));overlap =0;
[PxF fraw] = cpsd(x,F,win,overlap,length(win),Fs);
[PFF fRaw] = cpsd(F,F,win,overlap,length(win),Fs);
HRaw = (PxF)./(PFF);

%% Linear FRF
df = Fs/length(F);f2 = (0:df:Fs-df);
Hlinear = 1./(-m*((2*pi*f2).^2)+c*1i*2*pi*f2+k); 

Lg1 = [',',num2str(Frms),',',num2str(length(F)),...
            ',(',num2str(length(win)),')'];
       
Lg2 = [',',num2str(Frms),',',num2str(length(F))];

%% FRF magnitude
figure('Name','FRF Magnitude');
plot(fRaw,db(HRaw),'.g',f1,db(Hcpsd1),f1,db(Hcpsd2),f1,db(Hcpsd3),f1,db(Hcpsdv),f2,db(Hlinear));
legend(['Raw',Lg2],['H1',Lg1],['H2',Lg1],['H3',Lg1],['Hv',Lg1],'Linear');
xlabel('Frequency');ylabel('FRF(db)');
xlim([fmin+df fmax-df]);

%% phase : H1 and H2 and Hv angle are the same for welch and raw estimates
figure('Name','phase 1');
plot(fRaw,angle(HRaw),'.g',f1,angle(Hcpsd1),'.',f1,angle(Hcpsd2),'o',f1,angle(Hcpsdv),'d',f2,angle(Hlinear),'--');
legend(['Raw',Lg1],['H1',Lg1],['H2',Lg1],['Hv',Lg1],'Linear');
xlabel('Frequency(Hz)');ylabel('Phase(rad)');
title('Phase Plot');
xlim(myxlim);

%% Nyquist diagram
fmin=0;fmax=6;
% indexmax1 = length(find(f1<fmax));
indexmax1 = sum(f1<fmax);
indexmin1 = sum(f1<fmin)+1;
% f1 = f1(indexmin1:indexmax1);
Hcpsd1 = Hcpsd1(indexmin1:indexmax1);
Hcpsd2 = Hcpsd2(indexmin1:indexmax1);
Hcpsdv = Hcpsdv(indexmin1:indexmax1);

indexmax3 = sum(fRaw<fmax)-1;
indexmin3 = sum(fRaw<fmin)+2;
% fRaw = fRaw(indexmin3:indexmax3);
HRaw = HRaw(indexmin3:indexmax3);

figure('Name','Nyquist 1');
plot(real(HRaw),imag(HRaw),'.g',real(Hcpsd1),imag(Hcpsd1),'.-',real(Hcpsd2),imag(Hcpsd2),'o-',...
    real(Hcpsdv),imag(Hcpsdv),'-d',real(Hlinear),imag(Hlinear),'--');grid on;
xlabel('real');ylabel('Image');
title('Nyquist diagram');
legend(['Raw',Lg1],['H1',Lg1],['H2',Lg1],['Hv',Lg1],'Linear');

end