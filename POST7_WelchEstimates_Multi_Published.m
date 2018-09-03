function POST7_WelchEstimates_Multi_Published
% In the Name of ALLAH

clc;
clear all;
fprintf('In the Name of ALLAH');fprintf('\n');
%% inputs:
% winLenMat = 2.^[10 10 10 10];
winLenMat = 2.^[12 12 12 12];

% winLenMat = 2.^[14 13 11 10];
fmax = 6;
fmin = 0; 

%% Loading files
[inputFileName,PathName] = uigetfile('*.mat',...
    'Select the file in order of low Nsamples to high Nsamples'...
    ,'MultiSelect', 'on');
fprintf(['Input Files:',num2str(length(inputFileName))]);fprintf('\n');

h1 = figure('Name','FRF Magnitude for h1');
h2 = figure('Name','FRF Magnitude for h2');
ph = figure('Name','pahse plot');
Ny1 = figure('Name','Nyquist for h1');
Ny2 = figure('Name','Nyquist for h2');
hb_Ny = figure('Name','Nyquist for hb');
lengendMat_h1={};lengendMat_h2={};lengendMat_hb={};

for ii = 1:length(inputFileName)
    load([PathName,inputFileName{ii}]);
    disp(ud);
    Frms = ud.Frms;
    m=ud.m;k=ud.k;c=ud.c;knl=ud.knl;
    x=ud.x;F=ud.F;Fs=ud.Fs;
%     winlength = ud.BestBiased;
%     winlength = ud.BestSmoothed;
    winlength = winLenMat(ii);
    myudfields = fieldnames(ud);
    clear ud;
  
    % Reshaping
    F = reshape(F,length(F),1);
    x = reshape(x(:,1),length(x),1);
    win = hanning(winlength);overlap =length(win)/2;

    % printing inputs
    fprintf('\n');
    fprintf(['Nsample:',num2str(length(F)),',WinLen:',num2str(length(win))...
        ,',Overlap:',num2str(overlap)]);
    fprintf('\n');
    %% cpsd method
    [PFx f1] = cpsd(F,x,win,overlap,length(win),Fs);
    [PxF f1] = cpsd(x,F,win,overlap,length(win),Fs);
    [Pxx f1] = cpsd(x,x,win,overlap,length(win),Fs);
    [PFF f1] = cpsd(F,F,win,overlap,length(win),Fs);
    Hcpsd1 = (PxF)./(PFF);
    Hcpsd2 = (Pxx)./(PFx);
    Hcpsd3 = sqrt((Pxx)./(PFF));
    Hcpsdv = (Hcpsd1+Hcpsd2)/2;
    
    indexmax1 = sum(f1<fmax);
    
    indexmax1 = indexmax1-0;
    indexmin = sum(find(f1<fmin))+1;
    
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
    f2 = f2(1:indexmax2);
    Hlinear = Hlinear(1:indexmax2);
    H_hb = H_hb(1:indexmax2);
      
    lengendMat_h1 = [lengendMat_h1,['H1',',',num2str(Frms),',',num2str(length(F)),',(',num2str(length(win)),')']];
    lengendMat_h2 = [lengendMat_h2,['H2',',',num2str(Frms),',',num2str(length(F)),',(',num2str(length(win)),')']];
    lengendMat_hb = [lengendMat_hb,['hb',',',num2str(Frms),',',num2str(length(F)),',(',num2str(length(win)),')']];
             
    figure(h1);
    plot(f1,db(Hcpsd1));hold on; grid on;
    xlabel('Frequency');ylabel('FRF(db)');
    
    figure(h2);
    plot(f1,db(Hcpsd2));hold on; grid on;
    xlabel('Frequency');ylabel('FRF(db)');
    
    figure(ph);
    plot(f1,angle(Hcpsd1),'.'); hold on; grid on;
    xlabel('Frequency(Hz)');ylabel('Phase(rad)');
    title('Phase Plot');
        
    figure(Ny1);
    plot(real(Hcpsd1),imag(Hcpsd1)); hold on; grid on;
    xlabel('real');ylabel('Image');
    title('Nyquist diagram');
    
    figure(Ny2);
    plot(real(Hcpsd2),imag(Hcpsd2)); hold on; grid on;
    xlabel('real');ylabel('Image');
    title('Nyquist diagram');
    
    figure(hb_Ny);
    plot(real(H_hb),imag(H_hb)); hold on; grid on;
    xlabel('real');ylabel('Image');
    title('Nyquist diagram');

end

figure(h1);
legend(lengendMat_h1)
figure(h2);
legend(lengendMat_h2)
figure(Ny1);
legend(lengendMat_h1)
figure(Ny2);
legend(lengendMat_h2)
figure(ph);
legend(lengendMat_h1)
figure(hb_Ny);
legend(lengendMat_hb)

end


