function POST19_PDF_CDF_Multi_Published

clc;
clear all;
[inputFileName,PathName] = uigetfile('*.mat',...
    'Select the file in order of low Nsamples to high Nsamples'...
    ,'MultiSelect', 'on');
fprintf(['Input Files:',num2str(length(inputFileName))]);fprintf('\n');

% leg1={};leg2={};
MyInputLegend={};
InputPDF=[];
InputCDF=[];
Input_t=[];
MyOutputLegend={};
OutputPDF=[];
OutputCDF=[];
Output_t=[];
for i1=1:length(inputFileName)
    load([PathName,inputFileName{i1}]);
    x=ud.x;F=ud.F;Fs=ud.Fs;Frms=ud.Frms;fmax=ud.fmax;InputType = ud.InputType;
    disp(inputFileName{i1})
    fprintf(['Fs:',num2str(Fs)]);
    fprintf(['  ,RMS:',num2str(Frms)]);
    fprintf(['  ,Samples:',num2str(length(F))]);
    fprintf('\n');
    x=x(:,1);
    Xrms = norm(x)/sqrt(length(x));
    Nsample=length(F);
    
    N1 = Nsample/256;% choose 16,32,64,128 ,... 
    LegAdder1 = [];
    LegAdder1 = [LegAdder1,'Frms:',num2str(Frms),','];
%     LegAdder1 = [LegAdder1,['Nsample:',num2str(Nsample),',Fs:',num2str(Fs),',']];
%     LegAdder1 = [LegAdder1,['fmax=',num2str(fmax),',']];

    [t PDF CDF] = PDF_CDF(F,N1);
    InputPDF=[InputPDF PDF'];
    InputCDF=[InputCDF CDF'];
    Input_t=[Input_t t'];
    MyInputLegend = [MyInputLegend,LegAdder1];

    [t PDF] = Noramal_Distribution(F,N1);
    InputPDF=[InputPDF PDF'];
    InputCDF=[InputCDF NaN(length(t),1)];
    Input_t=[Input_t t'];
    MyInputLegend = [MyInputLegend,'Normal Distribution'];
    
    LegAdder2 = [];
    LegAdder2 = [LegAdder2,',Xrms:',num2str(Xrms)];
    LegAdder2 = [LegAdder2,[',Frms:',num2str(Frms)]];
%     LegAdder2 = [LegAdder2,[',Nsample:',num2str(Nsample),',Fs:',num2str(Fs)]];
%     LegAdder2 = [LegAdder2,[',fmax=',num2str(fmax)]];
     
    [t PDF CDF] = PDF_CDF(x,N1);
    OutputPDF=[OutputPDF PDF'];
    OutputCDF=[OutputCDF CDF'];
    Output_t=[Output_t t'];
    MyOutputLegend = [MyOutputLegend,LegAdder2];

    [t PDF] = Noramal_Distribution(x,N1);
    OutputPDF=[OutputPDF PDF'];
    OutputCDF=[OutputCDF NaN(length(t),1)];
    Output_t=[Output_t t'];
    MyOutputLegend = [MyOutputLegend,'Normal Distribution'];

    XrmsVec(i1)=Xrms;
    FrmsVec(i1)=Frms;
end

h1=figure('Name',[InputType,' PDF Plots For Input Signals']);
h2=figure('Name',[InputType,' CDF Plots For Input Signals']);
h3=figure('Name',[InputType,' PDF Plots For output Signals']);
h4=figure('Name',[InputType,' CDF Plots For output Signals']);

figure(h1);plot(Input_t,InputPDF);
legend(MyInputLegend)
figure(h2);plot(Input_t,InputCDF);
legend(MyInputLegend)

figure(h3);plot(Output_t,OutputPDF);
legend(MyOutputLegend)
figure(h4);plot(Output_t,OutputCDF);
legend(MyOutputLegend)

figure('Name','Output RMS Vs Input RMS');
plot(FrmsVec,XrmsVec,'-o');
xlabel('Input RMS');ylabel('Output RMS');

figure('Name','Input&Outpt PDF')
plot(Input_t(:,1:2:end),InputPDF(:,1:2:end));hold on;grid on;
% Mylegend = MyInputLegend{1:2:end};
% legend(MyInputLegend{1:2:end})
plot(Output_t(:,1:2:end),OutputPDF(:,1:2:end));
% Mylegend=[MyInputLegend{1:2:end},MyOutputLegend{1:2:end}];
legend(MyInputLegend{1:2:end},MyOutputLegend{1:2:end})
end


function [t PDF CDF] = PDF_CDF(signal,N1)
Nsample = length(signal);
%%PDF
N=floor(Nsample/N1);
A=min(signal); B=max(signal);
Delta=(B-A)/N;
t=A-Delta/2+[1:N]*Delta; 
PDF=hist(signal,t)/(Delta*Nsample);
%%CDF
p=hist(signal,t)/Nsample;
CDF=cumsum(p);
end
function [t Guassian_dist] = Noramal_Distribution(signal,N1)

Nsample = length(signal);
N=floor(Nsample/N1);
A=min(signal); B=max(signal);
Delta=(B-A)/N;
t=A-Delta/2+[1:N]*Delta; 
%Guassian pdf
M = mean(signal);S = std(signal);
Guassian_dist = 1/sqrt(2*pi)/S * exp(-(t-M).^2/2/S^2);
% Guassian_dist = pdf('Normal',x,M,S);
% figure(h1);plot(t,Guassian_dist,'r--');legend('PDF','Guassian Distribution');
% figure(h2);plot(t1,Guassian_dist,'r--');legend('PDF','Guassian Distribution')
end