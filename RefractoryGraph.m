clc
close all

%% parameters

landma=150;
m=1;
eta=.002;          %A typical value of the refractory period is 2 ms
delta=.2;
n_max=round(delta/eta)-1;
coefficient =1;

%% ratedgamma distribution

n_max_ratedgamma=coefficient*n_max;
n_ratedgamma=0:1:n_max_ratedgamma;
PDF_ratedgamma=zeros(1,n_max_ratedgamma+1);


for k=0:1:n_max_ratedgamma
    
    PDF_ratedgamma(k+1)=0;
    
    for kk=k*m:1:k*m+m-1
    PDF_ratedgamma(k+1)=PDF_ratedgamma(k+1)+(((landma^kk)*(delta^kk))/(factorial(kk)));    
    end   
    PDF_ratedgamma(k+1)=exp(-1*landma*delta)*PDF_ratedgamma(k+1);
 
end


%% truncated_ratedgamma distribution

n_max_truncated_ratedgamma=n_max;
n_truncated_ratedgamma=0:1:n_max_truncated_ratedgamma;
PDF_truncated_ratedgamma=zeros(1,n_max_truncated_ratedgamma+1); 

normalizationfactor=0;

for k=1:1:n_max_truncated_ratedgamma+1
    normalizationfactor=normalizationfactor+PDF_ratedgamma(k);
 end



 for k=0:1:n_max_truncated_ratedgamma
   PDF_truncated_ratedgamma(k+1)= PDF_ratedgamma(k+1)/normalizationfactor;
 end

 
 
 
 %% refractory_ratedgamma distribution 
 

n_max_refractory_ratedgamma=n_max;
n_refractory_ratedgamma=0:1:n_max_refractory_ratedgamma;

PDF_refractory_ratedgamma1=zeros(1,n_max_refractory_ratedgamma+1); 
PDF_refractory_ratedgamma2=zeros(1,n_max_refractory_ratedgamma+1); 
PDF_refractory_ratedgamma=zeros(1,n_max_refractory_ratedgamma+1); 
 

for k=0:1:n_max_refractory_ratedgamma
  
    
for kk=0:1:(k+1)*m-1
PDF_refractory_ratedgamma1(k+1)=PDF_refractory_ratedgamma1(k+1)+...
(((landma^kk)*((delta-((k+1)*eta))^kk))/(factorial(kk)));
end
PDF_refractory_ratedgamma1(k+1)=exp(-1*landma*(delta-((k+1)*eta)))*PDF_refractory_ratedgamma1(k+1);


for kk=0:1:k*m-1
PDF_refractory_ratedgamma2(k+1)=PDF_refractory_ratedgamma2(k+1)+...
(((landma^kk)*((delta-(k*eta))^kk))/(factorial(kk)));
end
PDF_refractory_ratedgamma2(k+1)=exp(-1*landma*(delta-(k*eta)))*PDF_refractory_ratedgamma2(k+1);


PDF_refractory_ratedgamma(k+1)=PDF_refractory_ratedgamma1(k+1)-PDF_refractory_ratedgamma2(k+1);

end

PDF_ratedgamma(PDF_ratedgamma<0)=0;
PDF_refractory_ratedgamma(PDF_refractory_ratedgamma<0)=0;

%%  Plot
figure;
plot(n_ratedgamma,PDF_ratedgamma,'-o')
hold on
%plot(n_truncated_ratedgamma,PDF_truncated_ratedgamma,'-*')
%hold on
plot(n_refractory_ratedgamma,PDF_refractory_ratedgamma,'--gs')


%legend('ratedgamma','truncated ratedgamma','refractory ratedgamma')
legend('rated-gamma','refractory rated-gamma')
xlabel('{n}')
ylabel('probability')
title('RefractoryGraph')

annotation('textbox',...
    [0.619959914217186 0.729750695801505 0.232981262253403 0.063076923076923],...
    'String','\lambda=100, \Delta=.2, \eta=.002',...
    'FitBoxToText','on');

%\lambda=100, \Delta=.2, \eta=.002
