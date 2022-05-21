clc
close all

%% parameters

landma=2:2:100;

eta=.02;          %A typical value of the refractory period is 2 ms
delta=.2;
n_max=(delta/eta)-1;


Samplesize=10000;

Mode=zeros(1,size(landma,2));
poisson_Mode=zeros(1,size(landma,2));
Appro_Mode=zeros(1,size(landma,2));


h = waitbar(0,'Please wait...');


for i=1:1:size(landma,2)

%%
LookupTable(1)=exp(-1*landma(i)*(delta-eta));
numberofspiking=zeros(1,n_max+1);
Samples=[]; 

length=1;


%% Generating random variates
for ii=1:1:Samplesize
    
    
rng('shuffle');
USample= rand;
 
j=1;    
while USample>LookupTable(j)

if j==length
     
 LookupTable(j+1)=0;   
  for jj=0:1:j   
LookupTable(j+1)=LookupTable(j+1)+ ...
   (((landma(i)^jj)/factorial(jj))*((delta-((j+1)*eta))^jj));
  end  
LookupTable(j+1)=LookupTable(j+1)*exp(-1*landma(i)*(delta-((j+1)*eta)));   
    
length=length+1;  
end  
 
j=j+1;    
end   

Samples=[Samples,j-1]; 
numberofspiking(j)=numberofspiking(j)+1;



    
   
end

[A,B]=max(numberofspiking);
probabiltyofmode=A/Samplesize;
Mode(i)=B-1;
poisson_Mode(i)=landma(i)*delta;
Appro_Mode(i)=(landma(i)*exp(landma(i)*eta)*(delta-eta))/(1+(landma(i)*eta*exp(landma(i)*eta)));

waitbar(i/size(landma,2)) 


end
close(h)

filename = 'd:\RefractoryMode.mat';
save(filename)


figure
plot(landma,Mode,'-o')
hold on
plot(landma,floor(poisson_Mode),'-*')
hold on
plot(landma,floor(Appro_Mode),'--gs')

legend('Monte Carlo Mode','Poisson Mode:\Delta*\lambda','Our approximate Mode')
xlabel('{\lambda}')

