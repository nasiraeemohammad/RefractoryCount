clc;
close all;

%% parameter

global eta delta P lambda_vector m N

% biological parameter

delta=1*(10^-1);  
m=1;

% setup parameter
prob_lb=.5*10^-2;             %lower bound for probability  pruning
max_pruning=3;

eta_set=[1*(10^-3),2*(10^-3),4*(10^-3)];


lambda_min=5;     
lambda_max_L=20;
lambda_max_H=100;
lambda_max_selected=80;
eta_selected=2*(10^-3);
lambda_max_set=lambda_max_L:5:lambda_max_H;

capacity_varyinglambda=zeros(length(eta_set),length(lambda_max_set)); 
etaindex=0;

for eta=eta_set

etaindex=etaindex+1;

%%  cardinality of channel output's alphabet

N=(delta/eta)-1;

% lambda=lambda_max;
% N=0; 
% CDF=0;
% while (1-abs(CDF))>10^(-15)
%     
% S=0;
% for kk=0:1:N  
% S=S+((lambda^kk)*((delta-(N+1)*eta)^kk))/factorial(kk);    
% end
% CDF=S*exp(delta-(N+1)*eta)
% 
% N=N+1; %cardinality of channel output's alphabet
% end





lambdamaxindex=0;
Opt_PDF_varyinglambda=cell(length(lambda_max_set),1);

 
for lambda_max=lambda_max_set


 
capacity=0;
lambdamaxindex=lambdamaxindex+1;


lambda_vector=lambda_min:1:lambda_max;
lambda_vector_length=size(lambda_vector,2);





%%  pruning

for pruning=1:.5:max_pruning

%channel
P=zeros(lambda_vector_length,N+1);
for i=1:1:lambda_vector_length
for j=0:1:N
    
ConditionalPDF_T1=0;
for IPDF1=0:1:((j+1)*m-1)
ConditionalPDF_T1=ConditionalPDF_T1+(((lambda_vector(i)^IPDF1)*((delta-(j+1)*eta)^IPDF1))/factorial(IPDF1));  
end
ConditionalPDF_T1=ConditionalPDF_T1*exp(-1*lambda_vector(i)*(delta-(j+1)*eta));

ConditionalPDF_T2=0;
for IPDF2=0:1:((j)*m-1)
ConditionalPDF_T2=ConditionalPDF_T2+(((lambda_vector(i)^IPDF2)*((delta-(j)*eta)^IPDF2))/factorial(IPDF2));  
end
ConditionalPDF_T2=ConditionalPDF_T2*exp(-1*lambda_vector(i)*(delta-(j)*eta));

P(i,j+1)= ConditionalPDF_T1-ConditionalPDF_T2;      
end
end

P(P<0)=0; 
P(~isfinite(P))=0;
for i=1:1:lambda_vector_length
 P(i,:)=P(i,:)/sum(P(i,:));   
end 


lb =zeros(1,lambda_vector_length);
ub =ones(1,lambda_vector_length);  

             
Aeq=ones(1,lambda_vector_length);              %the coefficients of linear equalities
beq=1;                                           %the constant terms of linear equalities    
A=[];                                             %the coefficients of linear inequalities
b=[];                                             % the constant terms of linear inequalities
initial_0_Probability=(1/lambda_vector_length)*ones(1,lambda_vector_length);    
options=optimoptions('fmincon','SpecifyObjectiveGradient',true);      
fun=@informationwithgrad;  %Specify Objective Function and Objective Gradient
%nonlcon=@energyconstraintwithgrad;  %Specify Constraint Function and Constraint Gradient
nonlcon=[]; 


%optimazation solver 
[x,fval,exitflag,output,lambda,grad,hessian]=fmincon(fun,initial_0_Probability,A,b,Aeq,beq,lb,ub,nonlcon,options); 



lambda_vector(x<prob_lb)=[];
lambda_vector_length=size(lambda_vector,2);
end



%% cardinality of channel output's alphabet Optimization

maxpoint=size(lambda_vector,2);

lambda_vector_pruning=lambda_vector;
for k=1:1:maxpoint
    
    
subk=nchoosek(lambda_vector_pruning,k);
Nsubk=size(subk,1);   
    
    
    
for kk=1:1:Nsubk
    
lambda_vector=subk(kk,:);     
%channel
P=zeros(k,N+1);
for i=1:1:k
for j=0:1:N

ConditionalPDF_T1=0;
for IPDF1=0:1:((j+1)*m-1)
ConditionalPDF_T1=ConditionalPDF_T1+(((lambda_vector(i)^IPDF1)*((delta-(j+1)*eta)^IPDF1))/factorial(IPDF1));  
end
ConditionalPDF_T1=ConditionalPDF_T1*exp(-1*lambda_vector(i)*(delta-(j+1)*eta));

ConditionalPDF_T2=0;
for IPDF2=0:1:((j)*m-1)
ConditionalPDF_T2=ConditionalPDF_T2+(((lambda_vector(i)^IPDF2)*((delta-(j)*eta)^IPDF2))/factorial(IPDF2));  
end
ConditionalPDF_T2=ConditionalPDF_T2*exp(-1*lambda_vector(i)*(delta-(j)*eta));

P(i,j+1)= ConditionalPDF_T1-ConditionalPDF_T2;    
end
end

P(P<0)=0; 
P(~isfinite(P))=0;
for i=1:1:k
 P(i,:)=P(i,:)/sum(P(i,:));   
end     
    
    
    
 
lb=zeros(1,k);
ub =ones(1,k);    
Aeq=ones(1,k);                               %the coefficients of linear equalities
beq=1;                                             %the constant terms of linear equalities    
A=[];                                            %the coefficients of linear inequalities
b=[];                                            % the constant terms of linear inequalities
initial_0_Probability=(1/k)*ones(1,k);    
options=optimoptions('fmincon','SpecifyObjectiveGradient',true);      
fun=@informationwithgrad;  %Specify Objective Function and Objective Gradient
%nonlcon=@energyconstraintwithgrad_sub;  %Specify Constraint Function and Constraint Gradient
nonlcon=[];   

%optimazation solver 
[x,fval,exitflag,output,lambda,grad,hessian]=fmincon(fun,initial_0_Probability,A,b,Aeq,beq,lb,ub,nonlcon,options);   
    
    

if (-1*fval)>(capacity+0.01)
capacity=(-1*fval);    
capacity_Mass=lambda_vector;
capacity_PDF=x;
capacity_point=size(capacity_Mass,2);
%lagrange_multiplier=lambda.ineqnonlin;
capacity_varyinglambda(etaindex,lambdamaxindex)=capacity;
Opt_PDF_varyinglambda{lambdamaxindex}=[capacity_Mass,capacity_PDF];
end
    
    

end %kk=1:1:Nsubk
    
end %for k=1:1:maxpoint



%plot kkt condition for selected lambda_max :start
if lambda_max==lambda_max_selected  && eta==eta_selected 
    
Energy_lambda_max_selected=((delta*capacity_Mass)./((eta*capacity_Mass)+1))*(capacity_PDF');
    
%Determinde optimal output marginal PMF
P_selected=zeros(capacity_point,N+1);
for i=1:1:capacity_point
for j=0:1:N 
    
ConditionalPDF_T1=0;
for IPDF1=0:1:((j+1)*m-1)
ConditionalPDF_T1=ConditionalPDF_T1+(((capacity_Mass(i)^IPDF1)*((delta-(j+1)*eta)^IPDF1))/factorial(IPDF1));  
end
ConditionalPDF_T1=ConditionalPDF_T1*exp(-1*capacity_Mass(i)*(delta-(j+1)*eta));

ConditionalPDF_T2=0;
for IPDF2=0:1:((j)*m-1)
ConditionalPDF_T2=ConditionalPDF_T2+(((capacity_Mass(i)^IPDF2)*((delta-(j)*eta)^IPDF2))/factorial(IPDF2));  
end
ConditionalPDF_T2=ConditionalPDF_T2*exp(-1*capacity_Mass(i)*(delta-(j)*eta));

P_selected(i,j+1)= ConditionalPDF_T1-ConditionalPDF_T2;    
end
end

P_selected(P_selected<0)=0;
P_selected(~isfinite(P_selected))=0;


for i=1:1:capacity_point
P_selected(i,:)=P_selected(i,:)/sum(P_selected(i,:));   
end

%Determine selected optimal output pdf
y_selected=capacity_PDF*P_selected;
y_selected=y_selected/sum(y_selected);

%Determine selected optimal capacity no per cost
P_selected_capacity_nopercost_LOG=log2(P_selected./y_selected);
P_selected_capacity_nopercost_LOG(~isfinite(P_selected_capacity_nopercost_LOG))=0;

%Determine selected Conditional Mutual Information no per cost
cond_mutual_discrete_selected_nopercost=sum(P_selected.*P_selected_capacity_nopercost_LOG,2);    %is column vector
Information_selected=capacity_PDF*cond_mutual_discrete_selected_nopercost;


%Discretizing lambda for checking KKt condotion
lambdavector_discrete=lambda_min:.05:lambda_max;
lambdavector_size_discrete=size(lambdavector_discrete,2);

P_discrete=zeros(lambdavector_size_discrete,N+1);
for i=1:1:lambdavector_size_discrete
for j=0:1:N
    
ConditionalPDF_T1=0;
for IPDF1=0:1:((j+1)*m-1)
ConditionalPDF_T1=ConditionalPDF_T1+(((lambdavector_discrete(i)^IPDF1)*((delta-(j+1)*eta)^IPDF1))/factorial(IPDF1));  
end
ConditionalPDF_T1=ConditionalPDF_T1*exp(-1*lambdavector_discrete(i)*(delta-(j+1)*eta));

ConditionalPDF_T2=0;
for IPDF2=0:1:((j)*m-1)
ConditionalPDF_T2=ConditionalPDF_T2+(((lambdavector_discrete(i)^IPDF2)*((delta-(j)*eta)^IPDF2))/factorial(IPDF2));  
end
ConditionalPDF_T2=ConditionalPDF_T2*exp(-1*lambdavector_discrete(i)*(delta-(j)*eta));

P_discrete(i,j+1)= ConditionalPDF_T1-ConditionalPDF_T2;   
end
end

P_discrete(P_discrete<0)=0; 
P_discrete(~isfinite(P_discrete)) =0;

for i=1:1:lambdavector_size_discrete
 P_discrete(i,:)=P_discrete(i,:)/sum(P_discrete(i,:));   
end

%Determine KKT condition
P_discrete_LOG=log2(P_discrete./y_selected);
P_discrete_LOG(~isfinite(P_discrete_LOG))=0;

%Conditional Mutual Information
cond_mutual_discrete=sum(P_discrete.*P_discrete_LOG,2);    %is column vector

%Energy discrete
Energy_discrete=((delta*lambdavector_discrete)./((eta*lambdavector_discrete)+1))';          %is column vector

%kkt_condition
kkt_condition=(Information_selected*Energy_discrete)-(Energy_lambda_max_selected*cond_mutual_discrete);


%figure global KTT condition','Capacity-Achieving Distribution'
figure;
plot(lambdavector_discrete,kkt_condition,'LineWidth',.8)
hold on;
stem(capacity_Mass,capacity_PDF,'LineWidth',.8)  
xlabel('{\lambda}')
legend('KTT condition','Capacity-Achieving Distribution','fontsize', 10)

end
%kkt condition plot for selected lambda_max :end





end   %for lambda_max=lambda_max_set


if eta==eta_selected

%% Plot optimal PDF as a function of energybudjet
figure;
 for i = 1:length(lambda_max_set)
    L=length(Opt_PDF_varyinglambda{i})/2;
    Opt_dist_varyinglambda=Opt_PDF_varyinglambda{i};
    for j = 1:L
    pause(.01)
    stem3(lambda_max_set(i),Opt_dist_varyinglambda(1,j),Opt_dist_varyinglambda(1,j+L),'LineWidth',.6);
    hold on;
    end
end
xlabel('b:upper bound on input intensity')
ylabel('Positions of mass points')
zlabel('Probability of mass points')
set(get(gca,'xlabel'),'rotation',20);  %where angle is in degrees
set(get(gca,'ylabel'),'rotation',-10); %where angle is in degrees



end  %end of   %for lambda_max=lambda_max_set




end %end of eta=eta_set

figure;
legend('AutoUpdate','on')
for iii=1:1:length(eta_set)
    
 plot(lambda_max_set,capacity_varyinglambda(iii,:),'DisplayName',[' \eta= ' num2str(eta_set(iii)) ])  
 legend('show'); 
 hold on;   
    
    
end

xlabel('b:upper bound on input intensity')
ylabel('Capacity per Unit cost')









%% Function part


function [f,g] = informationwithgrad(x)

global  P lambda_vector delta eta m


y_Pr=x*P;
y_Pr=y_Pr/sum(y_Pr);

P_LOG=log2(P./y_Pr);
P_LOG(~isfinite(P_LOG))=0;

%determinde  mutul information
condi_mutual=sum(P.*P_LOG,2);          %is column vector
Information=x*condi_mutual;

ConsumedEnergy=((delta*lambda_vector)./((eta*lambda_vector)+m))*(x');

%determinde Objective function
f=-(Information/ConsumedEnergy);     %Objective function
   

  
%gradient function
D_Information=condi_mutual-(1/log(2));
D_ConsumedEnergy=((delta*lambda_vector)./((eta*lambda_vector)+m))';

g= -((D_Information*ConsumedEnergy)-(D_ConsumedEnergy*Information))/(ConsumedEnergy^2);   %garadian of Objective function


 
end
