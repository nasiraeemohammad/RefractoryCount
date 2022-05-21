clc
close all

%% parameters

landma=1:2:100;
eta=.004;          %A typical value of the refractory period is 2 ms
delta=.2;
m=1;
n_max=(delta/eta)-1;
mean_refractory_ratedgamma=(n_max+1)+zeros(1,size(landma,2));

for i=1:1:size(landma,2)



for n=0:1:n_max
for k=0:1:(((n+1)*m)-1)
mean_refractory_ratedgamma(i)=mean_refractory_ratedgamma(i)...
-((exp(-1*landma(i)*(delta-((n+1)*eta)))*(landma(i)^k)*((delta-((n+1)*eta))^k))/(factorial(k)));
end
end


end

legend('AutoUpdate','on')
plot(landma,mean_refractory_ratedgamma,'DisplayName',[' \eta= ' num2str(eta) ' , \Delta= ' num2str(delta)  ])
legend('show'); 

xlabel('{\lambda}')
ylabel('E[N|\Lambda=\lambda]')
%title('RefractoryGraph')
