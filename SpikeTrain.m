clc
close all

%landma, eta, delta, n
parameter=[
150, .002, .3, 1;
90, .002, .3, 4
];

ISItrain=cell(size(parameter,1),1) ;




for i=1:1:size(parameter,1)

s=0;
theta=1/parameter(i,1);
spiketime=[];
j=1;
while s<=parameter(i,3)
%s=s+gamrnd( n,theta)+eta;
rng('shuffle');
%ISI=random('Gamma',n,theta)+eta
ISI=gamrnd(parameter(i,4),theta)+parameter(i,2);
s=s+ISI;
spiketime(j)=s;
j=j+1;   
end


ISItrain{i}=spiketime;

end


subplot(2,1,1);
value=ones(1,length(ISItrain{1}));
p1=stem(ISItrain{1},value,'Marker','none','Color','blue');
 p1(1).LineWidth = 2;
xlabel('time[sec]') 
xlim([0 delta])
set(gca,'ytick',[])
%axis off
subplot(2,1,2); 
value=ones(1,length(ISItrain{2}));
p2=stem(ISItrain{2},value,'Marker','none','Color','blue');
 p2(1).LineWidth = 2;
xlim([0 delta])
set(gca,'ytick',[])
xlabel('time[sec]')




% value=ones(1,length(spiketime1));
% p=stem(spiketime1,value,'Marker','none','Color','blue')
% xlim([0 delta])
% set(gca,'ytick',[])
% p(1).LineWidth = 1;
% %axis off
% xlabel('time[sec]')