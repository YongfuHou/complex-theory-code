function BAP1(S)                      %求点强分布
T=520;
PB=zeros(1,1);                      
PC=zeros(1,1);
Q=zeros(1,1);
for i=1:T
   % if S(1,i)~=0
 Q(1,i)=round(S(1,i));
  %  end
end
PD=tabulate(Q(:));

for i=1:length(find(PD(:,1)))-1
% PB(1,i)=log10(PD(i+1,1));
  PB(1,i)=PD(i+1,1);
end

for i=1:length(find(PD(:,1)))-1
 %  PC(1,i)=log10(PD(i+1,3)/100);
  PC(1,i)=PD(i+1,3)/100;
end

figure(1)
axis auto;
plot(PB,PC,'ro');
hold on;

savefile1='PD.mat'; %保存数据的文件名
save(savefile1,'PD');








% x=data(:)
% x=sort(x);
% d=diff([x;max(x)+1]);
% count = diff(find([1;d])) ;
% y =[x(find(d)) count]
