function W_fenbu(W)
T=1020;
Q=zeros(1,1);
sum_t=0;
for i=1:T
    for j=i+1:T
        if W(i,j)~=0
       sum_t=sum_t+1;
       Q(1, sum_t)=W(i,j);
        end
    end
end
for i=1:length(Q)
   % if S(1,i)~=0
 Q(1,i)=round(Q(1,i));
  %  end
end
PD=tabulate(Q(:));
for i=1:length(find(PD(:,1)))-1
%PB(1,i)=log10(PD(i+1,1));
  PB(1,i)=PD(i+1,1);
end

for i=1:length(find(PD(:,1)))-1
   PC(1,i)=log10(PD(i+1,3)/100);
 % PC(1,i)=PD(i+1,3)/100;
end

figure(1)
axis auto;
plot(PB,PC,'*');
hold on;
savefile1='Q.mat'; %保存数据的文件名
save(savefile1,'Q');