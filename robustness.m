function re_num=cycle_con_max(A)
re_num=zeros(1,10);   %最大子图节点个数
matrix_len=length(A);    %总的节点个数
for i=1:10
    for k=1:10
    select_num=round(matrix_len*i*0.1);        %失效节点个数
    select_node=zeros(1,1);                 %失效节点记录
    select_node=randperm(matrix_len,select_num);        %选出失效选出的节点
     B=A;
   for j=1:select_num
     B(:,select_node(1,j))=0;
     B(select_node,:)=0;
   end
   C=con_max(B);
  re_num(1,i)=length(C)/matrix_len+re_num(1,i);
    end
    re_num(1,i)= re_num(1,i)/10;
end
savefile1='re_num.mat'; %保存每次剩余的节点数
save(savefile1,'re_num');


