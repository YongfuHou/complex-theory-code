function re_num=cycle_con_max(A)
re_num=zeros(1,10);   %�����ͼ�ڵ����
matrix_len=length(A);    %�ܵĽڵ����
for i=1:10
    for k=1:10
    select_num=round(matrix_len*i*0.1);        %ʧЧ�ڵ����
    select_node=zeros(1,1);                 %ʧЧ�ڵ��¼
    select_node=randperm(matrix_len,select_num);        %ѡ��ʧЧѡ���Ľڵ�
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
savefile1='re_num.mat'; %����ÿ��ʣ��Ľڵ���
save(savefile1,'re_num');


