%% -----�����е�n0���ڵ�����翪ʼ�����������������������ӵĻ�������BA�ޱ������----- %%

clear
clc
T=500;%��������ʱ�����ݻ���������ģN��
A=zeros(1,1);%�����ڽӾ���
Edge=zeros(1,1);%���ӱ߾���
C=zeros(1,1);%����ڵ��������õ��ŵ�
Neigb_Record=zeros(1,1);%��¼�ڵ��ھӽڵ�
N_Record=zeros(1,1);%��¼�ڵ��ڽӽڵ�
Neigb_count=zeros(1,1);%��¼�ڵ��ھӽڵ���� 
A_c=zeros(1,1);%��¼���ھӽڵ����������ŵ�
S=zeros(1,1); %�����и��ڵ�ǿ�ȷֲ�
K=zeros(1,1); %�����и��ڵ�ȷֲ�
 W= zeros(1);%��Ȩ����
Power0=zeros(1,1);%N���ڵ㽨�����˺�ĳ�ʼ����
X_Node=zeros(1,1);%��ʼ���ڵ������
Y_Node=zeros(1,1);%��ʼ���ڵ�������
Netnode_record=zeros(1,1);%��¼�����д��ڵĽڵ���
 Dis_twonode=zeros(1,1);%����ڵ��ľ������ 
%% ----------------------------��ʼ������---------------------------- %%
n0=20;%δ����ǰ������ڵ����n0
N_flag=n0;
node_num=n0;%��ǰ�����д��ڵĽڵ�������ÿ��һ����������и���
Netnode_record=[1:n0];%��¼�����д��ڵĽڵ���
L=50;%���񳤿�
cov=200;%�ڵ㸲�Ƿ�Χ��������
Area=500;%�������帲�Ǵ�С
E=500;%�����½�����
Power_down=1/E;%ÿ��ʱ�������½�ֵ
num=Area/L;%��������
S_sum_local=zeros(num,num);%�����ڽڵ�ǿ�ȣ���num*num������
S_sum_whole=0;%���������нڵ�Ľڵ�ǿ��֮��
S_P=zeros(num,num);%�ڵ�ǿ�ȷֲ��ܶȾ���
S_P_Array=zeros(1,num*num);%����ѡ���������
Grid_num=[1:num*num];%����������
Node_num_local=zeros(1,num*num);%��¼ÿ�������нڵ����֮��
Range=100;%�����ڵ��ѡλ�õķ�Χֵ 
M=1000;
C_p=zeros(1,1);%���û���ʹ�õ��ŵ�
Node_c=zeros(1,1);%�ڵ������ŵ��ļ���
C_record=[1:M];%��¼�������ŵ��ı��
%------��ʼ���ڵ�ֲ�������ֲ���-------%
for i=1:n0
    X_Node(1,i)=100+300*rand(1);
    Y_Node(1,i)=100+300*rand(1);
end
%----------��ʼ���ڵ����---------%
for i=1:n0
    for j=1:n0
     Dis_twonode(i,j)=sqrt((X_Node(1,i)-X_Node(1,j))^2+(Y_Node(1,i)-Y_Node(1,j))^2);
    end
end
%----------��ʼ���ڵ�����---------%
for i=1:n0
    Power0(1,i)=1;%�ڵ�����ʱ�ĳ�ʼ����Ϊ1
end
%-----��ʼ�ھ�/�ڽӾ���A(����ڴ��䷶Χ�ڣ���ȫ��ͨ)-------%
for i=1:n0
    for j=1:n0
        if (i~=j&& Dis_twonode(i,j)<=cov)      
           A(i,j)=1; %�ھӾ���
          Edge(i,j)=1; %�ڽӾ���
        else
             A(i,j)=0; %�ھӾ���
              Edge(i,j)=0; %�ڽӾ���
        end
    end
end
%-----��ʼ��Ȩ����-------%
W=Edge;
%------��ʼ����¼�ھӽڵ��------%
CC=0;
for u=1:n0
    for v=1:n0
        if A(u,v)==1
            CC=CC+1;
            Neigb_Record(u,CC)=v;
        end
    end
    CC=0;
end
%------��ʼ����¼�ڽӽڵ��------%
CC=0;
for u=1:n0
    for v=1:n0
        if Edge(u,v)==1
            CC=CC+1;
            N_Record(u,CC)=v;
        end
    end
    CC=0;
end
%------��ʼ���ŵ��ֲ�------%
  a=0.01;        %�ŵ���on-off�ĸ���
b=0.001;        %�ŵ�off-on�ĸ���
H=M;       %�ŵ�����
time_on=zeros(1,H);     %ÿ���ŵ��Ĳ�����ʱ��
time_off=zeros(1,H);     %ÿ���ŵ��Ŀ���ʱ��
cycle=zeros(1,H);       %��¼ÿ���ŵ�������ֵ
time_initial=zeros(1,H);    %��¼ÿ���ŵ���ʼ״ֵ̬
H_state_before=zeros(1,H);       %��¼�ŵ�ǰһʱ�̵�״̬
H_state=zeros(1,H);        %��¼�ŵ���ǰʱ�̵�״̬,�ŵ���״̬��0��ʾ����off״̬���ã�1��ʾ����on״̬������
time_on=exprnd(1/a,1,H);     %����ÿ���ŵ��Ŀ���ʱ��
time_off=exprnd(1/b,1,H);    %����ÿ���ŵ��Ĳ�����ʱ��
cycle=time_on+time_off;       %����ÿ���ŵ�������ֵ
num_migration=zeros(1,T);    %��¼�ŵ��ɳ�Ĭ��Ϊ��Ծ������
for i=1:H
    time_initial(1,i)=cycle(1,i)*rand(1);         
    if  time_initial(1,i)<time_on(1,i)
        H_state(1,i)=1;
    end
end
%------��¼��ʼ���ڵ�ǿ�ȷֲ�-------%
for i=1:n0
   S(1,i)=sum(W(i,:));%����ڵ�ڵ��ǿ
end
for i=1:n0
   if  S(1,i)~=0
       S_sum_whole=S_sum_whole+S(1,i);
   end
end

%% ---------------------��������+���������ڽڵ�ǿ��+������������ڽڵ���֮��---------------- %%
 for i=1:num %��num*num������ ��
        Y_D=(i-1)*L;%ȷ�������귶Χ
        Y_U=Y_D+L;
        for j=1:num  %��u 
            X_L=(j-1)*L;%ȷ�������귶Χ
            X_R=X_L+L;
         %����һ�������ڵĽڵ�ǿ��
  for u=1:N_flag
                if ((X_L<X_Node(1,u)) && (X_Node(1,u)<=X_R))
                    if ((Y_D<Y_Node(1,u)) && (Y_Node(1,u)<=Y_U))
                        if  S(1,u)~=0
                           S_sum_local(i,j)=S_sum_local(i,j)+S(1,u);%�����ڽڵ㣬�ڵ�ǿ�����
                        end
                        Node_num_local(1,((i-1)*num+j))=Node_num_local(1,((i-1)*num+j))+1;%��������нڵ����
                    end
                end          
   end
            if Node_num_local(1,((i-1)*num+j))==0
                S_P(i,j)=0;
            else        
      %��������ڵ�ȷֲ��ܶȸ���
               S_P(i,j)=S_sum_local(i,j)/S_sum_whole;
            end
            %��������ʼ�¼��һ����������
            S_P_Array(1,((i-1)*num+j))=S_P(i,j);
        end  
 end
 
%  node_n=zeros(1,1);
%  edge_n=zeros(1,1);
%  for aaa=1:10
%% ---------*****************---------��ѭ��--------*****************------------ %%
%  T=100*aaa;
for ttt=1:20
    ttt
    N_flag=N_flag+1; %��Ƿ�   
    node_num=node_num+1;%�����½ڵ��������ʵ�ʽڵ���������ɾ���㡢�߹��̣�
    Power0(1,N_flag)=1;%�����ӽڵ�������ʼ��
    %��¼�����������ڵ���
    Netnode_record(1,N_flag)=N_flag; 
    %% -----ÿ��ʱ�������нڵ�������1/E+�����ӽڵ�������ʼ��------ %%
    for i=1:(N_flag-1)
        if Power0(1,i)~=0
            Power0(1,i)=Power0(1,i)-1/E;
        end
    end
    %% -Part 1-ÿ��ʱ��������������һ���ڵ㣬ȷ����������---------- %%
    flag=1;%flag=1��ʾ�ڵ����������⣬����ѡ�����ꣻ��֮�����������ڣ���������ѡ��
    tab=1;%tab=1��ʾ�ڵ㸲�Ƿ�Χ���޽ڵ㣬����ѡ�����ꣻ��֮�����Ƿ�Χ���нڵ㣬��������ѡ��   %20130206�޸�
    while (flag || tab)
        %������ѡ��һ������
        sel=randsrc(1,1,[Grid_num(1,:);S_P_Array(1,:)]);
        %ѡ��������ڸ����������λ�ô�����һ���ڵ�
        Node_num_local(1,sel)=Node_num_local(1,sel)+1;
        %ȡ�������������λ��ֵ
        i=fix(sel/10)+1;
        j=sel-(i-1)*10;
        %ȷ��������������귶Χ
        X_center=(L/2+(j-1)*L);%ȷ���������ĵ�����
        Y_center=(L/2+(i-1)*L);
     %ȷ�������귶Χ
        X_sel_L=X_center-(Range/2);
        X_sel_R=X_center+(Range/2);
        %ȷ�������귶Χ
        Y_sel_D=Y_center-(Range/2);
        Y_sel_U=Y_center+(Range/2);
    
        %part1���ж��Ƿ��ڸ��Ƿ�Χ��
        X_Node(1,N_flag)=X_sel_L+rand(1)*Range;
        Y_Node(1,N_flag)=Y_sel_D+rand(1)*Range;
        if (0<=X_Node(1,N_flag) && X_Node(1,N_flag)<=500)
            if (0<=Y_Node(1,N_flag) && Y_Node(1,N_flag)<=500)
               flag=0;%����Ҫ��
            end         
        end
      %part2���жϽڵ㸲�Ƿ�Χ���Ƿ��������ڵ�
        local_neigber_count=0;%��¼�����ӽڵ�ĸ��Ƿ�Χ�ڽڵ����
        X_min=X_Node(1,N_flag)-cov/2; %ȷ�������ڵ�ĸ�������
        X_max=X_Node(1,N_flag)+cov/2; 
        Y_min=Y_Node(1,N_flag)-cov/2;
        Y_max=Y_Node(1,N_flag)+cov/2;
        for j=1:N_flag
            if N_flag~=j           
                if ((X_min<=X_Node(1,j)) &&(X_Node(1,j)<=X_max))
                    if ((Y_min<=Y_Node(1,j)) && (Y_Node(1,j)<=Y_max))
                        A(N_flag,j)=1;%����ͼ
                        A(j,N_flag)=1;
                        local_neigber_count=local_neigber_count+1;%�����ӽڵ�ĸ��Ƿ�Χ�ڽڵ����                   
                    end
                end
            end
        end
        if local_neigber_count~=0
            tab=0;%��������
        end
    end
    %--���������ڵ��ھӽڵ�ż�¼����--%
    CC=0;
    Neigb_Record(:,:)=0;
    for u=1:N_flag
        for v=1:N_flag
            if A(u,v)==1
                CC=CC+1;
                Neigb_Record(u,CC)=v;
            end
        end
        CC=0;
    end
%% --Part 2--ÿ��ʱ���������ӽڵ����������Ͻڵ�����m�����ڽڵ㸲�Ƿ�Χ�����ߣ�------ %%
 m=12;%���ӱ������̶�һ��ֵ��      �ɱ����
 p=0.3;% ��ȥ��һ�����⣬ʣ��ı��Ը���p��������      �ɱ����
 w0=1;%����ߵĳ�ʼȨֵ
 deta=1;%ÿ�α����ӣ��ڵ��ǿ���ӵ���ֵ
 alf=0.5;%   �������ŵ��Ŀɵ�����
mt=(m-1)*p;%�������ӵıߵĸ���
mp=(m-mt);%�������ӵıߵĸ���
 Pt=1;%(W)���书��30dBm    
 M=1000;%���û��Լ��ŵ�����
bandwidth_sum=40*10^7;%ϵͳ�ܴ���10MHz
W0=bandwidth_sum/M;%ÿ�����ز�����
 WGN=10^(-17.4)*W0*10^-3;%�������ʣ�W��  
An=zeros(1,1);%��¼�ڵ������ŵ�����
 if local_neigber_count~=0 %�ھӽڵ�����Ϊ0
        if (1<=local_neigber_count && local_neigber_count<=m)%��ȫ��
            for i=1:local_neigber_count
                out=Neigb_Record(N_flag,i);
                Edge(N_flag,out)=1; %���ӱ߾���
                Edge(out,N_flag)=1; 
            end
        else
            %��¼ѡ�������ߵĽڵ��
            selout_node=zeros(1,m);
            %%========step4:�������ڵ㸲�Ƿ�Χ�ڽڵ��������ڵ������Լ�Ǳ����·����
            Dis_TwoNode=zeros(1,local_neigber_count);%���������ڵ㵽���Ƿ�Χ�ڽڵ�ľ������
            PL_TwoNode=zeros(1,local_neigber_count);%���������ڵ㵽���Ƿ�Χ�ڽڵ���ŵ�˥����󣡣���
            SNR_TwoNode=zeros(1,local_neigber_count);%���������ڵ㵽���Ƿ�Χ�ڽڵ����·����Ⱦ��󣡣���
            Cap_TwoNode=zeros(1,local_neigber_count);%���������ڵ㵽���Ƿ�Χ�ڽڵ����·�������󣡣�����
            for i=1:local_neigber_count
                Dis_TwoNode(1,i)=sqrt((X_Node(1,N_flag)-X_Node(1,Neigb_Record(N_flag,i)))^2+(Y_Node(1,N_flag)-Y_Node(1,Neigb_Record(N_flag,i)))^2);
                PL_TwoNode(1,i)=30.6+36.7*log10(Dis_TwoNode(1,i));%��߶�˥��������йأ�С�߶�˥����������ֲ�
                SNR_TwoNode(1,i)=Pt*PL_TwoNode(1,i)./WGN;
                for j=1:M
                   if H_state(1,j)~=1
                     Cap_TwoNode(1,i)=Cap_TwoNode(1,i)+W0*log2(1+SNR_TwoNode(1,i)); 
                   else
                        Cap_TwoNode(1,i)=Cap_TwoNode(1,i)+0;
                   end
                end               
            end
           %%=========step5:�����ӱߵĸ���
             Edge_P_Array=zeros(1,local_neigber_count); %�ڵ㱻����ѡ�����ߵĸ���
             Edge_P_mid=zeros(1,local_neigber_count);%������ѡ������м�����������������ŵ�����
             Neigb_Choose=zeros(1,local_neigber_count);%�����ڵ��ھӽڵ�
              Edge_P_mid_sum=0;%������ѡ������м�������
              %���м���ʲ���
            for i=1:local_neigber_count
               Edge_P_mid(1,i)=(Power0(1,Neigb_Record(N_flag,i)).^alf)*(Cap_TwoNode(1,i).^(1-alf));
                Edge_P_mid_sum=Edge_P_mid_sum+(Edge_P_mid(1,i)*S(1,Neigb_Record(N_flag,i)));
            end
             %�����ո���
            for i=1:local_neigber_count
                Edge_P_Array(1,i)=(Edge_P_mid(1,i)*S(1,Neigb_Record(N_flag,i)))/Edge_P_mid_sum;
                Neigb_Choose(1,i)=Neigb_Record(N_flag,i);
            end
             %%===========step6:���������ӱ�
             %��������ѡ��ı�
              mp_num=1;%��ʼֵ
            while mp_num<=mp %���ӱ�ֱ�����ӵ�mp��Ϊֹ
                out=randsrc(1,1,[Neigb_Choose(1,:);Edge_P_Array(1,:)]);%
                %�ж�ѡ���Ľڵ��Ƿ��ظ�
                if length(find(selout_node==out))==0 %���ظ���
                    selout_node(1,mp_num)=out;%��¼
                    Edge(N_flag,out)=1; %���ӱ߾���
                    Edge(out,N_flag)=1;
                    mp_num=mp_num+1;%������һ
                end
            end
            %�����������ӵı�
            selout_node2=zeros(1,1);
             mt_possible_num=0;
              mt_possible_num_p_sum=0;
               mt_possible_num_p=zeros(1,1);
                 for i=1:mp
                  out=selout_node(1,i);   
                   for j=1:length(find( N_Record(out,:)))
                       if sqrt((X_Node(1,N_flag)-X_Node(1,N_Record(out,j)))^2+(Y_Node(1,N_flag)-Y_Node(1,N_Record(out,j)))^2)<=cov
                           mt_possible_num=mt_possible_num+1;
                            selout_node2(1,mt_possible_num)=N_Record(out,j);
                            mt_possible_num_p_sum= mt_possible_num_p_sum+W(out,N_Record(out,j))/S(1,out);  %Ϊ�˹�һ���������ֵ���ۺ�
                       end
                   end
                 end
                 DD=0;
                for i=1:mp
                  out=selout_node(1,i);   
                   for j=1:length(find( N_Record(out,:)))
                       if sqrt((X_Node(1,N_flag)-X_Node(1,N_Record(out,j)))^2+(Y_Node(1,N_flag)-Y_Node(1,N_Record(out,j)))^2)<=cov
                           DD=DD+1;
                            mt_possible_num_p(1,DD)= (W(out,N_Record(out,j))/S(1,out))/mt_possible_num_p_sum;
                       end
                   end
                 end
                 
            if  mt_possible_num<=mt;%����ѡ��Ľڵ�����ӽڵ����С��Ҫ�������ӵĽڵ��������ȫ����
                  for i=1: length(find(selout_node2))
                     out=selout_node2(1,i);                      
                       Edge(N_flag,out)=1; 
                       Edge(out,N_flag)=1; 
                  end
            else  
              mt_num=1;
                   while mt_num<=mt %���ӱ�ֱ�����ӵ�mt��Ϊֹ
                   out=randsrc(1,1,[selout_node2(1,:); mt_possible_num_p(1,:)]);
                %�ж�ѡ���Ľڵ��Ƿ��ظ�
                  if length(find(selout_node==out))==0 %���ظ���
                    selout_node(1,mp_num)=out;%��¼
                    Edge(N_flag,out)=1; %���ӱ߾���
                    Edge(out,N_flag)=1;
                    mp_num=mp_num+1;%������һ
                    mt_num=mt_num+1;
                     end
                   end
           end
        end
 else
       Edge(N_flag,:)=0; %�����ڵ����@@@@@@@@
        Edge(:,N_flag)=0; 
 end


  %% --����1���ڵ��m���߹��̺󣬸�������ڵ�ǿ�����Ȩ--- %%
 for i=1:N_flag-1
       if Edge(N_flag,i)==1
           W(N_flag,i)=w0;
            W(i,N_flag)= W(N_flag,i);
       end
 end
 for i=1:N_flag-1
     if Edge(N_flag,i)==1
         for j=1:length(find( N_Record(i,:)))
           W(i,N_Record(i,j))= W(i,N_Record(i,j))+(W(i,N_Record(i,j))./S(1,i))*deta;
           W(N_Record(i,j),i)=W(i,N_Record(i,j));
         end
     end
 end
for i=1:N_flag
   S(1,i)=sum(W(i,:));%����ڵ�ڵ��ǿ
end
%% --�����ڽӽڵ��--- %%
CC=0;
for u=1:N_flag
    for v=1:N_flag
        if Edge(u,v)==1
            CC=CC+1;
            N_Record(u,CC)=v;
        end
    end
    CC=0;
end
 %% --Part 3--ÿ��ʱ���Է����Ÿ���ɾ������ɾ������n������ȫ����Χѡ��------ %%
 n=4;%�趨ɾ����·��Ŀ
    %��¼ѡ����ɾ����·һ�˽ڵ��node1&node2
    selout_node1=zeros(1,n);
    selout_node2=zeros(1,n);
 %%====step1:��ѡ��ɾ����·��1�˵ĸ���
    Node1_P_Array=zeros(1,1);%��̬����ѡ��ɾ����·��1�˸���
    Node1_P_mid=zeros(1,1);%ѡ��ɾ����·��1�˸����м�ֵ
    Neigb1_Choose=zeros(1,1);%��1�˽ڵ����߽ڵ�
    Node1_P_mid_sum=0;%ѡ��ɾ����·��1�˸����м�ֵ���
     %���ѡ�ڵ㼯��
    k=0;
    for i=1:N_flag
        if Netnode_record(1,i)~=0
            k=k+1;
            Neigb1_Choose(1,k)=Netnode_record(1,i);
        end
    end
    %--��Node1�м����
    for i=1:k
        Node1_P_mid(1,i)=(Power0(1,Neigb1_Choose(1,i))*S(1,Neigb1_Choose(1,i)))^-1;
        Node1_P_mid_sum=Node1_P_mid_sum+Node1_P_mid(1,i);
    end
 %---��Node1���ո���
    for i=1:k
        Node1_P_Array(1,i)=Node1_P_mid(1,i)./Node1_P_mid_sum;
    end
  %%===step2:������ɾ����
   n_num=1;%��ʼֵ    
   while n_num<=n %ɾ����ֱ��n��Ϊֹ
     out1=randsrc(1,1,[Neigb1_Choose(1,:); Node1_P_Array(1,:)]);%**ONE:������ѡ����·��һ�˽ڵ�node1
    node_linked_count=sum(Edge(out1,:));%��¼ѡ���ڵ�1�����߽ڵ����
     if node_linked_count~=0
        selout_node1(1,n_num)=out1;%��¼
         %�����߽ڵ����
                linked_array=zeros(1,node_linked_count);
                t=0;
                for i=1:N_flag
                    if Edge(out1,i)~=0
                        t=t+1;%����
                        linked_array(1,t)=i;
                    end
                end
      %%========step3:��node1�ڵ����߽ڵ���node1�ڵ������Լ�Ǳ����·����
            Dis_TwoNode1=zeros(1,node_linked_count);%������out1���ӽڵ�ľ���
            PL_TwoNode1=zeros(1,node_linked_count);%������out1�ڵ���ŵ�˥����󣡣���
            SNR_TwoNode1=zeros(1,node_linked_count);%������out1�ڵ����·����Ⱦ��󣡣���
            Cap_TwoNode1=zeros(1,node_linked_count);%���������ڵ㵽���Ƿ�Χ�ڽڵ����·�������󣡣�����
            for i=1:node_linked_count
               Dis_TwoNode1(1,i)=sqrt((X_Node(1,out1)-X_Node(1,linked_array(1,i)))^2+(Y_Node(1,out1)-Y_Node(1,linked_array(1,i)))^2);
               PL_TwoNode1(1,i)=30.6+36.7*log10(Dis_TwoNode1(1,i));%��߶�˥��������йأ�С�߶�˥����������ֲ�
               SNR_TwoNode1(1,i)=Pt*PL_TwoNode1(1,i)./WGN;
                for j=1:M
                   if H_state(1,j)==0
                      Cap_TwoNode1(1,i)=Cap_TwoNode1(1,i)+W0*log2(1+SNR_TwoNode1(1,i));
                   else
                        Cap_TwoNode1(1,i)=Cap_TwoNode1(1,i)+0;
                   end
                end               
            end
             %%=========step4:��ɾ���ߵĸ���
                Link_P_Array=zeros(1,node_linked_count);%��̬����ɾ����ѡ�����
                Link_P_mid=zeros(1,node_linked_count);%��ɾ��ѡ������м�ֵ
                %Link_Choose=zeros(1,node_linked_count);%�����ڵ��ھӽڵ�
                Link_P_mid_sum=0;%������ѡ������м�ֵ���
           %���м����
                for i=1:node_linked_count
                    Link_P_mid(1,i)=(Power0(1,linked_array(1,i)).^alf)*(Cap_TwoNode1(1,i).^(1-alf));
                    Link_P_mid_sum=Link_P_mid_sum+(Link_P_mid(1,i)*S(1,linked_array(1,i)))^-1;
                end
        %�����ո���
                for i=1:node_linked_count
                    Link_P_Array(1,i)=((Link_P_mid(1,i)*S(1,linked_array(1,i)))^-1)/Link_P_mid_sum;
                end
        %%===========step5:������ɾ��һ����node2��
                out2=randsrc(1,1,[linked_array(1,:);Link_P_Array(1,:)]);%������ѡ��ɾ���ߵ���һ�˽ڵ�node2
                selout_node2(1,n_num)=out2;%��¼
                Edge(out1,out2)=0; %���ӱ߾���ɾ��������
                Edge(out2,out1)=0;
                W(out1,out2)=0;
                W(out2,out1)=0;
 %%%%%%%%%%%%%%%%%%%%%%%%%������Ҫ���ڵĲ���%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
 n_num=n_num+1;%ɾ��������һ
     end 
   end
    %% -----------���µ�ǿ----------- %%
  for i=1:N_flag
   S(1,i)=sum(W(i,:));%����ڵ�ڵ��ǿ
  end
%-------------�����ڽӽڵ��----------------%

CC=0;
for u=1:N_flag
    for v=1:N_flag
        if Edge(u,v)==1
            CC=CC+1;
            N_Record(u,CC)=v;
        end
    end
    CC=0;
end
  %% -----------------������ʣ����------------------ %%
sum_edge1=0;
for i=1:N_flag
   sum_edge1=sum_edge1+sum(Edge(i,:))/2;
end
   %% --Part4:ÿһʱ�����ڵ����ŵ��Ը���AP��ռ�ã��ڵ����Խ����ŵ��л��ĸ���ΪCP -- %%
    H_state_before=H_state;           %�����ڵ��ŵ�״̬��ֵ��before
    for j=1:H
        time_initial(1,j)=rem(time_initial(1,j)+1,cycle(1,j));              %�����ŵ�״̬
     if  time_initial(1,j)<=time_on(1,j)
        H_state(1,j)=1;
     end
      if  time_initial(1,j)>time_on(1,j)
        H_state(1,j)=0;
      end
    end
 num_migration(1,ttt)= numel( find( H_state_before-H_state==-1));     %��¼�ŵ�״̬��0��Ϊ1����Ŀ
 AP=num_migration(1,ttt)/sum_edge1 ;                   %���û��Ļ�Ծ��
   % AP=num_migration(1,ttt)/M ;
   cp=0.5;
   for i=1:N_flag
       for j=i+1:N_flag
          if Edge(i,j)==1
              Edge(i,j)=randsrc(1,1,[0 1;AP 1-AP]);
               Edge(j,i)= Edge(i,j);
               if  Edge(j,i)==0
                   W(i,j)=0;
                   W(j,i)=0;
               end
          end
       end
   end
   % -----------���µ�ǿ----------- %%
  for i=1:N_flag
   S(1,i)=sum(W(i,:));%����ڵ�ڵ��ǿ
  end
  % -----------���½ڵ��----------- %%
  for i=1:N_flag
   K(1,i)=sum(Edge(i,:));%����ڵ�ڵ��ǿ
  end
  %-------------�����ڽӽڵ��----------------%

CC=0;
for u=1:N_flag
    for v=1:N_flag
        if Edge(u,v)==1
            CC=CC+1;
            N_Record(u,CC)=v;
        end
    end
    CC=0;
end
  %% -----------Part5: �ڵ��Ϊ0�Ĺ����ڵ�ɾ��----------- %%
    k=0;
    for i=1:N_flag
       if  S(1,i)==0         
           A(i,:)=0;%�����ڽӾ���
           A(:,i)=0;           
          % Power0(1,i)=0;%N���ڵ㽨�����˺�ĳ�ʼ����
           X_Node(1,i)=0;%�ڵ���������
           Y_Node(1,i)=0;%�ڵ����������
           Netnode_record(1,i)=0;%��¼�����д��ڵĽڵ������  
           k=k+1;%ȫ���ڵ�������
       end
    end
    node_num=N_flag-k;
    
end
 %% -----------------�ų������ڵ㣬�����ߵ����Ӿ���------------------ %%
  node_final=zeros(1,node_num);
  node_tmp=0;
  for i=1:N_flag
      if sum(Edge(i,:))~=0;
          node_tmp=node_tmp+1;
          node_final(1,node_tmp)=i;               %%��ѡ���������Ľڵ�
      end
  end
   Edge_f=zeros(1,1);
  for i=1:(node_num-1) 
     for j=i+1:node_num
        if Edge( node_final(1,i), node_final(1,j))~=0              %%  �����������ڵ�֮��ı߾���
           Edge_f(i,j)=1;
           Edge_f(j,i)=1; 
        end
     end
  end  
%% -----------------������ʣ����------------------ %%
sum_edge=0;
for i=1:N_flag
   sum_edge=sum_edge+sum(Edge(i,:))/2;
end
%% -----------------�����ܵ�ǿ------------------ %%
sum_s=0;
for i=1:N_flag
   sum_s=sum_s+sum(W(i,:));
end
%% - ----------------ģ�⹥�����ɾ��һ���ڵ�------------------ %%
out_s=randsrc(1,1,1:node_num);
Edge_s=Edge_f;
Edge_s(:,out_s)=[];
Edge_s(out_s,:)=[];
%% - ----------------ģ����⹥��ɾ��һ���ڵ�------------------ %%
[zd,wz]=max(S);
out_k=find(node_final==wz);
Edge_k=Edge_f;
Edge_k(:,out_k)=[];
Edge_k(out_k,:)=[];
%% -----------------������������ͼ------------------ %%
figure(1)
axis auto;
X_tmp=X_Node(X_Node~=0);
Y_tmp=Y_Node(Y_Node~=0);
plot(X_tmp,Y_tmp,'ro','MarkerEdgeColor','g','MarkerFaceColor','r','markersize',8);
hold on;
for i=1:(node_num-1) 
     for j=i+1:node_num
        if Edge( node_final(1,i), node_final(1,j))~=0 
             line([X_Node(1, node_final(1,i)),X_Node(1, node_final(1,j))],[Y_Node(1, node_final(1,i)),Y_Node(1, node_final(1,j))],'LineStyle','-','linewidth',1.2); 
             hold on;      %% ����BA�ޱ������ͼ
         end
     end
end
 axis equal;
grid on
hold off 
%  %% -----------���ݴ洢����------------ %%
%  
% eval(['save S',num2str(T)]); %%%����ÿ��ѭ���ĵ�ǿ�ֲ�
% saveas(gca,num2str(T),'jpg');   %%%����ÿ��ѭ��������ͼƬ
% close(gcf);
% node_n(1,aaa)= node_num;    %%%%  ��¼ÿ��ʣ��Ľڵ���
%  edge_n(1,aaa)=sum_edge;
% eval(['save  Edge_f',num2str(T)]); %%%����ÿ��ѭ�������Ӿ���
%  end
% savefile1='node_n.mat'; %����ÿ��ʣ��Ľڵ���
% save(savefile1,'node_n');
% savefile1='edge_n.mat'; %����ÿ��ʣ����ܱ���
% save(savefile1,'edge_n');



