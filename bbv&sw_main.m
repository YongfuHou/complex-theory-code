%% -----从已有的n0个节点的网络开始，采用增长机制与优先连接的机制生成BA无标度网络----- %%

clear
clc
T=500;%网络总总时步（演化后的网络规模N）
A=zeros(1,1);%网络邻接矩阵
Edge=zeros(1,1);%连接边矩阵
C=zeros(1,1);%网络节点连接所用的信道
Neigb_Record=zeros(1,1);%记录节点邻居节点
N_Record=zeros(1,1);%记录节点邻接节点
Neigb_count=zeros(1,1);%记录节点邻居节点个数 
A_c=zeros(1,1);%记录与邻居节点相连所用信道
S=zeros(1,1); %网络中各节点强度分布
K=zeros(1,1); %网络中各节点度分布
 W= zeros(1);%边权矩阵
Power0=zeros(1,1);%N个节点建立拓扑后的初始能量
X_Node=zeros(1,1);%初始化节点横坐标
Y_Node=zeros(1,1);%初始化节点纵坐标
Netnode_record=zeros(1,1);%记录网络中存在的节点标号
 Dis_twonode=zeros(1,1);%定义节点间的距离变量 
%% ----------------------------初始化变量---------------------------- %%
n0=20;%未增长前的网络节点个数n0
N_flag=n0;
node_num=n0;%当前网络中存在的节点总数，每做一个动作后进行更新
Netnode_record=[1:n0];%记录网络中存在的节点标号
L=50;%网格长宽
cov=200;%节点覆盖范围！！！！
Area=500;%网络整体覆盖大小
E=500;%能量下降幅度
Power_down=1/E;%每个时步能量下降值
num=Area/L;%网格数量
S_sum_local=zeros(num,num);%网格内节点强度，共num*num个网格
S_sum_whole=0;%网络中所有节点的节点强度之和
S_P=zeros(num,num);%节点强度分布密度矩阵
S_P_Array=zeros(1,num*num);%网格选择概率向量
Grid_num=[1:num*num];%网格编号向量
Node_num_local=zeros(1,num*num);%记录每个网格中节点个数之和
Range=100;%新增节点可选位置的范围值 
M=1000;
C_p=zeros(1,1);%主用户所使用的信道
Node_c=zeros(1,1);%节点所用信道的集合
C_record=[1:M];%记录网络中信道的编号
%------初始化节点分布（随机分布）-------%
for i=1:n0
    X_Node(1,i)=100+300*rand(1);
    Y_Node(1,i)=100+300*rand(1);
end
%----------初始化节点距离---------%
for i=1:n0
    for j=1:n0
     Dis_twonode(i,j)=sqrt((X_Node(1,i)-X_Node(1,j))^2+(Y_Node(1,i)-Y_Node(1,j))^2);
    end
end
%----------初始化节点能量---------%
for i=1:n0
    Power0(1,i)=1;%节点入网时的初始能量为1
end
%-----初始邻居/邻接矩阵A(如果在传输范围内，则全联通)-------%
for i=1:n0
    for j=1:n0
        if (i~=j&& Dis_twonode(i,j)<=cov)      
           A(i,j)=1; %邻居矩阵
          Edge(i,j)=1; %邻接矩阵
        else
             A(i,j)=0; %邻居矩阵
              Edge(i,j)=0; %邻接矩阵
        end
    end
end
%-----初始边权矩阵-------%
W=Edge;
%------初始化记录邻居节点号------%
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
%------初始化记录邻接节点号------%
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
%------初始化信道分布------%
  a=0.01;        %信道的on-off的概率
b=0.001;        %信道off-on的概率
H=M;       %信道总数
time_on=zeros(1,H);     %每个信道的不可用时间
time_off=zeros(1,H);     %每个信道的可用时间
cycle=zeros(1,H);       %记录每个信道的周期值
time_initial=zeros(1,H);    %记录每个信道初始状态值
H_state_before=zeros(1,H);       %记录信道前一时刻的状态
H_state=zeros(1,H);        %记录信道当前时刻的状态,信道的状态，0表示处于off状态可用，1表示处于on状态不可用
time_on=exprnd(1/a,1,H);     %产生每个信道的可用时间
time_off=exprnd(1/b,1,H);    %产生每个信道的不可用时间
cycle=time_on+time_off;       %产生每个信道的周期值
num_migration=zeros(1,T);    %记录信道由沉默变为活跃的数量
for i=1:H
    time_initial(1,i)=cycle(1,i)*rand(1);         
    if  time_initial(1,i)<time_on(1,i)
        H_state(1,i)=1;
    end
end
%------记录初始化节点强度分布-------%
for i=1:n0
   S(1,i)=sum(W(i,:));%求各节点节点点强
end
for i=1:n0
   if  S(1,i)~=0
       S_sum_whole=S_sum_whole+S(1,i);
   end
end

%% ---------------------划分网格+计算网格内节点强度+计算各个网格内节点数之和---------------- %%
 for i=1:num %共num*num个网格 行
        Y_D=(i-1)*L;%确定纵坐标范围
        Y_U=Y_D+L;
        for j=1:num  %列u 
            X_L=(j-1)*L;%确定横坐标范围
            X_R=X_L+L;
         %计算一个网格内的节点强度
  for u=1:N_flag
                if ((X_L<X_Node(1,u)) && (X_Node(1,u)<=X_R))
                    if ((Y_D<Y_Node(1,u)) && (Y_Node(1,u)<=Y_U))
                        if  S(1,u)~=0
                           S_sum_local(i,j)=S_sum_local(i,j)+S(1,u);%网格内节点，节点强度求和
                        end
                        Node_num_local(1,((i-1)*num+j))=Node_num_local(1,((i-1)*num+j))+1;%求该网格中节点个数
                    end
                end          
   end
            if Node_num_local(1,((i-1)*num+j))==0
                S_P(i,j)=0;
            else        
      %计算网格节点度分布密度概率
               S_P(i,j)=S_sum_local(i,j)/S_sum_whole;
            end
            %将网格概率记录在一个行向量内
            S_P_Array(1,((i-1)*num+j))=S_P(i,j);
        end  
 end
 
%  node_n=zeros(1,1);
%  edge_n=zeros(1,1);
%  for aaa=1:10
%% ---------*****************---------主循环--------*****************------------ %%
%  T=100*aaa;
for ttt=1:20
    ttt
    N_flag=N_flag+1; %标记符   
    node_num=node_num+1;%增加新节点后网络中实际节点数（包括删除点、边过程）
    Power0(1,N_flag)=1;%新增加节点能量初始化
    %记录网络中新增节点标号
    Netnode_record(1,N_flag)=N_flag; 
    %% -----每个时步网络中节点能量减1/E+新增加节点能量初始化------ %%
    for i=1:(N_flag-1)
        if Power0(1,i)~=0
            Power0(1,i)=Power0(1,i)-1/E;
        end
    end
    %% -Part 1-每个时步向网络中增加一个节点，确定横纵坐标---------- %%
    flag=1;%flag=1表示节点落在区域外，重新选择坐标；反之，落在区域内，结束坐标选择
    tab=1;%tab=1表示节点覆盖范围内无节点，重新选择坐标；反之，覆盖范围内有节点，结束坐标选择   %20130206修改
    while (flag || tab)
        %依概率选择一个网格
        sel=randsrc(1,1,[Grid_num(1,:);S_P_Array(1,:)]);
        %选出网格后，在该网格中随机位置处增加一个节点
        Node_num_local(1,sel)=Node_num_local(1,sel)+1;
        %取得网格横纵坐标位置值
        i=fix(sel/10)+1;
        j=sel-(i-1)*10;
        %确定该网格横纵坐标范围
        X_center=(L/2+(j-1)*L);%确定网格中心点坐标
        Y_center=(L/2+(i-1)*L);
     %确定横坐标范围
        X_sel_L=X_center-(Range/2);
        X_sel_R=X_center+(Range/2);
        %确定纵坐标范围
        Y_sel_D=Y_center-(Range/2);
        Y_sel_U=Y_center+(Range/2);
    
        %part1：判断是否在覆盖范围内
        X_Node(1,N_flag)=X_sel_L+rand(1)*Range;
        Y_Node(1,N_flag)=Y_sel_D+rand(1)*Range;
        if (0<=X_Node(1,N_flag) && X_Node(1,N_flag)<=500)
            if (0<=Y_Node(1,N_flag) && Y_Node(1,N_flag)<=500)
               flag=0;%满足要求
            end         
        end
      %part2：判断节点覆盖范围内是否有其他节点
        local_neigber_count=0;%记录新增加节点的覆盖范围内节点个数
        X_min=X_Node(1,N_flag)-cov/2; %确定新增节点的覆盖区域
        X_max=X_Node(1,N_flag)+cov/2; 
        Y_min=Y_Node(1,N_flag)-cov/2;
        Y_max=Y_Node(1,N_flag)+cov/2;
        for j=1:N_flag
            if N_flag~=j           
                if ((X_min<=X_Node(1,j)) &&(X_Node(1,j)<=X_max))
                    if ((Y_min<=Y_Node(1,j)) && (Y_Node(1,j)<=Y_max))
                        A(N_flag,j)=1;%无向图
                        A(j,N_flag)=1;
                        local_neigber_count=local_neigber_count+1;%新增加节点的覆盖范围内节点个数                   
                    end
                end
            end
        end
        if local_neigber_count~=0
            tab=0;%满足条件
        end
    end
    %--更新新增节点邻居节点号记录矩阵--%
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
%% --Part 2--每个时步网络增加节点依概率与老节点连边m条（在节点覆盖范围内连边）------ %%
 m=12;%增加边数（固定一个值）      可变参数
 p=0.3;% 除去第一条边外，剩余的边以概率p三角连接      可变参数
 w0=1;%加入边的初始权值
 deta=1;%每次被连接，节点点强增加的数值
 alf=0.5;%   能量与信道的可调参数
mt=(m-1)*p;%三角连接的边的个数
mp=(m-mt);%择优连接的边的个数
 Pt=1;%(W)发射功率30dBm    
 M=1000;%主用户以及信道个数
bandwidth_sum=40*10^7;%系统总带宽10MHz
W0=bandwidth_sum/M;%每个子载波带宽
 WGN=10^(-17.4)*W0*10^-3;%噪声功率（W）  
An=zeros(1,1);%记录节点间可用信道集合
 if local_neigber_count~=0 %邻居节点数不为0
        if (1<=local_neigber_count && local_neigber_count<=m)%边全连
            for i=1:local_neigber_count
                out=Neigb_Record(N_flag,i);
                Edge(N_flag,out)=1; %连接边矩阵
                Edge(out,N_flag)=1; 
            end
        else
            %记录选出来连边的节点号
            selout_node=zeros(1,m);
            %%========step4:求新增节点覆盖范围内节点与新增节点间距离以及潜在链路容量
            Dis_TwoNode=zeros(1,local_neigber_count);%定义新增节点到覆盖范围内节点的距离矩阵
            PL_TwoNode=zeros(1,local_neigber_count);%定义新增节点到覆盖范围内节点间信道衰落矩阵！！！
            SNR_TwoNode=zeros(1,local_neigber_count);%定义新增节点到覆盖范围内节点间链路信噪比矩阵！！！
            Cap_TwoNode=zeros(1,local_neigber_count);%定义新增节点到覆盖范围内节点间链路容量矩阵！！！！
            for i=1:local_neigber_count
                Dis_TwoNode(1,i)=sqrt((X_Node(1,N_flag)-X_Node(1,Neigb_Record(N_flag,i)))^2+(Y_Node(1,N_flag)-Y_Node(1,Neigb_Record(N_flag,i)))^2);
                PL_TwoNode(1,i)=30.6+36.7*log10(Dis_TwoNode(1,i));%大尺度衰落与距离有关，小尺度衰落服从瑞利分布
                SNR_TwoNode(1,i)=Pt*PL_TwoNode(1,i)./WGN;
                for j=1:M
                   if H_state(1,j)~=1
                     Cap_TwoNode(1,i)=Cap_TwoNode(1,i)+W0*log2(1+SNR_TwoNode(1,i)); 
                   else
                        Cap_TwoNode(1,i)=Cap_TwoNode(1,i)+0;
                   end
                end               
            end
           %%=========step5:求增加边的概率
             Edge_P_Array=zeros(1,local_neigber_count); %节点被择优选择连边的概率
             Edge_P_mid=zeros(1,local_neigber_count);%边增加选择概率中间参数，关于能量和信道质量
             Neigb_Choose=zeros(1,local_neigber_count);%新增节点邻居节点
              Edge_P_mid_sum=0;%边增加选择概率中间参数求和
              %求中间概率参数
            for i=1:local_neigber_count
               Edge_P_mid(1,i)=(Power0(1,Neigb_Record(N_flag,i)).^alf)*(Cap_TwoNode(1,i).^(1-alf));
                Edge_P_mid_sum=Edge_P_mid_sum+(Edge_P_mid(1,i)*S(1,Neigb_Record(N_flag,i)));
            end
             %求最终概率
            for i=1:local_neigber_count
                Edge_P_Array(1,i)=(Edge_P_mid(1,i)*S(1,Neigb_Record(N_flag,i)))/Edge_P_mid_sum;
                Neigb_Choose(1,i)=Neigb_Record(N_flag,i);
            end
             %%===========step6:依概率增加边
             %增加择优选择的边
              mp_num=1;%初始值
            while mp_num<=mp %增加边直到增加到mp条为止
                out=randsrc(1,1,[Neigb_Choose(1,:);Edge_P_Array(1,:)]);%
                %判断选出的节点是否重复
                if length(find(selout_node==out))==0 %无重复的
                    selout_node(1,mp_num)=out;%记录
                    Edge(N_flag,out)=1; %连接边矩阵
                    Edge(out,N_flag)=1;
                    mp_num=mp_num+1;%边数加一
                end
            end
            %增加三角连接的边
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
                            mt_possible_num_p_sum= mt_possible_num_p_sum+W(out,N_Record(out,j))/S(1,out);  %为了归一化，求概率值的综合
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
                 
            if  mt_possible_num<=mt;%若所选择的节点的连接节点个数小于要三角连接的节点个数，就全连接
                  for i=1: length(find(selout_node2))
                     out=selout_node2(1,i);                      
                       Edge(N_flag,out)=1; 
                       Edge(out,N_flag)=1; 
                  end
            else  
              mt_num=1;
                   while mt_num<=mt %增加边直到增加到mt条为止
                   out=randsrc(1,1,[selout_node2(1,:); mt_possible_num_p(1,:)]);
                %判断选出的节点是否重复
                  if length(find(selout_node==out))==0 %无重复的
                    selout_node(1,mp_num)=out;%记录
                    Edge(N_flag,out)=1; %连接边矩阵
                    Edge(out,N_flag)=1;
                    mp_num=mp_num+1;%边数加一
                    mt_num=mt_num+1;
                     end
                   end
           end
        end
 else
       Edge(N_flag,:)=0; %孤立节点加入@@@@@@@@
        Edge(:,N_flag)=0; 
 end


  %% --增加1个节点和m条边过程后，更新网络节点强度与边权--- %%
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
   S(1,i)=sum(W(i,:));%求各节点节点点强
end
%% --更新邻接节点号--- %%
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
 %% --Part 3--每个时步以反择优概率删除网络删除连边n条（在全网范围选择）------ %%
 n=4;%设定删除链路数目
    %记录选出来删除链路一端节点号node1&node2
    selout_node1=zeros(1,n);
    selout_node2=zeros(1,n);
 %%====step1:求选择删除链路第1端的概率
    Node1_P_Array=zeros(1,1);%动态定义选择删除链路第1端概率
    Node1_P_mid=zeros(1,1);%选择删除链路第1端概率中间值
    Neigb1_Choose=zeros(1,1);%第1端节点连边节点
    Node1_P_mid_sum=0;%选择删除链路第1端概率中间值求和
     %求候选节点集合
    k=0;
    for i=1:N_flag
        if Netnode_record(1,i)~=0
            k=k+1;
            Neigb1_Choose(1,k)=Netnode_record(1,i);
        end
    end
    %--求Node1中间概率
    for i=1:k
        Node1_P_mid(1,i)=(Power0(1,Neigb1_Choose(1,i))*S(1,Neigb1_Choose(1,i)))^-1;
        Node1_P_mid_sum=Node1_P_mid_sum+Node1_P_mid(1,i);
    end
 %---求Node1最终概率
    for i=1:k
        Node1_P_Array(1,i)=Node1_P_mid(1,i)./Node1_P_mid_sum;
    end
  %%===step2:依概率删除边
   n_num=1;%初始值    
   while n_num<=n %删除边直到n条为止
     out1=randsrc(1,1,[Neigb1_Choose(1,:); Node1_P_Array(1,:)]);%**ONE:依概率选择链路的一端节点node1
    node_linked_count=sum(Edge(out1,:));%记录选出节点1的连边节点个数
     if node_linked_count~=0
        selout_node1(1,n_num)=out1;%记录
         %求连边节点矩阵
                linked_array=zeros(1,node_linked_count);
                t=0;
                for i=1:N_flag
                    if Edge(out1,i)~=0
                        t=t+1;%计数
                        linked_array(1,t)=i;
                    end
                end
      %%========step3:求node1节点连边节点与node1节点间距离以及潜在链路容量
            Dis_TwoNode1=zeros(1,node_linked_count);%定义与out1连接节点的距离
            PL_TwoNode1=zeros(1,node_linked_count);%定义与out1节点间信道衰落矩阵！！！
            SNR_TwoNode1=zeros(1,node_linked_count);%定义与out1节点间链路信噪比矩阵！！！
            Cap_TwoNode1=zeros(1,node_linked_count);%定义新增节点到覆盖范围内节点间链路容量矩阵！！！！
            for i=1:node_linked_count
               Dis_TwoNode1(1,i)=sqrt((X_Node(1,out1)-X_Node(1,linked_array(1,i)))^2+(Y_Node(1,out1)-Y_Node(1,linked_array(1,i)))^2);
               PL_TwoNode1(1,i)=30.6+36.7*log10(Dis_TwoNode1(1,i));%大尺度衰落与距离有关，小尺度衰落服从瑞利分布
               SNR_TwoNode1(1,i)=Pt*PL_TwoNode1(1,i)./WGN;
                for j=1:M
                   if H_state(1,j)==0
                      Cap_TwoNode1(1,i)=Cap_TwoNode1(1,i)+W0*log2(1+SNR_TwoNode1(1,i));
                   else
                        Cap_TwoNode1(1,i)=Cap_TwoNode1(1,i)+0;
                   end
                end               
            end
             %%=========step4:求删除边的概率
                Link_P_Array=zeros(1,node_linked_count);%动态定义删除边选择概率
                Link_P_mid=zeros(1,node_linked_count);%边删除选择概率中间值
                %Link_Choose=zeros(1,node_linked_count);%新增节点邻居节点
                Link_P_mid_sum=0;%边增加选择概率中间值求和
           %求中间概率
                for i=1:node_linked_count
                    Link_P_mid(1,i)=(Power0(1,linked_array(1,i)).^alf)*(Cap_TwoNode1(1,i).^(1-alf));
                    Link_P_mid_sum=Link_P_mid_sum+(Link_P_mid(1,i)*S(1,linked_array(1,i)))^-1;
                end
        %求最终概率
                for i=1:node_linked_count
                    Link_P_Array(1,i)=((Link_P_mid(1,i)*S(1,linked_array(1,i)))^-1)/Link_P_mid_sum;
                end
        %%===========step5:依概率删除一条（node2）
                out2=randsrc(1,1,[linked_array(1,:);Link_P_Array(1,:)]);%依概率选择删除边的另一端节点node2
                selout_node2(1,n_num)=out2;%记录
                Edge(out1,out2)=0; %连接边矩阵删除边清零
                Edge(out2,out1)=0;
                W(out1,out2)=0;
                W(out2,out1)=0;
 %%%%%%%%%%%%%%%%%%%%%%%%%可能需要调节的参数%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
 n_num=n_num+1;%删除边数加一
     end 
   end
    %% -----------更新点强----------- %%
  for i=1:N_flag
   S(1,i)=sum(W(i,:));%求各节点节点点强
  end
%-------------更新邻接节点号----------------%

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
  %% -----------------计算所剩边数------------------ %%
sum_edge1=0;
for i=1:N_flag
   sum_edge1=sum_edge1+sum(Edge(i,:))/2;
end
   %% --Part4:每一时步，节点间的信道以概率AP被占用，节点间可以进行信道切换的概率为CP -- %%
    H_state_before=H_state;           %将现在的信道状态赋值给before
    for j=1:H
        time_initial(1,j)=rem(time_initial(1,j)+1,cycle(1,j));              %更新信道状态
     if  time_initial(1,j)<=time_on(1,j)
        H_state(1,j)=1;
     end
      if  time_initial(1,j)>time_on(1,j)
        H_state(1,j)=0;
      end
    end
 num_migration(1,ttt)= numel( find( H_state_before-H_state==-1));     %记录信道状态由0变为1的数目
 AP=num_migration(1,ttt)/sum_edge1 ;                   %主用户的活跃度
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
   % -----------更新点强----------- %%
  for i=1:N_flag
   S(1,i)=sum(W(i,:));%求各节点节点点强
  end
  % -----------更新节点度----------- %%
  for i=1:N_flag
   K(1,i)=sum(Edge(i,:));%求各节点节点点强
  end
  %-------------更新邻接节点号----------------%

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
  %% -----------Part5: 节点度为0的孤立节点删除----------- %%
    k=0;
    for i=1:N_flag
       if  S(1,i)==0         
           A(i,:)=0;%网络邻接矩阵
           A(:,i)=0;           
          % Power0(1,i)=0;%N个节点建立拓扑后的初始能量
           X_Node(1,i)=0;%节点横坐标清除
           Y_Node(1,i)=0;%节点纵坐标清除
           Netnode_record(1,i)=0;%记录网络中存在的节点标号清除  
           k=k+1;%全网节点数更新
       end
    end
    node_num=N_flag-k;
    
end
 %% -----------------排除孤立节点，建立边的连接矩阵------------------ %%
  node_final=zeros(1,node_num);
  node_tmp=0;
  for i=1:N_flag
      if sum(Edge(i,:))~=0;
          node_tmp=node_tmp+1;
          node_final(1,node_tmp)=i;               %%挑选出不孤立的节点
      end
  end
   Edge_f=zeros(1,1);
  for i=1:(node_num-1) 
     for j=i+1:node_num
        if Edge( node_final(1,i), node_final(1,j))~=0              %%  建立不孤立节点之间的边矩阵
           Edge_f(i,j)=1;
           Edge_f(j,i)=1; 
        end
     end
  end  
%% -----------------计算所剩边数------------------ %%
sum_edge=0;
for i=1:N_flag
   sum_edge=sum_edge+sum(Edge(i,:))/2;
end
%% -----------------计算总点强------------------ %%
sum_s=0;
for i=1:N_flag
   sum_s=sum_s+sum(W(i,:));
end
%% - ----------------模拟攻击随机删除一个节点------------------ %%
out_s=randsrc(1,1,1:node_num);
Edge_s=Edge_f;
Edge_s(:,out_s)=[];
Edge_s(out_s,:)=[];
%% - ----------------模拟刻意攻击删除一个节点------------------ %%
[zd,wz]=max(S);
out_k=find(node_final==wz);
Edge_k=Edge_f;
Edge_k(:,out_k)=[];
Edge_k(out_k,:)=[];
%% -----------------绘制网络拓扑图------------------ %%
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
             hold on;      %% 画出BA无标度网络图
         end
     end
end
 axis equal;
grid on
hold off 
%  %% -----------数据存储操作------------ %%
%  
% eval(['save S',num2str(T)]); %%%保存每次循环的点强分布
% saveas(gca,num2str(T),'jpg');   %%%保存每次循环产生的图片
% close(gcf);
% node_n(1,aaa)= node_num;    %%%%  记录每次剩余的节点数
%  edge_n(1,aaa)=sum_edge;
% eval(['save  Edge_f',num2str(T)]); %%%保存每次循环的连接矩阵
%  end
% savefile1='node_n.mat'; %保存每次剩余的节点数
% save(savefile1,'node_n');
% savefile1='edge_n.mat'; %保存每次剩余的总边数
% save(savefile1,'edge_n');



