clear
clc
tic
%% Parameter Input

 Pipe number  From  To  Length (m）  Inner diameter (m)  Heat transfer coefficient (kW/m·K) 
 Pipe_info_input=[1	1	2	15	0.606	1.052132423
 2	2	3	970	0.309	0.604567196
 3	3	4	70	0.309	0.604567196
 4	3	5	1550	0.309	0.604567196
 5	2	6	204	0.606	1.052132423
 6	6	7	153	0.357	0.78493095
 7	7	8	270	0.1	0.358623086
 8	7	9	300	0.357	0.78493095
 9	9	10	200	0.1	0.358623086
 10	9	11	300	0.357	0.78493095
 11	11	12	200	0.1	0.358623086
 12	11	13	450	0.357	0.78493095
 13	13	14	352	0.309	0.604567196
 14	13	15	1634	0.357	0.78493095
 15	6	16	362	0.357	0.78493095
 16	16	17	830	0.207	0.527514587
 17	16	18	723	0.357	0.78493095
 18	18	19	218	0.259	0.580986402
 19	18	20	2250	0.309	0.78493095
 ];
 
%Node_number Source_injected_flow_rate (kg/s)  Node_load (kg/s)  Node_pressure (Pa)  Source_temperature (℃)
 Node_info_input=[  1	-1000	0	981000	200
 2	0	0	0	0
 3	0	0	0	0
 4	0	4.166666667	0	0
5	0	4	0	0
 6	0	0	0	0
 7	0	0	0	0
 8	0	1.111111111	0	0
 9	0	0	0	0
 10	0	0.555555556	0	0
 11	0	0	0	0
 12	0	0.833333333	0	0
 13	0	0	0	0
 14	0	8.444444444	0	0
 15	0	3.6	0	0
 16	0	0	0	0
 17	0	2.777777778	0	0
 18	0	0	0	0
 19	0	6.666666667	0	0
 20	0	5	0	0
  ];         
 Ta=10; %Ambient temperature (℃)
 Ck=1.2;% Local heat dissipation loss coefficient
 eta=0.2;%Local fiction coefficient
 error_P0=50;% overall allowable error of pressure (Pa)
 error_T0=0.05;% overall allowable error of temperature (K)


%-----------------------------------------------------------------------------
Pipe_info=Pipe_info_input;
Node_info=Node_info_input;
n_node=size(Node_info,1);%Node number
n_pipe=size(Pipe_info,1);%pipe number
          
if sum(Node_info(:,2:3)==-1000)+sum(Node_info(:,4)==0)>n_node
    error('已知信息不足，无法计算');
else if  sum(Node_info(:,2:3)==-1000)+sum(Node_info(:,4)==0)<n_node
        error('已知信息过多，无法计算');
    end
end

%节点重新排序
Ps_orig=find(Node_info(:,4)~=0);
n_Ps=length(Ps_orig);
if sum(Node_info(:,5)~=0)~=n_Ps
    error('已知压强与温度数量不匹配');
end
for i=1:n_Ps
    if i~=Ps_orig(i)
    temp=Node_info(i,:);
    Node_info(i,:)=Node_info(Ps_orig(i),:);
    Node_info(Ps_orig(i),:)=temp;
    Pipe_info(find(Pipe_info_input(:,2)==Ps_orig(i)),2)=i;
    Pipe_info(find(Pipe_info_input(:,2)==i),2)=Ps_orig(i);
    Pipe_info(find(Pipe_info_input(:,3)==Ps_orig(i)),3)=i;
    Pipe_info(find(Pipe_info_input(:,3)==i),3)=Ps_orig(i);
    end
end
%生成关联矩阵
A=zeros(n_node,n_pipe);%节点管段关联矩阵
for i=1:n_pipe
    A(Pipe_info(i,2),i)=-1;
    A(Pipe_info(i,3),i)=1;
end
pipeline_info=zeros(n_pipe,2);
Ai=zeros(n_node,n_pipe);Ao=zeros(n_node,n_pipe);
for i = 1: n_node
    for j = 1: n_pipe
        if A(i,j)==1
            Ai(i,j)=1; pipeline_info(j,2)=i;
        end
        if A(i,j)==-1
            Ao(i,j)=1; pipeline_info(j,1)=i;
        end
    end
end
Qn0=Node_info(:,2)-Node_info(:,3);
Qn0(find(Node_info(:,2)==-1000))=-1000;
Qn0(find(Node_info(:,3)==-1000))=-1000;
P0=Node_info(:,4);%节点压强Pa,除了平衡节点外都为与平衡节点的差值
Ts=Node_info(:,5)+273.15;%摄氏度转化成K
l=Pipe_info(:,4);%管段当量长度m
% error_P=error_P0/sum(l);
% error_T=error_T0/sum(l);
d=Pipe_info(:,5);%管段内径m
%y=fsolve(@(x)(1/x+2*log10(0.00006757+2.51/122804/x)),[0.5])
lambda=0.11*(0.0002./d).^0.25;%管段内工质摩擦阻力系数+68/112974
K=Ck*1e-3*Pipe_info(:,6);%管段传热系数，换算至kW/m·K
Ta=Ta+273.15;
% h=xlsread('h.xlsx','A2:C21');
% h(:,1)=h(:,1)*1000;
% n_h=size(h,1);
%T=xlsread('h.xlsx','D2:D21');

%%生成求解顺序
Pipe_Sequence=zeros(n_pipe,1);
i=0;
for j=1:n_pipe
    if Ai(pipeline_info(j,1),:)==0
        i=i+1;
        Pipe_Sequence(i)=j;
    end
end
while i<n_pipe
    for j=1:n_pipe
        if ~sum(Pipe_Sequence==j)
            temp=find(pipeline_info(:,2)==pipeline_info(j,1));
            mark=1;
            for j2=1:length(temp)
                if ~sum(Pipe_Sequence==temp(j2))
                    mark=0;
                end
            end
            if mark==1;
                i=i+1;
                Pipe_Sequence(i)=j;
            end        
        end
    end
end



%%水力计算
P_init=max(P0);
P=P0;
P(find(P==0))=P_init;
Z=ones(n_pipe,1)*mean(GetZ(P0(find(P0~=0))*1e-6,Ts(find(Ts~=273.15))));
Tin=ones(n_pipe,1)*mean(Ts(find(Ts~=273.15)));
Tin(find(Ts~=273.15))=Ts(find(Ts~=273.15));
Tout=Tin;

R=461.526;%J/(kg·K)
%转换成Pa,K
c1=0.285764132334731e-6;c2=1.979104677597431;
c3=-1.893766236738926e-4;c4=1.945912391397789e+03;
M=ones(n_pipe,1)*10;%管道流量初始值取10kg/s
Qn=Qn0;
Qn(find(Qn0==-1000))=-sum(Qn0(find(Qn0~=-1000)))/sum(Qn0==-1000);
Qn_in=Node_info(:,2);
deltaP=1000;
CountIteration_total=0;%总迭代次数
CountIteration1=0;%大循环次数

s=d.^5*pi^2/8./lambda./Z/R./(Tin+Tout)./l/(1+eta);
Mc=zeros(n_pipe,1);
% n_AddedNodes=[];

while CountIteration_total==0 || sum(n_AddedNodes)>0
CountIteration1=CountIteration1+1;
CountIteration2=0;%单次循环迭代次数
while CountIteration2==0 || (CountIteration2<10000000  && ((max(abs(P-P_last))>0.5) || (max(abs(M-M_last))>0.005) || (max(abs(T-T_last))>0.005)))
    CountIteration2=CountIteration2+1;
    P_last=P;
    G=diag(s./M);
    temp1=A*G*A';temp2=Qn-(Ai+Ao)*Mc/2-A*G*A'*P.^2;
    deltaP=inv(temp1(n_Ps+1:n_node,n_Ps+1:n_node))*temp2(n_Ps+1:n_node);
    P(n_Ps+1:n_node)=real((P(n_Ps+1:n_node).^2+deltaP).^0.5);
%    P(find(P0~=0))=P0(find(P0~=0));
    if CountIteration2>1 
        Z=0.5*(GetZ(Ai'*P*1e-6,Tout)+GetZ(Ao'*P*1e-6,Tin));
    end
    Z(find(Z<0.9))=0.9;
    s=d.^5*pi^2/8./lambda./Z/R./(Tin+Tout)./l/(1+eta);
    M=-s.*(A'*(P.^2));
    M=real(sign(M).*sqrt(abs(M)));
    for j=1:n_pipe
        if sign(M)==-1
            M=-M;
            temp=pipeline_info(j,1);
            pipeline_info(j,1)= pipeline_info(j,2);
            pipeline_info(j,2)= pipeline_info(j,1);
            A(pipeline_info(j,1),j)=-1;A(pipeline_info(j,2),j)=1;
        end
    end
    Qn=-Ai*(M-Mc/2)+Ao*(M+Mc/2);
    Qn(find(Qn0~=-1000))=Qn0(find(Qn0~=-1000));
    Qn_in(find(Qn_in==-1000))=Qn(find(Qn_in==-1000))-Node_info(find(Qn_in==-1000),3);
    M_diag=diag(M);
    if CountIteration2==1 
        T_last=zeros(n_pipe*2+n_node,1);
    else
        T_last=T;
    end
    Jacob_T=[diag(K)*diag(l)/2-diag(c1.*M_diag*Ao'*P)-c2*M_diag        diag(K)*diag(l)/2+diag(c1*M_diag*Ai'*P)+c2*M_diag               zeros(n_pipe,n_node)
                                 zeros(n_node,n_pipe)                            -Ai*M_diag                                            diag(Ai*M+Qn_in)
                                 eye(n_pipe)                                     zeros(n_pipe)                                            -Ao' ];

    T=inv(Jacob_T)*[Ta*K.*l-c3*M_diag*A'*P ;   Qn_in.*Ts     ; zeros(n_pipe,1)       ];
    Tin=T(1:n_pipe);
    Tout=T(n_pipe+1:2*n_pipe);
    Tn=T(n_pipe*2+1:n_pipe*2+n_node);

   
%     若不考虑冷凝水，则注释掉下面部分
    for j=1:n_pipe
        Tin(Pipe_Sequence(j))=(Ai(pipeline_info(Pipe_Sequence(j),1),:)*diag(M-Mc/2)*Tout+Qn_in(pipeline_info(Pipe_Sequence(j),1))*Ts(pipeline_info(Pipe_Sequence(j),1)))/(Ai(pipeline_info(Pipe_Sequence(j),1),:)*(M-Mc/2)+Qn_in(pipeline_info(Pipe_Sequence(j),1)));
        ts=202.7651*(P(pipeline_info(Pipe_Sequence(j),2))*1e-6).^(0.2148)-22.907+273.15;
        if Tout(Pipe_Sequence(j))<ts
            Tout(Pipe_Sequence(j))=ts; 
            P_av=(P(pipeline_info(Pipe_Sequence(j),1))+P(pipeline_info(Pipe_Sequence(j),2)))/2*1e-6;
            hc=3704.12*P_av^(0.6)/(1-2.22271*P_av^(0.2)+7.55636*P_av^(0.4)-1.476876*P_av^(0.6));
            Mc1=M(Pipe_Sequence(j))*((c1*P(pipeline_info(Pipe_Sequence(j),2))+c2)*ts+c3*P(pipeline_info(Pipe_Sequence(j),2))-(c1*P(pipeline_info(Pipe_Sequence(j),1))+c2)*Tin(Pipe_Sequence(j))-c3*P(pipeline_info(Pipe_Sequence(j),1)))+K(Pipe_Sequence(j))/2*(Tin(Pipe_Sequence(j))+ts-2*Ta)*l(Pipe_Sequence(j));
            Mc2=(0.5*((c1*P(pipeline_info(Pipe_Sequence(j),2))+c2)*ts+c3*P(pipeline_info(Pipe_Sequence(j),2))+(c1*P(pipeline_info(Pipe_Sequence(j),1))+c2)*Tin(Pipe_Sequence(j))+c3*P(pipeline_info(Pipe_Sequence(j),1)))+c4-hc);
            Mc(Pipe_Sequence(j))=Mc1/Mc2;
        else
            Mc(Pipe_Sequence(j))=0;
        end
    end
    Tn=inv(diag(Ai*M+Qn_in))*(Ai*M_diag*Tout+  Qn_in.*Ts);
    T=[Tin; Tout; Tn];
    if CountIteration2>1 
        Z=0.5*(GetZ(Ai'*P*1e-6,Tout)+GetZ(Ao'*P*1e-6,Tin));
    end
    Z(find(Z<0.9))=0.9;
    M_last=M;
    s=d.^5*pi^2/8./lambda./Z/R./(Tin+Tout)./l/(1+eta);
    M=-s.*(A'*(P.^2));
    M=real(sign(M).*sqrt(abs(M)));
    Qn=-Ai*(M-Mc/2)+Ao*(M+Mc/2);
    Qn(find(Qn0~=-1000))=Qn0(find(Qn0~=-1000));
end
CountIteration_total=CountIteration_total+CountIteration2;%P(2)*1E-3 %Tout-273.15
%(P(20)^2-P(18)^2)+8*lambda(19)*GetZ(0.5*(P(20)+P(18))*1e-6,0.5*(Tin(19)+Tout(19)))*(Tin(19)+Tout(19))*R*l(19)/pi/pi/(d(19)^5)*M(19)*M(19)
%(P(18)^2-P(16)^2)+8*lambda(17)*GetZ(0.5*(P(16)+P(18))*1e-6,0.5*(Tin(17)+Tout(17)))*(Tin(17)+Tout(17))*R*l(17)/pi/pi/(d(17)^5)*M(17)*M(17)
%步长校正
n_AddedNodes=zeros(n_pipe,1);
% error_P=error_P0/n_pipe;error_T=error_T0/n_pipe;
for j=1:n_pipe
    h=l(j);
    error_P=error_P0*(P(pipeline_info(j,1))-P(pipeline_info(j,2)))/(max(P)-min(P));
    error_T=error_T0*(Tin(j)-Tout(j))/(max(Tin)-min(Tout));
    if Mc(j)==0
        K1=inv(GetA3(P(pipeline_info(j,1)),Tin(j),M(j),d(j)))*GetB3(P(pipeline_info(j,1)),Tin(j),M(j),d(j),lambda(j)*(1+eta),K(j),Ta);
        K2=inv(GetA3(P(pipeline_info(j,1))+K1(2)*h/2,Tin(j)+K1(1)*h/2,M(j),d(j)))*GetB3(P(pipeline_info(j,1))+K1(2)*h/2,Tin(j)+K1(1)*h/2,M(j),d(j),lambda(j)*(1+eta),K(j),Ta);
        K3=inv(GetA3(P(pipeline_info(j,1))+K2(2)*h/2,Tin(j)+K2(1)*h/2,M(j),d(j)))*GetB3(P(pipeline_info(j,1))+K2(2)*h/2,Tin(j)+K2(1)*h/2,M(j),d(j),lambda(j)*(1+eta),K(j),Ta);
        K4=inv(GetA3(P(pipeline_info(j,1))+K3(2)*h,Tin(j)+K3(1)*h,M(j),d(j)))*GetB3(P(pipeline_info(j,1))+K3(2)*h,Tin(j)+K3(1)*h,M(j),d(j),lambda(j)*(1+eta),K(j),Ta);
        dK=(K1+2*K2+2*K3+K4)*h/6;
        T_RK4(j)=Tin(j)+dK(1);
        P_RK4(j)=P(pipeline_info(j,1))+dK(2); 
%         abs(P_RK4(j)-P(pipeline_info(j,2)))
%         abs(T_RK4(j)-Tout(j))
        if abs(P_RK4(j)-P(pipeline_info(j,2)))>error_P || abs(T_RK4(j)-Tout(j))>error_T
            n_AddedNodes(j)=min(3,max(floor((abs(P_RK4(j)-P(pipeline_info(j,2)))/error_P)^(1/3)),floor((abs(T_RK4(j)-Tout(j))/error_T)^(1/3))));
        end
    else
        K1=inv(GetA3(P(pipeline_info(j,1)),Tin(j),M(j),d(j)))*GetB3(P(pipeline_info(j,1)),Tin(j),M(j),d(j),lambda(j)*(1+eta),K(j),Ta);
        if 202.7651*((P(pipeline_info(j,1))+K1(2)*h/2)*1e-6).^(0.2148)-22.907+273.15>Tin(j)+K1(1)*h/2
            K1(2)=GetB4(P(pipeline_info(j,1)),M(j),d(j),lambda(j)*(1+eta));
            K1(1)=(202.7651*((P(pipeline_info(j,1))+K1(2)*h/2)*1e-6).^(0.2148)-22.907+273.15-Tin(j))/h*2;
        end
        K2=inv(GetA3(P(pipeline_info(j,1))+K1(2)*h/2,Tin(j)+K1(1)*h/2,M(j),d(j)))*GetB3(P(pipeline_info(j,1))+K1(2)*h/2,Tin(j)+K1(1)*h/2,M(j),d(j),lambda(j)*(1+eta),K(j),Ta);
        if 202.7651*((P(pipeline_info(j,1))+K2(2)*h/2)*1e-6).^(0.2148)-22.907+273.15>Tin(j)+K2(1)*h/2
            K2(2)=GetB4(P(pipeline_info(j,1))+K1(2)*h/2,M(j),d(j),lambda(j)*(1+eta));
            K2(1)=(202.7651*((P(pipeline_info(j,1))+K2(2)*h/2)*1e-6).^(0.2148)-22.907+273.15-Tin(j))/h*2;
        end
        K3=inv(GetA3(P(pipeline_info(j,1))+K2(2)*h/2,Tin(j)+K2(1)*h/2,M(j),d(j)))*GetB3(P(pipeline_info(j,1))+K2(2)*h/2,Tin(j)+K2(1)*h/2,M(j),d(j),lambda(j)*(1+eta),K(j),Ta);
        if 202.7651*((P(pipeline_info(j,1))+K3(2)*h)*1e-6).^(0.2148)-22.907+273.15>Tin(j)+K3(1)*h
            K3(2)=GetB4(P(pipeline_info(j,1))+K2(2)*h/2,M(j),d(j),lambda(j)*(1+eta));
            K3(1)=(202.7651*((P(pipeline_info(j,1))+K3(2)*h)*1e-6).^(0.2148)-22.907+273.15-Tin(j))/h;
        end
        K4=inv(GetA3(P(pipeline_info(j,1))+K3(2)*h,Tin(j)+K3(1)*h,M(j),d(j)))*GetB3(P(pipeline_info(j,1))+K3(2)*h,Tin(j)+K3(1)*h,M(j),d(j),lambda(j)*(1+eta),K(j),Ta);
        if 202.7651*((P(pipeline_info(j,1))+K4(2)*h)*1e-6).^(0.2148)-22.907+273.15>Tin(j)+K3(1)*h
            K4(2)=GetB4(P(pipeline_info(j,1))+K3(2)*h/2,M(j),d(j),lambda(j)*(1+eta));
            K4(1)=(202.7651*((P(pipeline_info(j,1))+K4(2)*h)*1e-6).^(0.2148)-22.907+273.15-Tin(j))/h;
        end
        dK=(K1+2*K2+2*K3+K4)*h/6;
        P_RK4(j)=P(pipeline_info(j,1))+dK(2); 
        T_RK4(j)=202.7651*(P_RK4(j)*1e-6).^(0.2148)-22.907+273.15;
        if abs(P_RK4(j)-P(pipeline_info(j,2)))>error_P 
            n_AddedNodes(j)=min(3,floor((abs(P_RK4(j)-P(pipeline_info(j,2)))/error_P)^(1/3)));
        end        
    end
end
if sum(n_AddedNodes)>0
for j=1:n_pipe
    Pipe_info(j,4)=Pipe_info(j,4)/(n_AddedNodes(j)+1);
    Pipe_Sequence(find(Pipe_Sequence==j)+n_AddedNodes(j)+1:end+n_AddedNodes(j))=Pipe_Sequence(find(Pipe_Sequence==j)+1:end);
    for jj=1:n_AddedNodes(j)
        Node_info(end+1,:)=[Node_info(end,1)+1 0 0 0 0];
        Pipe_info(end+1,:)=[Pipe_info(end,1)+1 Node_info(end,1)  Node_info(end,1)+1  Pipe_info(j,4:6)];
        M(Pipe_info(end,1),1)=M(j);
        Mc(Pipe_info(end,1),1)=0;
        P(Node_info(end,1),1)=P(pipeline_info(j,1))*(n_AddedNodes(j)+1-jj)/(n_AddedNodes(j)+1)+P(pipeline_info(j,2))*jj/(n_AddedNodes(j)+1);
        Tn(Node_info(end,1),1)=Tn(pipeline_info(j,1))*(n_AddedNodes(j)+1-jj)/(n_AddedNodes(j)+1)+Tn(pipeline_info(j,2))*jj/(n_AddedNodes(j)+1);
        Tin(Pipe_info(end,1),1)=Tn(Node_info(end,1));
        if jj>1
            Tout(Pipe_info(end,1)-1,1)=Tn(Node_info(end,1));
        end
        Pipe_Sequence(find(Pipe_Sequence==j)+jj)=Pipe_info(end,1);
    end
    if n_AddedNodes(j)>0
        Tout(Pipe_info(end,1),1)=Tout(j);
        Tout(j,1)=Tn(Node_info(end,1)-n_AddedNodes(j)+1);
        Pipe_info(end,3)=Pipe_info(j,3);
        Pipe_info(j,3)=Node_info(end,1)-n_AddedNodes(j)+1;
        Mc(Pipe_info(end,1),1)=0;
    end
end
Qn0=[Qn0; zeros(sum(n_AddedNodes),1)];
Qn=[Qn; zeros(sum(n_AddedNodes),1)];
Qn_in=Node_info(:,2);
Ts=Node_info(:,5)+273.15;%摄氏度转化成K
l=Pipe_info(:,4);%管段当量长度m
d=Pipe_info(:,5);%管段内径m
lambda=0.11*(0.0002./d).^0.25;%管段内工质摩擦阻力系数+68/112974
K=Ck*1e-3*Pipe_info(:,6);%管段传热系数，换算至kW/m·K
n_node=n_node+sum(n_AddedNodes);
n_pipe=n_pipe+sum(n_AddedNodes);
%生成关联矩阵
A=zeros(n_node,n_pipe);%节点管段关联矩阵
for i=1:n_pipe
    A(Pipe_info(i,2),i)=-1;
    A(Pipe_info(i,3),i)=1;
end
pipeline_info=zeros(n_pipe,2);
Ai=zeros(n_node,n_pipe);Ao=zeros(n_node,n_pipe);
for i = 1: n_node
    for j = 1: n_pipe
        if A(i,j)==1
            Ai(i,j)=1; pipeline_info(j,2)=i;
        end
        if A(i,j)==-1
            Ao(i,j)=1; pipeline_info(j,1)=i;
        end
    end
end
Z=0.5*(GetZ(Ai'*P*1e-6,Tout)+GetZ(Ao'*P*1e-6,Tin));
Z(find(Z<0.9))=0.9;
s=d.^5*pi^2/8./lambda./Z/R./(Tin+Tout)./l/(1+eta);
end

end
    


%%蒸汽网计算结果显示
P=P*1e-3;%单位换算:压强 kPa
disp('蒸汽网计算结果：');
disp('蒸汽网迭代次数N=');
disp(CountIteration_total);
disp('各节点汽压 各节点注入流量');
disp('   (kPa)       (kg/s)');
disp([P Qn ]);
disp(' 各支路平均    冷凝水量 首端温度  末端温度');
disp(' 流量（kg/s）  （kg/s ）  (℃)      (℃)');
disp([M Mc Tin-273.15    Tout-273.15] );

toc
P(2)
min(Tout)-273.15
sum(Mc)

%单管道结果展示
X(1)=0;
for i=1:n_pipe
    X(i+1)=X(i)+l(Pipe_Sequence(i));
end
figure;
plot(X,[P(1);P(Pipe_info(Pipe_Sequence,3))]);
function A=GetA3(P,T,M,D) %P:Mpa,T:K
    A=zeros(2,2);
    R=461.526;%J/(kg·K)
    %转换成Pa,K
    c1=0.285764132334731e-6;c2=1.979104677597431;
    c3=-1.893766236738926e-4;c4=1.945912391397789e+03;
    Z=max(GetZ(P*1e-6,T),0.9);
    v=4*M*Z*R*T/pi/P/D/D;
    A(1,1)=0;
    A(1,2)=(Z*R*T)/P;
    A(2,1)=c1*P+c2;
    A(2,2)=c1*T+c3;
end
function B=GetB3(P,T,M,D,lambda,K,Ta) %P:Mpa,T:K %若计及局部阻力系数eta,则在调用时输入lambda*(1+eta）即可
    B=zeros(2,1);R=461.526;%J/(kg·K)
    Z=max(GetZ(P*1e-6,T),0.9);
    v=4*M*Z*R*T/pi/P/D/D;
    B=[-lambda*v*v/2/D;-K*(T-Ta)/M];
end

function B=GetB4(P,M,D,lambda) %P:Mpa,T:K %若计及局部阻力系数eta,则在调用时输入lambda*(1+eta）即可
    B=0;R=461.526;%J/(kg·K)
    T=202.7651*(P*1e-6).^(0.2148)-22.907+273.15;
    Z=max(GetZ(P*1e-6,T),0.9);
%     v=4*M*Z*R*T/pi/P/D/D;
    B=-8*lambda*M*M*Z*R*T/pi/pi/P/D^5;
end

                                   
