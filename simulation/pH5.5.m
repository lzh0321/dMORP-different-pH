%% Pgluc
vmaxg=22.32;% maximum glucose specific consumption rate (mmol/gDCW/h)
Kmg=25;% saturation constant (mmol/L)
Sg=zeros(12,1);
Sg(1)=110.596;% initial glucose concentration for the first step (mmol/L)
Slg=zeros(12,1);
Slg(1)=0;% initial lactate concentration for the first step (mmol/L)
dcwg=zeros(12,1);
dcwg(1)=0.159;% initial biomass for the first step (g/L)
v1=zeros(12,1);
v1(1)=vmaxg*Sg(1)/(Kmg+Sg(1));% glucose specific consumption rate for the first step (mmol/gDCW/h)
miug=zeros(12,1);
vlg=zeros(12,1);
FBAg=zeros(2460,12);
for j=1:12
    model_glu(j) = ecModel;
    model_glu(j) = changeRxnBounds(model_glu(j),'EX0003',0.01*v1(j),'b');%CO2
    model_glu(j) = changeRxnBounds(model_glu(j),'EX0004',-0.0044*v1(j),'b');%O2
    model_glu(j) = changeRxnBounds(model_glu(j),'EX0068',0.0562*v1(j),'b');%acetate
    model_glu(j) = changeRxnBounds(model_glu(j),'EX0069',0.0077*v1(j),'b');%citrate
    model_glu(j) = changeRxnBounds(model_glu(j),'EX0070',0.0015*v1(j),'b');%pyruvate
    model_glu(j) = changeRxnBounds(model_glu(j),'EX0027',0.0188*v1(j),'b');%acetoin
    model_glu(j) = changeRxnBounds(model_glu(j),'EX0002',1.86*v1(j),'b');%lactate
    % AA
    model_glu(j) = changeRxnBounds(model_glu(j),'EX0014',0.0141*v1(j),'b');%ala
    model_glu(j) = changeRxnBounds(model_glu(j),'EX0015',-0.1,'l');%arg
    model_glu(j) = changeRxnBounds(model_glu(j),'EX0016',-0.1,'l');%cys
    model_glu(j) = changeRxnBounds(model_glu(j),'EX0017',-0.0216*v1(j),'l');%glu
    model_glu(j) = changeRxnBounds(model_glu(j),'EX0018',-0.0013*v1(j),'l');%his
    model_glu(j) = changeRxnBounds(model_glu(j),'EX0019',-0.0153*v1(j),'l');%ile
    model_glu(j) = changeRxnBounds(model_glu(j),'EX0020',-0.0092*v1(j),'l');%leu
    model_glu(j) = changeRxnBounds(model_glu(j),'EX0021',-0.0038*v1(j),'l');%met
    model_glu(j) = changeRxnBounds(model_glu(j),'EX0022',-0.0326*v1(j),'l');%thr
    model_glu(j) = changeRxnBounds(model_glu(j),'EX0023',-0.0059*v1(j),'l');%val
    model_glu(j) = changeRxnBounds(model_glu(j),'EX0024',-0.0023*v1(j),'l');%tyr
    model_glu(j) = changeRxnBounds(model_glu(j),'T0036',0,'b');
    model_glu(j) = changeRxnBounds(model_glu(j),'T0037',0,'b');

    model_glu(j) = changeRxnBounds(model_glu(j),'EX0001',-v1(j),'b');
    model_glu(j) = changeObjective(model_glu(j),'EXBiomass');
    FBAsolution = optimizeCbModel(model_glu(j),'max','one');
    FBAg(:,j)=FBAsolution.x;
    vlg(j)=FBAsolution.x(3);
    miug(j)=FBAsolution.x(86);
    dcwg(j+1)=exp(miug(j)*1+log(dcwg(j)));% initial biomass fo the next step (g/L)
    Sg(j+1)=Sg(j)-v1(j)/miug(j)*(dcwg(j+1)-dcwg(j));% initial glucose concentration for the next step (mmol/L)
    Slg(j+1)=Slg(j)+vlg(j)/miug(j)*(dcwg(j+1)-dcwg(j));% initial lactate concentration for the next step (mmol/L)
    v1(j+1)=vmaxg*Sg(j+1)/(Kmg+Sg(j+1));% glucose specific consumption rate for the next step (mmol/gDCW/h)
end
%% 二次模拟
vmaxg=22.32;% maximum glucose specific consumption rate (mmol/gDCW/h)
Kmg=25;% saturation constant (mmol/L)
Sg=zeros(12,1);
Sg(1)=110.596;% initial glucose concentration for the first step (mmol/L)
Slg=zeros(12,1);
Slg(1)=0;% initial lactate concentration for the first step (mmol/L)
dcwg=zeros(12,1);
dcwg(1)=0.159;% initial biomass for the first step (g/L)
v1=zeros(12,1);
v1(1)=vmaxg*Sg(1)/(Kmg+Sg(1));% glucose specific consumption rate for the first step (mmol/gDCW/h)
vace=zeros(12,1);
vcit=zeros(12,1);
vpyr=zeros(12,1);
vacetoin=zeros(12,1);
Sace(1)=25.806;
Scit(1)=3.547;
Spyr(1)=0.0029;
Sacetoin(1)=0;
pro=zeros(12,1);
vlg=zeros(12,1);
FBAg=zeros(2460,12);
for j=1:12
    model_glu(j) = ecModel;
    model_glu(j) = changeRxnBounds(model_glu(j),'EX0003',0.01*v1(j),'b');%CO2
    model_glu(j) = changeRxnBounds(model_glu(j),'EX0004',-0.0044*v1(j),'b');%O2
    model_glu(j) = changeRxnBounds(model_glu(j),'EX0068',0.0562*v1(j),'b');%acetate
    model_glu(j) = changeRxnBounds(model_glu(j),'EX0069',0.0077*v1(j),'b');%citrate
    model_glu(j) = changeRxnBounds(model_glu(j),'EX0070',0.0015*v1(j),'b');%pyruvate
    model_glu(j) = changeRxnBounds(model_glu(j),'EX0027',0.0188*v1(j),'b');%acetoin
    model_glu(j) = changeRxnBounds(model_glu(j),'EX0002',1.86*v1(j),'b');
    % AA
    model_glu(j) = changeRxnBounds(model_glu(j),'EX0014',0.0141*v1(j),'b');%ala
    model_glu(j) = changeRxnBounds(model_glu(j),'EX0015',-0.1,'l');%arg
    model_glu(j) = changeRxnBounds(model_glu(j),'EX0016',-0.1,'l');%cys
    model_glu(j) = changeRxnBounds(model_glu(j),'EX0017',-0.0216*v1(j),'l');%glu
    model_glu(j) = changeRxnBounds(model_glu(j),'EX0018',-0.0013*v1(j),'l');%his
    model_glu(j) = changeRxnBounds(model_glu(j),'EX0019',-0.0153*v1(j),'l');%ile
    model_glu(j) = changeRxnBounds(model_glu(j),'EX0020',-0.0092*v1(j),'l');%leu
    model_glu(j) = changeRxnBounds(model_glu(j),'EX0021',-0.0038*v1(j),'l');%met
    model_glu(j) = changeRxnBounds(model_glu(j),'EX0022',-0.0326*v1(j),'l');%thr
    model_glu(j) = changeRxnBounds(model_glu(j),'EX0023',-0.0059*v1(j),'l');%val
    model_glu(j) = changeRxnBounds(model_glu(j),'EX0024',-0.0023*v1(j),'l');%tyr
    model_glu(j) = changeRxnBounds(model_glu(j),'T0036',0,'b');
    model_glu(j) = changeRxnBounds(model_glu(j),'T0037',0,'b');

    model_glu(j) = changeRxnBounds(model_glu(j),'EX0001',-v1(j),'b');
    model_glu(j) = changeRxnBounds(model_glu(j),'EXBiomass',miug(j),'b');
    model_glu(j) = changeObjective(model_glu(j),'prot_pool_exchange');
    FBAsolution = optimizeCbModel(model_glu(j),'max','one');
    FBAg(:,j)=FBAsolution.x;
    pro(j)=FBAsolution.x(2459);
    vlg(j)=FBAsolution.x(3);
    vace(j)=FBAsolution.x(69);
    vcit(j)=FBAsolution.x(70);
    vpyr(j)=FBAsolution.x(71);
    vacetoin(j)=FBAsolution.x(28);
    dcwg(j+1)=exp(miug(j)*1+log(dcwg(j)));% initial glucose concentration for the next step (mmol/L)
    Sg(j+1)=Sg(j)-v1(j)/miug(j)*(dcwg(j+1)-dcwg(j));% initial biomass fo the next step (g/L)
    Slg(j+1)=Slg(j)+vlg(j)/miug(j)*(dcwg(j+1)-dcwg(j));% initial lactate concentration for the next step (mmol/L)
    v1(j+1)=vmaxg*Sg(j+1)/(Kmg+Sg(j+1));% glucose specific consumption rate for the next step (mmol/gDCW/h)
end
result_glc2.rxns = ecModel.rxns;
result_glc2.x = FBAsolution.x;
%% ecMOMA
model_tre(1) = ecModel;
vmaxt=0.53;%mmol/gDCW/h
Kmt=15;%mmol/L海藻糖米氏常数 
St=zeros(69,1);
St(1)=57.58;%mmol/L海藻糖浓度
Slt=zeros(69,1);
Slt(1)=0;%mmol/L乳酸浓度
dcwt=zeros(69,1);
dcwt(1)=1.95;%g/L菌体干重 
v2=zeros(69,1);
v2(1)=vmaxt*St(1)/(Kmt+St(1));%mmol/gDCW/h海藻糖消耗速率
miut=zeros(69,1);
vlt=zeros(69,1);
vacet=zeros(69,1);
Sacet=zeros(69,1);
Sacet(1)=32.567;
MOMAt=zeros(2460,69);
model_tre(1) = changeRxnBounds(model_tre(1),'EX0064',-v2(1),'b');%tre

model_tre(1) = changeRxnBounds(model_tre(1),'EX0014',0.0417*v2(1),'l');%ala
model_tre(1) = changeRxnBounds(model_tre(1),'EX0015',-0.1,'l');%arg
model_tre(1) = changeRxnBounds(model_tre(1),'EX0016',-0.1,'l');%cys
model_tre(1) = changeRxnBounds(model_tre(1),'EX0017',-0.0114*v2(1),'l');%glu
model_tre(1) = changeRxnBounds(model_tre(1),'EX0018',-0.0069*v2(1),'l');%his
model_tre(1) = changeRxnBounds(model_tre(1),'EX0019',-0.0083*v2(1),'l');%ile
model_tre(1) = changeRxnBounds(model_tre(1),'EX0020',-0.0052*v2(1),'l');%leu
model_tre(1) = changeRxnBounds(model_tre(1),'EX0021',0.042*v2(1),'l');%met
model_tre(1) = changeRxnBounds(model_tre(1),'EX0022',-0.0351*v2(1),'l');%thr
model_tre(1) = changeRxnBounds(model_tre(1),'EX0023',-0.0081*v2(1),'l');%val
model_tre(1) = changeRxnBounds(model_tre(1),'EX0024',-0.0089*v2(1),'l');%tyr
model_tre(1) = changeRxnBounds(model_tre(1),'T0036',0,'b');
model_tre(1) = changeRxnBounds(model_tre(1),'T0037',0,'b');
model_tre(1) = changeRxnBounds(model_tre(1),'T0060',0,'b');
model_tre(1) = changeRxnBounds(model_tre(1),'T0061',0,'b');
model_tre(1) = changeRxnBounds(model_tre(1),'T0062_EXP_1',0,'b');
model_tre(1) = changeRxnBounds(model_tre(1),'T0062_EXP_2',0,'b');
model_tre(1) = changeRxnBounds(model_tre(1),'T0063_EXP_1',0,'b');
model_tre(1) = changeRxnBounds(model_tre(1),'T0063_EXP_2',0,'b');
model_tre(1) = changeRxnBounds(model_tre(1),'T0063_EXP_3',0,'b');
model_tre(1) = changeRxnBounds(model_tre(1),'T0063_EXP_4',0,'b');
model_tre(1) = changeRxnBounds(model_tre(1),'T0064_EXP_1',0,'b');
model_tre(1) = changeRxnBounds(model_tre(1),'T0064_EXP_2',0,'b');
model_tre(1) = changeRxnBounds(model_tre(1),'T0064_EXP_3',0,'b');
model_tre(1) = changeRxnBounds(model_tre(1),'T0064_EXP_4',0,'b');
model_tre(1) = changeRxnBounds(model_tre(1),'T0064_REV_EXP_1',0,'b');
model_tre(1) = changeRxnBounds(model_tre(1),'T0064_REV_EXP_2',0,'b');
model_tre(1) = changeRxnBounds(model_tre(1),'T0064_REV_EXP_3',0,'b');
model_tre(1) = changeRxnBounds(model_tre(1),'T0064_REV_EXP_4',0,'b');
model_tre(1) = changeRxnBounds(model_tre(1),'T0065',0,'b');
model_tre(1) = changeRxnBounds(model_tre(1),'T0066_EXP_1',0,'b');
model_tre(1) = changeRxnBounds(model_tre(1),'T0066_EXP_2',0,'b');
model_tre(1) = changeRxnBounds(model_tre(1),'T0066_EXP_3',0,'b');
model_tre(1) = changeRxnBounds(model_tre(1),'T0066_EXP_3',0,'b');
model_tre(1) = changeRxnBounds(model_tre(1),'T0007_REV',0,'b');
model_tre(1) = changeRxnBounds(model_tre(1),'T0066_EXP_4',v2(1),'b');
model_tre(1) = changeRxnBounds(model_tre(1),'NGAM',0,'l');
model_tre(1) = changeRxnBounds(model_tre(1),'EX0001',0,'b');

model_tre(1) = changeRxnBounds(model_tre(1),'EXBiomass',0,'b');
solMOMA = ecMOMA(model_tre(1),result_glc2);
MOMAt(:,1)=solMOMA.x;
result_tre(1).rxns = ecModel.rxns;
result_tre(1).x = solMOMA.x;
vacet(1)=solMOMA.x(69);
vlt(1)=solMOMA.x(3);
miut(1)=solMOMA.x(86);
St(2)=St(1)-v2(1)*dcwt(1)*1;%第一间隔的海藻糖初浓度
dcwt(2)=exp(miut(1)*1+log(dcwt(1)));%第一间隔的初始菌浓
Slt(2)=Slt(1)+vlt(1)*dcwt(1);%第一间隔的乳酸初浓度
Sacet(2)=Sacet(1)+vacet(1)*dcwt(1);
v2(2)=vmaxt*St(2)/(Kmt+St(2));
for i=1:67
    model_tre(i+1) = ecModel;
    model_tre(i+1) = changeRxnBounds(model_tre(i+1),'EX0064',-v2(i+1),'b');%tre

    model_tre(i+1) = changeRxnBounds(model_tre(i+1),'EX0014',0.0417*v2(i+1),'l');%ala
    model_tre(i+1) = changeRxnBounds(model_tre(i+1),'EX0015',-0.1,'l');%arg
    model_tre(i+1) = changeRxnBounds(model_tre(i+1),'EX0016',-0.1,'l');%cys
    model_tre(i+1) = changeRxnBounds(model_tre(i+1),'EX0017',-0.0114*v2(i+1),'l');%glu
    model_tre(i+1) = changeRxnBounds(model_tre(i+1),'EX0018',-0.0069*v2(i+1),'l');%his
    model_tre(i+1) = changeRxnBounds(model_tre(i+1),'EX0019',-0.0083*v2(i+1),'l');%ile
    model_tre(i+1) = changeRxnBounds(model_tre(i+1),'EX0020',-0.0052*v2(i+1),'l');%leu
    model_tre(i+1) = changeRxnBounds(model_tre(i+1),'EX0021',0.042*v2(i+1),'l');%met
    model_tre(i+1) = changeRxnBounds(model_tre(i+1),'EX0022',-0.0351*v2(i+1),'l');%thr
    model_tre(i+1) = changeRxnBounds(model_tre(i+1),'EX0023',-0.0081*v2(i+1),'l');%val
    model_tre(i+1) = changeRxnBounds(model_tre(i+1),'EX0024',-0.0089*v2(i+1),'l');%tyr
    model_tre(i+1) = changeRxnBounds(model_tre(i+1),'T0036',0,'b');
    model_tre(i+1) = changeRxnBounds(model_tre(i+1),'T0037',0,'b');
    model_tre(i+1) = changeRxnBounds(model_tre(i+1),'T0060',0,'b');
    model_tre(i+1) = changeRxnBounds(model_tre(i+1),'T0061',0,'b');
    model_tre(i+1) = changeRxnBounds(model_tre(i+1),'T0062_EXP_1',0,'b');
    model_tre(i+1) = changeRxnBounds(model_tre(i+1),'T0062_EXP_2',0,'b');
    model_tre(i+1) = changeRxnBounds(model_tre(i+1),'T0063_EXP_1',0,'b');
    model_tre(i+1) = changeRxnBounds(model_tre(i+1),'T0063_EXP_2',0,'b');
    model_tre(i+1) = changeRxnBounds(model_tre(i+1),'T0063_EXP_3',0,'b');
    model_tre(i+1) = changeRxnBounds(model_tre(i+1),'T0063_EXP_4',0,'b');
    model_tre(i+1) = changeRxnBounds(model_tre(i+1),'T0064_EXP_1',0,'b');
    model_tre(i+1) = changeRxnBounds(model_tre(i+1),'T0064_EXP_2',0,'b');
    model_tre(i+1) = changeRxnBounds(model_tre(i+1),'T0064_EXP_3',0,'b');
    model_tre(i+1) = changeRxnBounds(model_tre(i+1),'T0064_EXP_4',0,'b');
    model_tre(i+1) = changeRxnBounds(model_tre(i+1),'T0064_REV_EXP_1',0,'b');
    model_tre(i+1) = changeRxnBounds(model_tre(i+1),'T0064_REV_EXP_2',0,'b');
    model_tre(i+1) = changeRxnBounds(model_tre(i+1),'T0064_REV_EXP_3',0,'b');
    model_tre(i+1) = changeRxnBounds(model_tre(i+1),'T0064_REV_EXP_4',0,'b');
    model_tre(i+1) = changeRxnBounds(model_tre(i+1),'T0066_EXP_1',0,'b');
    model_tre(i+1) = changeRxnBounds(model_tre(i+1),'T0066_EXP_2',0,'b');
    model_tre(i+1) = changeRxnBounds(model_tre(i+1),'T0066_EXP_3',0,'b');
    model_tre(i+1) = changeRxnBounds(model_tre(i+1),'T0007_REV',0,'b');
    model_tre(i+1) = changeRxnBounds(model_tre(i+1),'T0066_EXP_4',v2(i+1),'b');
    model_tre(i+1) = changeRxnBounds(model_tre(i+1),'T0065',0,'b');
    model_tre(i+1) = changeRxnBounds(model_tre(i+1),'NGAM',0,'l');
    model_tre(i+1) = changeRxnBounds(model_tre(i+1),'EX0001',0,'b');

    model_tre(i+1) = changeRxnBounds(model_tre(i+1),'EXBiomass',0,'b');
    solMOMA = ecMOMA(model_tre(i+1),result_tre(i));
    MOMAt(:,i+1)=solMOMA.x;
    result_tre(i+1).rxns = ecModel.rxns;
    result_tre(i+1).x = solMOMA.x;
    vacet(i+1)=solMOMA.x(69);
    vlt(i+1)=solMOMA.x(3);
    miut(i+1)=solMOMA.x(86);
    St(i+2)=St(i+1)-v2(i+1)*dcwt(i+1)*1;%下一间隔的海藻糖初浓度
    dcwt(i+2)=exp(miut(i+1)*1+log(dcwt(i+1)));%下一间隔的初始菌浓
    Slt(i+2)=Slt(i+1)+vlt(i+1)*dcwt(i+1);%下一间隔的乳酸初浓度
    Sacet(i+2)=Sacet(i+1)+vacet(i+1)*dcwt(i+1);
    v2(i+2)=vmaxt*St(i+2)/(Kmt+St(i+2));
end