%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   3. SIMULATION RUN PLOT FILE:                          %
%               FUNCTION FILE FOR RUNNING MULTIPLE CASES                  %                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Subroutine Explanation
% This function runs each case and calls the main program. It gives proper
% inputs to the program, gets the outputs, compresses the outputs from
% multiple runs into one vector. It returns these to the first program, as
% well as, making a table for .dat/.txt file or outputting to an .xlsx file

function Run=SimulationRunPlot(I,K,MT,ML,D_k,CCPFail,RunNum,Scenario,Collusion,CCPTrade,DISTTrade)
%% Function Call
    for i=1:14
    
    if RunNum==1
        %Run 1: Collusion, Increasing Predators, Decreasing Distressed Banks
        Pred = i;
        Dist = I-Pred;
    elseif RunNum==2
        %Run 2: No Collusion, Increasing Predators, Decreasing Distressed Banks
        Pred = i;
        Dist = I-Pred;
    elseif RunNum==3
        %Run 3: Stable Distressed Banks, Increasing Predators Banks
        Pred = i;
        Dist = 2;
    elseif RunNum==4
        %Run 4: Stable Predator Banks, Increasing Distressed Banks
        Pred = 2;
        Dist = i;
    elseif RunNum==5
        %Run 5: CCP cannot sell in buyback period, Increasing Predators, Decreasing Distressed Banks
        Pred = i;
        Dist = I-Pred;
    elseif RunNum==6
       %Run 6: Distressed banks can sell in buyback period, Increasing Predators, Decreasing Distressed Banks
        Pred = i;
        Dist = I-Pred;
    elseif RunNum==7
        %Run 7: CCP cannot sell and distressed banks can in buyback period, Increasing Predators, Decreasing Distressed Banks
        Pred = i;
        Dist = I-Pred;
    end
            
    
    Values=SimulationAutomated(I,K,MT,ML,D_k,CCPFail,Pred,Dist,Collusion,CCPTrade,DISTTrade)
   
    %Want to capture a array of values for each run (over k)
    Xmarket(:,:,i)=Values(1); 
    N(:,:,i)=Values(2);
    Predator(:,:,i)=Values(3);
    Distress(:,:,i)=Values(4);
    DfltCount(:,:,i)=Values(5);
    %Liquidation Values
    LHCCPLoss(:,:,i)=Values(6);
    LPCCPLoss(:,:,i)=Values(7);
    lell(:,:,i)=Values(8);
    pell(:,:,i)=Values(9);
    hell(:,:,i)=Values(10);
    LDflt(:,:,i)=Values(11);
    LTHC_0(:,:,i)=Values(12);
    LTPC_0(:,:,i)=Values(13);
    LTHC_i(:,:,i)=Values(14);
    LTPC_i(:,:,i)=Values(15);
    LHLiquitationLoss(:,:,i)=Values(16);
    LPLiquidationLoss(:,:,i)=Values(17);
    LHLiquitationGain(:,:,i)=Values(18);
    LPLiquidationGain(:,:,i)=Values(19);
    LEHG_tot(:,:,i)=Values(20);
    LEPG_tot(:,:,i)=Values(21);
    %Buyback Values
    BHCCPLoss(:,:,i)=Values(22);
    BPCCPLoss(:,:,i)=Values(23);
    bhell(:,:,i)=Values(24);
    bpell(:,:,i)=Values(25);
    BDefault(:,:,i)=Values(26);
    AvgBuybackProfitsVF(:,:,i)=Values(27);
    AvgBuybackProfitsV0(:,:,i)=Values(28);
    bell(:,:,i)=Values(29);
    sell(:,:,i)=Values(30);
    BTHC_0(:,:,i)=Values(31);
    BTPC_0(:,:,i)=Values(32);
    BTHC_i(:,:,i)=Values(33);
    BTPC_i(:,:,i)=Values(34);
    BEHG_tot(:,:,i)=Values(35);
    BEPG_tot(:,:,i)=Values(36);
    BEHD_tot(:,:,i)=Values(37);
    BEPD_tot(:,:,i)=Values(38);
            % ValuesB3=[BHLiquidationLoss;BPLiquidationLoss;BHLiquidationGain;BPLiquidationGain;AvgPGR;AvgPCPlusLoss];
            % ValuesR1=[RTHC_0;RTPC_0;RTHC_i;RTPC_i;AvgRTHC; AvgRTPC;REHG_tot;REPG_tot;REHD_tot;REPD_tot];
            % ValuesR2=[RNewTHC_i;RNewTPC_i;AvgRNewTHC_i;AvgRNewTPC_i;LSurplus;BSurplus];
            % ValuesR3=[PredHProfitLoss;PredPProfitLoss];
            % Values=[Values0;Values1;ValuesL1;ValuesL2;ValuesB1;ValuesB2;ValuesB3;ValuesR1;ValuesR2;ValuesR3];
    BHLiquitationLoss(:,:,i)=Values(39);
    BPLiquidationLoss(:,:,i)=Values(40);
    BHLiquitationGain(:,:,i)=Values(41);
    BPLiquidationGain(:,:,i)=Values(42);
    AvgHGR(:,:,i)=Values(43);
    AvgHCPlusLoss(:,:,i)=Values(44)
    %Recovery Values
    RTHC_0(:,:,i)=Values(45);
    RTPC_0(:,:,i)=Values(46);
    RTHC_i(:,:,i)=Values(47);
    RTPC_i(:,:,i)=Values(48);
    AvgRTHC_i(:,:,i)=Values(49); 
    AvgRTPC_i(:,:,i)=Values(50);
    REHG_tot(:,:,i)=Values(51);
    REPG_tot(:,:,i)=Values(52);
    REHD_tot(:,:,i)=Values(53);
    REPD_tot(:,:,i)=Values(54);
    RNewTHC_i(:,:,i)=Values(55)
    RNewTPC_i(:,:,i)=Values(56)
    AvgRNewTHC_i(:,:,i)=Values(57);
    AvgRNewTPC_i(:,:,i)=Values(58);
    LSurplus(:,:,i)=Values(59);
    BSurplus(:,:,i)=Values(60);
    PredHProfitLoss(:,:,i)=Values(61);
    PredPProfitLoss(:,:,i)=Values(62);
    end
    
    %SQUEEZED VARIABLES FOR PLOTTING
    Xmarket=squeeze(Xmarket(1,1,:)); 
    N=squeeze(N(1,1,:));
    Predator=squeeze(Predator(1,1,:));
    Distress=squeeze(Distress(1,1,:));
    DfltCount=squeeze(DfltCount(1,1,:));
    %Liquidation Values
    LHCCPLoss=squeeze(LHCCPLoss(1,1,:));
    LPCCPLoss=squeeze(LPCCPLoss(1,1,:));
    lell=squeeze(lell(1,1,:));
    pell=squeeze(pell(1,1,:));
    hell=squeeze(hell(1,1,:));
    LDflt=squeeze(LDflt(1,1,:));
    LTHC_0=squeeze(LTHC_0(1,1,:));
    LTPC_0=squeeze(LTPC_0(1,1,:));
    LTHC_i=squeeze(LTHC_i(1,1,:));
    LTPC_i=squeeze(LTPC_i(1,1,:));
    LHLiquitationLoss=squeeze(LHLiquitationLoss(1,1,:));
    LPLiquidationLoss=squeeze(LPLiquidationLoss(1,1,:));
    LHLiquitationGain=squeeze(LHLiquitationGain(1,1,:));
    LPLiquidationGain=squeeze(LPLiquidationGain(1,1,:));
    LEHG_tot=squeeze(LEHG_tot(1,1,:));
    LEPG_tot=squeeze(LEPG_tot(1,1,:));
    %Buyback Values
    BHCCPLoss=squeeze(BHCCPLoss(1,1,:));
    BPCCPLoss=squeeze(BPCCPLoss(1,1,:));
    bhell=squeeze(bhell(1,1,:));
    bpell=squeeze(bpell(1,1,:));
    BDefault=squeeze(BDefault(1,1,:));
    AvgBuybackProfitsVF=squeeze(AvgBuybackProfitsVF(1,1,:));
    AvgBuybackProfitsV0=squeeze(AvgBuybackProfitsV0(1,1,:));
    bell=squeeze(bell(1,1,:));
    sell=squeeze(sell(1,1,:));
    BTHC_0=squeeze(BTHC_0(1,1,:));
    BTPC_0=squeeze(BTPC_0(1,1,:));
    BTHC_i=squeeze(BTHC_i(1,1,:));
    BTPC_i=squeeze(BTPC_i(1,1,:));
    BEHG_tot=squeeze(BEHG_tot(1,1,:));
    BEPG_tot=squeeze(BEPG_tot(1,1,:));
    BEHD_tot=squeeze(BEHD_tot(1,1,:));
    BEPD_tot=squeeze(BEPD_tot(1,1,:));
    BHLiquitationLoss=squeeze(BHLiquitationLoss(1,1,:));
    BPLiquidationLoss=squeeze(BPLiquidationLoss(1,1,:));
    BHLiquitationGain=squeeze(BHLiquitationGain(1,1,:));
    BPLiquidationGain=squeeze(BPLiquidationGain(1,1,:));
    AvgHGR=squeeze(AvgHGR(1,1,:))
    AvgHCPlusLoss=squeeze(AvgHCPlusLoss(1,1,:))
    %Recovery Values
    RTHC_0=squeeze(RTHC_0(1,1,:));
    RTPC_0=squeeze(RTPC_0(1,1,:));
    RTHC_i=squeeze(RTHC_i(1,1,:));
    RTPC_i=squeeze(RTPC_i(1,1,:));
    AvgRTHC_i=squeeze(AvgRTHC_i(1,1,:));
    AvgRTPC_i=squeeze(AvgRTPC_i(1,1,:));
    REHG_tot=squeeze(REHG_tot(1,1,:));
    REPG_tot=squeeze(REPG_tot(1,1,:));
    REHD_tot=squeeze(REHD_tot(1,1,:));
    REPD_tot=squeeze(REPD_tot(1,1,:));
    RNewTHC_i=squeeze(RNewTHC_i(1,1,:));
    RNewTPC_i=squeeze(RNewTPC_i(1,1,:));
    AvgRNewTHC_i=squeeze(AvgRNewTHC_i(1,1,:));
    AvgRNewTPC_i=squeeze(AvgRNewTPC_i(1,1,:));
    AggHProfitLoss_i=(RTHC_i-RNewTHC_i)./RTHC_i; %loss calculation
    AggPProfitLoss_i=(RTPC_i-RNewTPC_i)./RTPC_i;
    AvgHProfitLoss_i=(AvgRTHC_i-AvgRNewTHC_i)./AvgRTHC_i;
    AvgPProfitLoss_i=(AvgRTPC_i-AvgRNewTPC_i)./AvgRTPC_i;
    LSurplus=squeeze(LSurplus(1,1,:));
    BSurplus=squeeze(BSurplus(1,1,:));
    PredHProfitLoss=squeeze(PredHProfitLoss(1,1,:));
    PredPProfitLoss=squeeze(PredPProfitLoss(1,1,:));
    
    %This is the output for the function returned to program 1
    Run=[Xmarket,N,Predator,Distress,DfltCount,LHCCPLoss,LPCCPLoss,lell,pell,hell,...
        LDflt,LTHC_0,LTPC_0,LTHC_i,LTPC_i,LHLiquitationLoss,LPLiquidationLoss,LHLiquitationGain,LPLiquidationGain,LEHG_tot,...
        LEPG_tot,BHCCPLoss,BPCCPLoss,bhell,bpell,BDefault,AvgBuybackProfitsVF,AvgBuybackProfitsV0,bell,sell,...
        BTHC_0,BTPC_0,BTHC_i,BTPC_i,BEHG_tot,BEPG_tot,BEHD_tot,BEPD_tot,BHLiquitationLoss,BPLiquidationLoss,...
        BHLiquitationGain,BPLiquidationGain,AvgHGR,AvgHCPlusLoss,RTHC_0,RTPC_0,RTHC_i,RTPC_i,AvgRTHC_i,AvgRTPC_i,...
        REHG_tot,REPG_tot,REHD_tot,REPD_tot,RNewTHC_i,RNewTPC_i,AvgRNewTHC_i,AvgRNewTPC_i,AggHProfitLoss_i,AggPProfitLoss_i,...
        AvgHProfitLoss_i,AvgPProfitLoss_i,LSurplus,BSurplus,PredHProfitLoss,PredPProfitLoss]
    %Make data table
    fulldata=table(Xmarket,N,Predator,Distress,DfltCount,LHCCPLoss,LPCCPLoss,lell,pell,hell,...
        LDflt,LTHC_0,LTPC_0,LTHC_i,LTPC_i,LHLiquitationLoss,LPLiquidationLoss,LHLiquitationGain,LPLiquidationGain,LEHG_tot,...
        LEPG_tot,BHCCPLoss,BPCCPLoss,bhell,bpell,BDefault,AvgBuybackProfitsVF,AvgBuybackProfitsV0,bell,sell,...
        BTHC_0,BTPC_0,BTHC_i,BTPC_i,BEHG_tot,BEPG_tot,BEHD_tot,BEPD_tot,BHLiquitationLoss,BPLiquidationLoss,...
        BHLiquitationGain,BPLiquidationGain,AvgHGR,AvgHCPlusLoss,RTHC_0,RTPC_0,RTHC_i,RTPC_i,AvgRTHC_i,AvgRTPC_i,...
        REHG_tot,REPG_tot,REHD_tot,REPD_tot,RNewTHC_i,RNewTPC_i,AvgRNewTHC_i,AvgRNewTPC_i,AggHProfitLoss_i,AggPProfitLoss_i,...
        AvgHProfitLoss_i,AvgPProfitLoss_i,LSurplus,BSurplus,PredHProfitLoss,PredPProfitLoss)
    
    %Writes to excel file
    Name='fulldataexcel';
    End='.xlsx';
    Name=strcat(Name,Scenario); %concatenates the string to create an appropriate name for the file.
    filename = strcat(Name,End);
    writetable(fulldata,filename,'Sheet',RunNum,'Range','A1')
    
end
    
  
    
    
      

