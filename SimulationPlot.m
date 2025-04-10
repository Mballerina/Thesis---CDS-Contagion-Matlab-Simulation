%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    1. SIMULATION PLOT FILE:                             %
%                     FILE FOR PLOTTING CASES                             %                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Subroutine Explanation
% This subroutine calls the other subroutines and gives it the proper
% inputs in order to get back the information it needs to plot. First it
% calls SimulationCallPlot.m to see if it should run all market depth
% scenarios or just one. Then it calls SimulationRunPlot.m which will run
% all the cases and call the actually program to give outputs on each case.

%Clear all previous data to prevent program crashing
clear all

%% Determination of Total Auto Run or Separate Scenario Runs
%Prompts user to see if all cases should be run or just one particular
%market depth scenario.
prompt = 'Do you want to run all scenarios...(Y/N): ';
choice = input(prompt, 's');
    if isempty(choice)
        choice = 'Y'; %Runs all cases if there is no input
    end
    
scenario=SimulationCallPlot(choice) %Calls for four letter string
for index=1:4
    if scenario(index)~='X' %won't call a scenario if sees X
       Scenario=scenario(index); %sets the current scenario setting
        
        %% Setup Values:
        I=14;
        K=100;
        
        %% Case 1 : Stable Market Depth
        if Scenario == 'S'
            %Values for a stable market
            D_k= 221e9*ones(I,1,K); %stable financial market
            ML=0; %empty value or stable non-decreasing market
            MT=0*ones(I,1,K); %decreasing market liquidity
        end
        
        %% Case 2 : Decreasing Market Depth
        if Scenario == 'D'
            %Values for a stable and/or decreasing market
            D_k= 221e9*ones(I,1,K); %stable financial market
            ML=0; %empty value or stable non-decreasing market
            MT=38e9*ones(I,1,K); %decreasing market liquidity
        end
        
        %% Case 3 : Crisis Market Depth
        if Scenario == 'C'
            %Values for crisis market
            D_k= 12e9*ones(I,1,K); %financial crisis market
            ML='C'; %crisis value
            MT=0; %no-decreasing market liquidity
        end
        
        
        %% Run for 7 Scenarios:
        %Run 1: Collusion, Increasing Predators, Decreasing Distressed Banks
        RunNum=1;
        Collusion = 'Y'; %Collusion between multiple predators
        CCPFail='Y'; %stop when CCP fails
        CCPTrade='Y'; %CCP is allowed to traded into buyback period
        DISTTrade='N'; %distressed traders trading permission in buyback period.
        Run1=SimulationRunPlot(I,K,MT,ML,D_k,CCPFail,RunNum,Scenario,Collusion,CCPTrade,DISTTrade);
        
        %Run 2: No Collusion, Increasing Predators, Decreasing Distressed Banks
        RunNum=2;
        Collusion = 'N'; %Collusion between multiple predators
        CCPFail='Y'; %stop when CCP fails
        CCPTrade='Y'; %CCP is allowed to traded into buyback period
        DISTTrade='N'; %distressed traders trading permission in buyback period.
        Run2=SimulationRunPlot(I,K,MT,ML,D_k,CCPFail,RunNum,Scenario,Collusion,CCPTrade,DISTTrade);
        
        %Run 3: Stable Distressed Banks, Increasing Predators Banks
        RunNum=3;
        Collusion = 'N'; %Collusion between multiple predators
        CCPFail='Y'; %stop when CCP fails
        CCPTrade='Y'; %CCP is allowed to traded into buyback period
        DISTTrade='N'; %distressed traders trading permission in buyback period.
        Run3=SimulationRunPlot(I,K,MT,ML,D_k,CCPFail,RunNum,Scenario,Collusion,CCPTrade,DISTTrade);
        
        %Run 4: Stable Predator Banks, Increasing Distressed Banks
        RunNum=4;
        Collusion = 'N'; %Collusion between multiple predators
        CCPFail='Y'; %stop when CCP fails
        CCPTrade='Y'; %CCP is allowed to traded into buyback period
        DISTTrade='N'; %distressed traders trading permission in buyback period.
        Run4=SimulationRunPlot(I,K,MT,ML,D_k,CCPFail,RunNum,Scenario,Collusion,CCPTrade,DISTTrade);
        
        %Run 5: CCP cannot sell in buyback period, Increasing Predators, Decreasing Distressed Banks
        RunNum=5;
        Collusion = 'N'; %Collusion between multiple predators
        CCPFail='Y'; %stop when CCP fails
        CCPTrade='N'; %CCP is allowed to traded into buyback period
        DISTTrade='N'; %distressed traders trading permission in buyback period.
        Run5=SimulationRunPlot(I,K,MT,ML,D_k,CCPFail,RunNum,Scenario,Collusion,CCPTrade,DISTTrade);
        
        %Run 6: Distressed banks can sell in buyback period, Increasing Predators, Decreasing Distressed Banks
        RunNum=6;
        Collusion = 'N'; %Collusion between multiple predators
        CCPFail='Y'; %stop when CCP fails
        CCPTrade='Y'; %CCP is allowed to traded into buyback period
        DISTTrade='Y'; %distressed traders trading permission in buyback period.
        Run6=SimulationRunPlot(I,K,MT,ML,D_k,CCPFail,RunNum,Scenario,Collusion,CCPTrade,DISTTrade);
        
        %Run 7: CCP cannot sell and distressed banks can in buyback period, Increasing Predators, Decreasing Distressed Banks
        RunNum=7;
        Collusion = 'N'; %Collusion between multiple predators
        CCPFail='Y'; %stop when CCP fails
        CCPTrade='N'; %CCP is allowed to traded into buyback period
        DISTTrade='Y'; %distressed traders trading permission in buyback period.
        Run7=SimulationRunPlot(I,K,MT,ML,D_k,CCPFail,RunNum,Scenario,Collusion,CCPTrade,DISTTrade);
        
        
        %% Plot Case 1 : Stable Market Depth
        if Scenario == 'S'
            %% Figure2: Default Distribution Based On No. of Predatory vs. Distressed Banks
            figure
            hold off
            bardata2=[Run2(:,5)];
            bardata1=[Run2(:,4)];
            bardata=[bardata1, bardata2];
            b1=bar(bardata1);
            grid on
            %Data labels
            x1 = get(b1,'xdata');
            y1 = get(b1,'ydata');
            ygap1 = 0.6;  %// Specify vertical gap between the bar and label
            ylimits = get(gca,'ylim');
            set(gca,'ylim',[ylimits(1),ylimits(2)+0.2*max(y1)]); %extends y-axis
            labels1 = cellstr(num2str(y1')); %sets data labels
            for i = 1:length(x1) %// Loop over each bar
                %xpos1 = x1(i);
                xpos1 = i; %// Set x position for the text label
                ypos1 = y1(i) + ygap1; %// Set y position, including gap
                htext1 = text(xpos1,ypos1,labels1{i});          %// Add text label
                set(htext1,'VerticalAlignment','top', 'HorizontalAlignment','center');
            end
            
            hold on
            b2=bar(bardata2,'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1.5);
            %Data labels
            x2 = get(b2,'xdata');
            y2 = get(b2,'ydata');
            ygap2 =-0.2;  %// Specify vertical gap between the bar and label
            ylimits = get(gca,'ylim');
            set(gca,'ylim',[ylimits(1),ylimits(2)+0.2*max(y2)]);
            labels2 = cellstr(num2str(y2'));
            for i = 1:length(x2) %// Loop over each bar
                xpos2 = i;        %// Set x position for the text label
                ypos2 = y2(i) + ygap2; %// Set y position, including gap
                htext2 = text(xpos2,ypos2,labels2{i});          %// Add text label
                set(htext2,'VerticalAlignment','top', 'HorizontalAlignment','center');
            end
            %b1(1).LineWidth = 2;
            %b(1).EdgeColor = 'red';
            title('Default Distribution Based On No. of Predatory vs. Distressed Banks')
            x=xlabel('No. of Predatory Banks')
            y=ylabel('No. of Distressed or Defaulted Banks')
            %legend('show')
            legend('Distressed Banks', 'Defaulted Banks', 'Predatory Banks')
            legend('Location','north')
            
            print -depsc fig2.eps
            hold off
            
            %% Figure3: Default Distribution Based On No. of Two Predatory vs. Distressed Banks
            figure
            hold off
            bardata2=[Run4(:,3)]; %Distressed bank scenario is output 4
            bardata1=[Run4(:,5)];
            bardata=[bardata1, bardata2];
            b1=bar(bardata1);
            grid on
            %Data labels
            x1 = get(b1,'xdata');
            y1 = get(b1,'ydata');
            ygap1 = 0.4;  %// Specify vertical gap between the bar and label
            ylimits = get(gca,'ylim');
            set(gca,'ylim',[ylimits(1),ylimits(2)+0.2*max(y1)]);
            labels1 = cellstr(num2str(y1'));
            for i = 1:length(x1) %// Loop over each bar
                %xpos1 = x1(i);
                xpos1 = i; %// Set x position for the text label
                ypos1 = y1(i) + ygap1; %// Set y position, including gap
                htext1 = text(xpos1,ypos1,labels1{i});          %// Add text label
                set(htext1,'VerticalAlignment','top', 'HorizontalAlignment','center');
            end
            
            hold on
            b2=bar(bardata2,'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1.5);
            b2.FaceAlpha=0.5;
            grid on
            %Data labels
            x2 = get(b2,'xdata');
            y2 = get(b2,'ydata');
            ygap2 =-0.2;  %// Specify vertical gap between the bar and label
            ylimits = get(gca,'ylim');
            set(gca,'ylim',[ylimits(1),ylimits(2)+0.2*max(y2)]);
            labels2 = cellstr(num2str(y2'));
            for i = 1:length(x2) %// Loop over each bar
                xpos2 = i;        %// Set x position for the text label
                ypos2 = y2(i) + ygap2; %// Set y position, including gap
                htext2 = text(xpos2,ypos2,labels2{i});          %// Add text label
                set(htext2,'VerticalAlignment','top', 'HorizontalAlignment','center');
            end
            %b1(1).LineWidth = 2;
            %b(1).EdgeColor = 'red';
            title('Default Distribution Based On No. of Two Predatory vs. Distressed Banks');
            x=xlabel('No. of Predatory Banks');
            y=ylabel('No. of Defaulted Banks');
            %legend('show')
            legend('Defaulted Banks', 'Distressed Banks', 'Predatory Banks');
            legend('Location','north');
            
            print -depsc fig3.eps;
            hold off
            
            %% Figure4: Final CCP Loss under Predation and Distress
            figure
            hold off
            bardata4=[Run4(:,45)];
            bardata3=[Run3(:,45)];
            bardata2=[Run2(:,45)];
            bardata1=[Run1(:,45)];
            bardata=[bardata1,bardata2,bardata3,bardata4];
            b=bar(bardata);
            grid on
            title('Final CCP Loss under Predation and Distress');
            x=xlabel('No. of Banks');
            y=ylabel('Final Value CCP (USD)');
            %legend('show')
            legend('Collusion, Increasing Predation', 'No Collusion, Increasing Predation', 'Stable Distress, Increasing Predation', 'Stable Predation, Increasing Distress');
            legend('boxoff');
            legend('Location','south');
            
            print -depsc fig4.eps
            
            %% Figure5: Final Predation Profits
            figure
            hold off
            bardata4=[Run4(:,28)];
            bardata3=[Run3(:,28)];
            bardata2=[Run2(:,28)];
            bardata1=[Run1(:,28)];
            bardata=[bardata1,bardata2,bardata3,bardata4];
            b=bar(bardata);
            grid on
            title('Final Predation Profits');
            x=xlabel('No. of Banks');
            y=ylabel('Final Value CCP (USD)');
            %legend('show')
            legend('Collusion, Increasing Predation', 'No Collusion, Increasing Predation', 'Stable Distress, Increasing Predation', 'Stable Predation, Increasing Distress');
            legend('boxoff');
            legend('Location','southwest');
            
            print -depsc fig5.eps %black and white save is -deps2
            
        end
        
        %% Plot Case 2 : Decreasing Market Depth
        if Scenario == 'D'
            %% Figure6: Default Distribution Based on No. of Predatory vs. Distressed Banks
            figure
            hold off
            bardata2=[Run2(:,5)];
            bardata1=[Run2(:,4)];
            bardata=[bardata1, bardata2];
            b1=bar(bardata1);
            %Data labels
            x1 = get(b1,'xdata');
            y1 = get(b1,'ydata');
            ygap1 = 0.8;  %// Specify vertical gap between the bar and label
            ylimits = get(gca,'ylim');
            set(gca,'ylim',[ylimits(1),ylimits(2)+0.2*max(y1)]);
            labels1 = cellstr(num2str(y1'));
            for i = 1:(length(x1)-1) %// Loop over each bar
                %xpos1 = x1(i);
                xpos1 = i; %// Set x position for the text label
                ypos1 = y1(i) + ygap1; %// Set y position, including gap
                htext1 = text(xpos1,ypos1,labels1{i});          %// Add text label
                set(htext1,'VerticalAlignment','top', 'HorizontalAlignment','center');
            end
            
            hold on
            b2=bar(bardata2,'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1.5);
            b2.FaceAlpha=0.5;
            %Data labels
            x2 = get(b2,'xdata');
            y2 = get(b2,'ydata');
            ygap2 =0.8;  %// Specify vertical gap between the bar and label
            ylimits = get(gca,'ylim');
            set(gca,'ylim',[ylimits(1),ylimits(2)+0.2*max(y2)]);
            labels2 = cellstr(num2str(y2'));
            for i = 1:length(x2) %// Loop over each bar
                xpos2 = i;        %// Set x position for the text label
                ypos2 = y2(i) + ygap2; %// Set y position, including gap
                htext2 = text(xpos2,ypos2,labels2{i},'Color', 'r');          %// Add text label
                set(htext2,'VerticalAlignment','top', 'HorizontalAlignment','center');
            end
            title('Default Distribution Based On No. of Predatory vs. Distressed Banks');
            x=xlabel('No. of Predatory Banks');
            y=ylabel('No. of Distressed or Defaulted Banks');
            %legend('show')
            legend('Distressed Banks', 'Defaulted Banks', 'Predatory Banks');
            legend('Location','north');
            grid on
            
            print -depsc fig6.eps
            hold off 
            
            %% Figure7: Default Distribution Based on Two Predatory vs. Distressed Banks
            figure
            hold off
            bardata1=[Run4(:,5)];
            bardata2=[Run4(:,3)];
            bardata=[bardata1, bardata2];
            b1=bar(bardata1);
            %Data labels
            x1 = get(b1,'xdata');
            y1 = get(b1,'ydata');
            ygap1 = 0.8;  %// Specify vertical gap between the bar and label
            ylimits = get(gca,'ylim');
            set(gca,'ylim',[ylimits(1),ylimits(2)+0.2*max(y1)]);
            labels1 = cellstr(num2str(y1'));
            for i = 1:length(x1) %// Loop over each bar
                %xpos1 = x1(i);
                xpos1 = i; %// Set x position for the text label
                ypos1 = y1(i) + ygap1; %// Set y position, including gap
                htext1 = text(xpos1,ypos1,labels1{i});          %// Add text label
                set(htext1,'VerticalAlignment','top', 'HorizontalAlignment','center');
            end
            
            hold on
            b2=bar(bardata2,'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1.5);
            b2.FaceAlpha=0.5; %opacity parameter
            %Data labels
            x2 = get(b2,'xdata');
            y2 = get(b2,'ydata');
            ygap2 =0.8;  %// Specify vertical gap between the bar and label
            ylimits = get(gca,'ylim');
            set(gca,'ylim',[ylimits(1),ylimits(2)+0.2*max(y2)]);
            labels2 = cellstr(num2str(y2'));
            for i = 1:(length(x2)) %// Loop over each bar (except last because zero)
                xpos2 = i;        %// Set x position for the text label
                ypos2 = y2(i) + ygap2; %// Set y position, including gap
                htext2 = text(xpos2,ypos2,labels2{i}, 'color', 'r');          %// Add text label
                set(htext2,'VerticalAlignment','top', 'HorizontalAlignment','center');
            end
            title('Default Distribution Based On No. of Predatory vs. Distressed Banks');
            x=xlabel('No. of Distressed Banks');
            y=ylabel('No. of Predatory or Defaulted Banks');
            legend('Defaulted Banks','Predatory Banks');
            legend('Location','north');
            grid on 
            
            print -depsc fig7.eps
            hold off
            
            %% Figure8: Default Distribution Based on Predatory vs. Two Distressed Banks
            figure
            hold off
            bardata1=[Run3(:,4)];
            bardata2=[Run3(:,5)];
            bardata=[bardata1, bardata2];
            b1=bar(bardata1);
            %Data labels
            x1 = get(b1,'xdata');
            y1 = get(b1,'ydata');
            ygap1 = 0.1;  %// Specify vertical gap between the bar and label
            ylimits = get(gca,'ylim');
            set(gca,'ylim',[ylimits(1),ylimits(2)+0.2*max(y1)]);
            labels1 = cellstr(num2str(y1'));
            for i = 1:length(x1) %// Loop over each bar
                %xpos1 = x1(i);
                xpos1 = i; %// Set x position for the text label
                ypos1 = y1(i) + ygap1; %// Set y position, including gap
                htext1 = text(xpos1,ypos1,labels1{i});          %// Add text label
                set(htext1,'VerticalAlignment','top', 'HorizontalAlignment','center');
            end
            
            hold on
            b2=bar(bardata2,'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1.5);
            b2.FaceAlpha=0.5; %opacity parameter
            %Data labels
            x2 = get(b2,'xdata');
            y2 = get(b2,'ydata');
            ygap2 =0.1;  %// Specify vertical gap between the bar and label
            ylimits = get(gca,'ylim')
            set(gca,'ylim',[ylimits(1),ylimits(2)+0.2*max(y2)]);
            labels2 = cellstr(num2str(y2'));
            for i = 1:length(x2) %// Loop over each bar
                xpos2 = i;        %// Set x position for the text label
                ypos2 = y2(i) + ygap2; %// Set y position, including gap
                htext2 = text(xpos2,ypos2,labels2{i},'color','r');          %// Add text label
                set(htext2,'VerticalAlignment','top', 'HorizontalAlignment','center')
            end
            %b1(1).LineWidth = 2;
            %b(1).EdgeColor = 'red';
            title('Default Distribution Based On No. of Predatory vs. Two Distressed Banks');
            x=xlabel('No. of Predatory Banks');
            y=ylabel('No. of Distressed and Defaulted Banks');
            %legend('show')
            legend('Defaulted Banks', 'Distressed Banks', 'Predatory Banks');
            legend('Location','north');
            grid on
            
            print -depsc fig8.eps
            hold off 
            
            %% Figure9: Final Loss for CCP Under Predation and Distress
            %(Loss at CCP Failue Due to No. Distressed Banks)
            figure
            hold off
            subplot(1,2,1)
            yyaxis left
            bardata1=[Run1(:,45)];
            lb121=bar(bardata1,'FaceColor',[.5 0 0]);
            %lb121.FaceAlpha=0.5
            y=ylabel('Final Value CCP (USD) with Collusion');
            yyaxis right
            bardata2=[Run2(:,45)];
            rb121=bar(bardata2,'FaceColor',[0 0 .5]);
            %rb121.FaceAlpha=0.5
            y=ylabel('Final Value CCP (USD) with No Collusion');
            h=title('Final CCP Loss');
            %P = get(h,'Position')
            %set(h,'Position',[P(1) P(2)+0.3 P(3)])
            x=xlabel('No. of Predatory Banks');
            legend('Collusion', 'No Collusion');
            legend('boxoff');
            legend('Location','southoutside');
            legend('Orientation','horizontal');
            subplot(1,2,2)
            yyaxis left
            bardata3=[Run3(:,45)];
            lb122=bar(bardata3,'FaceColor',[.5 0 .5]);%,'EdgeColor',[.9 0 .9],'LineWidth',1.5);
            lb122.FaceAlpha=0.5;
            y=ylabel('Final Value CCP (USD), Predation Stable');
            yyaxis right %secondary axis
            bardata4=[Run4(:,45)];
            rb122=bar(bardata4,'FaceColor',[0 .5 .5]);%,'EdgeColor',[0 .9 .9],'LineWidth',1.5);
            rb122.FaceAlpha=0.5;
            y=ylabel('Final Value CCP (USD), Distressed Stable');
            h=title('Final CCP Loss');
            x=xlabel('No. of Banks');
            %legend('show')
            legend('Stable Distress', 'Stable Predation');
            legend('boxoff');
            legend('Location','southoutside');
            legend('Orientation','horizontal');
            
            print -depsc fig9.eps
            
            %% EXTENDED PLOTS
            %% Figure13: Average Profit or Loss of All Banks in CCP
            figure
            hold off
            subplot(2,3,1)
            x=[Run1(:,3),Run2(:,3),Run3(:,3),Run4(:,3)];
            y=[Run1(:,4),Run2(:,4),Run3(:,4),Run4(:,4)];
            z=[Run1(:,61)*100,Run2(:,61)*100,Run3(:,61)*100,Run4(:,61)*100];
            p1=plot(z,'LineWidth',2);
            p1(1).LineWidth=6;
            p1(3).LineWidth=4;
            p1(1).Color='blue';
            p1(2).Color='cyan';
            p1(3).Color='magenta';
            p1(4).Color='green';
            grid on
            x=xlabel('No. of Banks');
            y=ylabel({'Percentage Profit or';'Loss of Final Income'});
            
            subplot(2,3,2)
            x=[Run1(:,3),Run2(:,3),Run3(:,3),Run4(:,3)];
            y=[Run1(:,4),Run2(:,4),Run3(:,4),Run4(:,4)];
            z=[Run1(:,61)*100,Run2(:,61)*100];
            p2=plot(z,'LineWidth',2);
            grid on
            p2(1).LineWidth=6;
            p2(1).Color='blue';
            p2(2).Color='cyan';
            title({'Average Profit or Loss of All Banks in CCP After Margin Call'; '       '}); %two line title
            x=xlabel('No. of Predatory Banks');
            y=ylabel({'Percentage Profit or';'Loss of Final Income'});
            %legend('show')
            h=legend('Collusion', 'No Collusion');
            set(h,'FontSize',8);
            legend('Location','southeast');
            legend('boxoff');
            
            subplot(2,3,3)
            x=[Run1(:,3),Run2(:,3),Run3(:,3),Run4(:,3)];
            y=[Run1(:,4),Run2(:,4),Run3(:,4),Run4(:,4)];
            z=[Run3(:,61)*100,Run4(:,61)*100];
            %plot3(x,y,z,'LineWidth',2)
            p3=plot(z,'LineWidth',2);
            p3(1).LineWidth=4;
            p3(1).Color='magenta';
            p3(2).Color='green';
            grid on
            x=xlabel('No. of Banks');
            y=ylabel({'Percentage Profit or';'Loss of Final Income'});
            %legend('show')
            h=legend('Stable Distressed', 'Stable Predators');
            set(h,'FontSize',8);
            legend('Location','southwest');
            legend('boxoff');
            
            %-----------------
            
            subplot(2,3,4)
            x=[Run1(:,3),Run2(:,3),Run3(:,3),Run4(:,3)];
            y=[Run1(:,4),Run2(:,4),Run3(:,4),Run4(:,4)];
            z=[Run1(:,61)*100,Run2(:,61)*100,Run3(:,61)*100,Run4(:,61)*100];
            %surf(x,y,z)
            b1=bar(z);
            grid on
            b1(1).FaceColor='blue';
            b1(2).FaceColor='cyan';
            b1(3).FaceColor='magenta';
            b1(4).FaceColor='green';
            %title({'Average Profit or Loss of All Banks in CCP'; '       '});
            x=xlabel('No. of Banks');
            y=ylabel({'Percentage Profit or';'Loss of Final Income'});
            xticks([2,4,6,8,10,12,14]);
            
            subplot(2,3,5)
            x=[Run1(:,3),Run2(:,3),Run3(:,3),Run4(:,3)];
            y=[Run1(:,4),Run2(:,4),Run3(:,4),Run4(:,4)];
            z=[Run1(:,61)*100,Run2(:,61)*100];
            %surf(x,y,z)
            b2=bar(z);
            grid on
            b2(1).FaceColor='blue';
            b2(2).FaceColor='cyan';
            %title({'Average Profit or Loss of All Banks in CCP'; '       '});
            x=xlabel('No. of Predatory Banks');
            y=ylabel({'Percentage Profit or';'Loss of Final Income'});
            xticks([2,4,6,8,10,12,14]);
            %legend('show')
            h=legend('Collusion', 'No Collusion');
            set(h,'FontSize',8);
            legend('Location','southeast');
            legend('boxoff');
            
            subplot(2,3,6)
            x=[Run1(:,3),Run2(:,3),Run3(:,3),Run4(:,3)];
            y=[Run1(:,4),Run2(:,4),Run3(:,4),Run4(:,4)];
            z=[Run3(:,61)*100,Run4(:,61)*100];
            %surf(x,y,z)
            b3=bar(z);
            grid on
            b3(1).FaceColor='magenta';
            b3(2).FaceColor='green';
            %title({'Average Profit or Loss of All Banks in CCP'; '       '});
            x=xlabel('No. of Banks');
            y=ylabel({'Percentage Profit or';'Loss of Final Income'});
            xticks([2,4,6,8,10,12,14]);
            %legend('show')
            h=legend('Stable Distressed', 'Stable Predators');
            set(h,'FontSize',8); 
            legend('Location','southwest');
            legend('boxoff');
            
            print -depsc fig13.eps
            
            %% Figure14: Final Percentage Profit Loss: Average Total Banks vs. Predators
            figure
            hold off
            subplot(1,3,1)
            linedata2=[Run2(:,61)*100,Run2(:,65)*100];
            l2=plot(linedata2,'LineWidth',3);
            l2(1).Color=[0 0.6 0.7];
            l2(2).LineStyle='-.'
            grid on
            title({'No Collusion,';'Decreasing Distressed'});
            x=xlabel('No. of Predatory Banks');
            y=ylabel('Percentage of Final Income (USD)');
            
            
            subplot(1,3,2)
            linedata3=[Run3(:,61)*100,Run3(:,65)*100];
            l3=plot(linedata3,'LineWidth',3);
            l3(1).Color=[0 0.6 0.7];
            l3(2).LineStyle='-.'
            ylim([31,36]);
            grid on
            title({'Final Percentage Profit/Loss: Bank Average vs. Predators';'';'Stable Distressed,'});
            x=xlabel('No. of Predator Banks');
            y=ylabel('Percentage of Final Income (USD)');
            %legend('show')
            legend('Avg. Bank Profits', 'Predator Profits');
            legend('boxoff');
            legend('Location','best');
            
            subplot(1,3,3)
            linedata4=[Run4(:,61)*100,Run4(:,65)*100];
            l4=plot(linedata4,'LineWidth',3);
            l4(1).Color=[0 0.6 0.7];
            l4(2).LineStyle='-.'
            ylim([-50,100]);
            yticks([-50 -25 0 25 50 75 100]);
            grid on
            title({''; 'Stable Predators'});
            x=xlabel('No. of Distressed Banks');
            y=ylabel('Percentage of Final Income (USD)');
            
            print -depsc fig14.eps
            
            %% Figure15: Buyback Profits vs. Margin Call and Buyback Profits vs. Percentage Income 
            figure
            hold off
            subplot(2,3,1)
            yyaxis left
            l1=plot(Run2(:,28),'LineWidth',3);
            grid on
            y=ylabel({'Buyback Earnings';' on Positions (USD)'});
            yyaxis right
            l2=plot(Run2(:,43),'LineWidth',3);
            l2.Color=[1 0.5 0];
            grid on
            y=ylabel('Margin Payment (USD)');
            title({'';'No Collusion'});
            x=xlabel('No. of Predatory Banks');
            
            subplot(2,3,2)
            yyaxis left
            l1=plot(Run3(:,28),'LineWidth',3);
            grid on
            y=ylabel({'Buyback Earnings';' on Positions (USD)'});
            yyaxis right
            l2=plot(Run3(:,43),'LineWidth',3);
            l2.Color=[1 0.5 0];
            grid on
            y=ylabel('Margin Payment (USD)');
            title({'Predator Profits vs. Margin Call:';'Stable Distressed'});
            x=xlabel('No. of Predatory Banks');
            
            subplot(2,3,3)
            linedata4=[Run4(:,28),Run4(:,43)];
            yyaxis left
            l1=plot(Run4(:,28),'LineWidth',3);
            grid on
            y=ylabel({'Buyback Earnings';' on Positions (USD)'});
            yyaxis right
            l2=plot(Run4(:,43),'LineWidth',3);
            l2.Color=[1 0.5 0];
            grid on
            y=ylabel('Margin Payment (USD)');
            title({''; 'Stable Predators'});
            x=xlabel('No. of Distressed Banks');
            
            subplot(2,3,4)
            yyaxis left
            l1=plot(Run2(:,28),'LineWidth',3);
            grid on
            y=ylabel({'Buyback Earnings';' on Positions (USD)'});
            yyaxis right
            l2=plot(Run2(:,65)*100,'LineWidth',3);
            l2.Color=[1 0.7 0];
            grid on
            y=ylabel('Percent Income Lost to Margin');
            title({'';'No Collusion'});
            x=xlabel('No. of Predatory Banks');
            
            subplot(2,3,5)
            yyaxis left
            l1=plot(Run3(:,28),'LineWidth',3);
            grid on
            y=ylabel({'Buyback Earnings';' on Positions (USD)'});
            yyaxis right
            l2=plot(Run3(:,65)*100,'LineWidth',3);
            l2.Color=[1 0.7 0];
            grid on
            y=ylabel('Percent Income Lost to Margin');
            title({'Predator Profits vs. Percentage Income Lost:';'Stable Distressed'});
            x=xlabel('No. of Predatory Banks');
            
            subplot(2,3,6)
            yyaxis left
            l1=plot(Run4(:,28),'LineWidth',3);
            grid on
            y=ylabel({'Buyback Earnings';' on Positions (USD)'});
            yyaxis right
            l2=plot(Run4(:,65)*100,'LineWidth',3);
            l2.Color=[1 0.7 0];
            grid on
            y=ylabel('Percent Income Lost to Margin');
            title({''; 'Stable Predators'});
            x=xlabel('No. of Distressed Banks');
            
            print -depsc fig15.eps
            
            %% Figure16: Predation Buyback Profit vs. Loss - Original vs. Buyback Value of Positions
            figure
            hold off
            linedata2=[Run1(:,28),Run2(:,28),Run3(:,28),Run4(:,28)];
            l2=plot(linedata2,'LineWidth',2);
            grid on
            xlim([0 15]);
            ylim([-2.5e9 4e9]);
            l2(1).LineWidth = 3;
            l2(1).LineStyle = '-.';
            l2(2).LineWidth = 3;
            l2(3).LineWidth = 4;
            l2(3).LineStyle = '--';
            l2(4).LineWidth = 4;
            title({'Predation Buyback Profit/Loss:';'Original Value vs. Buyback Value of Positions'});
            x=xlabel('No. of Banks');
            y=ylabel('Predator Buyback Profit/Loss (USD)');
            %legend('show')
            legend('Collusion', 'No Collusion', 'Stable Distressed','Stable Predators');
            legend('Location','best');
            legend('boxoff');
            
            print -depsc fig16.eps
            
            %% Figure17: Margin Refill Required by CCP in Recovery Stage
            figure
            hold off
            linedata2=[Run1(:,43),Run2(:,43),Run3(:,43),Run4(:,43)];
            l2=plot(linedata2,'LineWidth',2);
            grid on
            l2(1).LineWidth = 3;
            l2(1).LineStyle = '-.';
            l2(2).LineWidth = 3;
            l2(3).LineWidth = 4;
            l2(3).LineStyle = '--';
            l2(4).LineWidth = 4;
            
            title('Margin Refill Required By CCP in Recovery Stage');
            x=xlabel('No. of Banks');
            y=ylabel('Margin Refill Amount for Predators (USD)');
            %legend('show')
            legend('Collusion, Increasing Predators', 'No Collusion, Increasing Predators', 'Stable Distressed, Increasing Predators','Stable Predators, Increasing Distressed');
            legend('Location','best');
            legend('boxoff');
            
            print -depsc fig17.eps
            
            %% Figure18: Difference in Liquidation Loss: Pure vs. Hybrid: Pure vs. Hybrid
            figure
            hold off
            linedata2=[abs(Run3(:,40))-abs(Run3(:,39)),abs(Run4(:,40))-abs(Run4(:,39)),abs(Run2(:,40))-abs(Run2(:,39))];
            l2=plot(linedata2,'LineWidth',4);
            grid on
            l2(3).LineStyle = '-.';
            l2(2).LineWidth = 5;
            l2(2).LineStyle = '--';
            title({'Difference in CCP Liquidation Loss:';' Pure vs. Hybrid'});
            x=xlabel('No. of Banks');
            y=ylabel('Difference in Liquidation Loss (USD)');
            %legend('show')
            legend('Stable Distressed, Increasing Predators','Stable Predators, Increasing Distressed','No Collusion, Increasing Predators');
            legend('Location','best');
            legend('boxoff');
            
            print -depsc fig18.eps
            
            %% Figure19: Difference in Aggregate & Avg Final Value for Banks - Pure vs. Hybrid
            figure
            hold off
            subplot(1,2,1)
            linedata21=[Run3(:,48),Run3(:,47)];
            l2=plot(linedata21,'LineWidth',2);
            %l2(1).LineStyle='-.';
            grid on
            title({'Difference in Final Aggregate Bank Value:';''});
            x=xlabel('No. of Distressed Banks');
            y=ylabel('Difference in Final Income of All Banks (USD)');
            legend('Pure Fund','Hybrid Fund');
            legend('Location','south');
            legend('boxoff');
            
            subplot(1,2,2)
            linedata22=[Run3(:,50),Run3(:,49)];
            l22=plot(linedata22,'LineWidth',2);
            grid on
            %l22(1).LineStyle='-.';
            title({'Pure vs. Hybrid, Stable Predator'; ''});
            x=xlabel('No. of Distressed Banks');
            y=ylabel('Difference in Final Income of All Banks (USD)');
            %legend('show')
            legend('Pure Fund','Hybrid Fund');
            legend('Location','south');
            legend('boxoff');
            
            print -depsc fig19.eps
            
            %% Figure20: Liquidation and Buyback Surplus for Banks: Hybrid vs. Pure
            figure
            hold off
            subplot(1,3,1)
            linedata2=[Run2(:,63),Run2(:,64)];
            l2=plot(linedata2,'LineWidth',4);
            grid on
            %l2(2).LineWidth = 4;
            title({''; 'No';'Collusion'});
            x=xlabel('No. of Predatory Banks');
            y=ylabel('Difference in Surplus of All Banks (USD)');
            
            subplot(1,3,2)
            linedata22=[Run3(:,63),Run3(:,64)];
            l22=plot(linedata22,'LineWidth',4);
            grid on
            ylim([0 14e10]);
            %l22(2).LineWidth = 4;
            title({'Liquidation/Buyback Bank Surplus: Hybrid vs. Pure'; 'Stable'; 'Distressed'});
            x=xlabel('No. of Predatory Banks');
            y=ylabel('Difference in Surplus of All Banks (USD)');
            %legend('show')
            legend('Liquidation Surplus', 'Buyback Surplus');
            legend('Location','best');
            legend('boxoff');
            
            subplot(1,3,3)
            linedata23=[Run4(:,63),Run4(:,64)]
            l23=plot(linedata23,'LineWidth',4);
            grid on
            %l22(2).LineWidth = 4;
            title({''; 'Stable';'Predators'});
            x=xlabel('No. of Distressed Banks');
            y=ylabel('Difference in Surplus of All Banks (USD)');
            
            print -depsc fig20.eps
            
            %% Figure21: Increasing Loss for CCP and Banks in Different Buyback Strategies
            figure
            hold off
            subplot(1,2,1)
            linedata2=[Run2(:,45),Run5(:,45),Run6(:,45)];
            l2=bar(linedata2);
            ylim([-9e13 1e13]);
            grid on
            title({'CCP Final Loss:'; 'Different Buyback Strategies';''});
            x=xlabel('No. of Banks');
            y=ylabel('Final Income (USD)');
            %legend('show')
            legend('No Collusion','No CCP Selling', 'Distressed Bank Selling');
            legend('Location','south');
            legend('boxoff');
            
            subplot(1,2,2)
            linedata22=[Run2(:,47),Run5(:,47),Run6(:,47)];
            l22=bar(linedata22);
            ylim([-9e13 1e13]);
            grid on
            %l22(2).LineWidth = 4;
            title({'Bank Final Loss:'; 'Different Buyback Strategies';''});
            x=xlabel('No. of Banks');
            y=ylabel('Final Income (USD)');
            %legend('show')
            legend('No Collusion','No CCP Selling', 'Distressed Bank Selling');
            legend('Location','south');
            legend('boxoff');
         
            print -depsc fig21.eps
            
        end
        
        %% Plot Case 3 : Crisis Market Depth
        if Scenario == 'C'
            %% Figure10: Loss on Liquidation Day 1-5
            figure
            hold off
            bardata4=[-1*Run4(:,16)];
            bardata3=[-1*Run3(:,16)];
            bardata2=[-1*Run2(:,16)];
            bardata1=[-1*Run1(:,16)];
            bardata=[bardata1,bardata2,bardata3,bardata4];
            b=bar(bardata);
            grid on 
            title('CCP Loss during Liquidation Window (Days 1-5)');
            x=xlabel('No. of Banks');
            y=ylabel('Liquidation Loss for CCP (USD)');
            %legend('show')
            legend('Collusion, Increasing Predation', 'No Collusion, Increasing Predation', 'Stable Distress, Increasing Predation', 'Stable Predation, Increasing Distress');
            legend('boxoff');
            legend('Location','south');
            %legend('Orientation','horizontal')
            
            print -depsc fig10.eps
            %saveas(gca, 'fig10.eps','epsc');
            
            %% Figure11: Final CCP Loss under Predation and Distress
            %May need a top axis for distressed banks
            figure
            hold off
            bardata4=[Run4(:,45)];
            bardata3=[Run3(:,45)];
            bardata2=[Run2(:,45)];
            bardata1=[Run1(:,45)];
            bardata=[bardata1,bardata2,bardata3,bardata4];
            b=bar(bardata);
            title('Final CCP Loss under Predation and Distress');
            x=xlabel('No. of Banks');;
            %legend('show')
            legend('Collusion, Increasing Predation', 'No Collusion, Increasing Predation', 'Stable Distress, Increasing Predation', 'Stable Predation, Increasing Distress');
            legend('boxoff');
            legend('Location','south');
            %legend('Orientation','horizontal')
            
            print -depsc fig11.eps
            
            %% Figure12: Difference in Max. Possible Loss Pure vs. Hybrid
            figure
            hold off
            linedata4=[abs(Run4(:,46))-abs(Run4(:,45)+Run4(:,43))];
            linedata3=[abs(Run3(:,46))-abs(Run3(:,45)+Run3(:,43))];
            linedata2=[abs(Run2(:,46))-abs(Run2(:,45)+Run2(:,43))];
            linedata1=[abs(Run1(:,46))-abs(Run1(:,45)+Run1(:,43))];
            linedata=[linedata1,linedata2,linedata3,linedata4];
            l=bar(linedata);
            title('Difference in Max. CCP Possible Loss Pure vs. Hybrid');
            x=xlabel('No. of Banks');
            y=ylabel('Difference (Pure - Hybrid) in Final Loss for CCP (USD)');
            %legend('show')
            legend('Collusion, Increasing Predation', 'No Collusion, Increasing Predation', 'Stable Distress, Increasing Predation', 'Stable Predation, Increasing Distress');
            legend('boxoff');
            legend('Location','north');
            
            print -depsc fig12.eps
            
        end
    end
end
display('Simulation and Plotting is finished!') %lets you know when plots achieved

%% Reference Values:
%Run=[Xmarket-1,N-2,Predator-3,Distress-4,DfltCount-5,LHCCPLoss-6,LPCCPLoss-7,lell-8,pell-9,hell-10,LDflt-11,LTHC_0-12,LTPC_0-13,LTHC_i,LTPC_i,...
%        LHLiquitationLoss,LPLiquidationLoss,LHLiquitationGain,LPLiquidationGain,LEHG_tot,LEPG_tot,BHCCPLoss,BPCCPLoss,...
%        bhell,bpell,BDefault,AvgBuybackProfitsVF,AvgBuybackProfitsV0,bell,sell,BTHC_0,BTPC_0,BTHC_i,BTPC_i,...
%        BEHG_tot,BEPG_tot,BEHD_tot,BEPD_tot,BHLiquitationLoss,BPLiquidationLoss,BHLiquitationGain,BPLiquidationGain,...
%        AvgPGR,AvgPCPlusLoss,RTHC_0,RTPC_0,RTHC_i,RTPC_i,AvgRTHC_i,AvgRTPC_i,REHG_tot,REPG_tot,REHD_tot,REPD_tot,...
%        RNewTHC_i,RNewTPC_i,AvgRNewTHC_i,AvgRNewTPC_i,AggHProfitLoss_i,AggPProfitLoss_i,AvgHProfitLoss_i,AvgPProfitLoss_i,...
%        LSurplus,BSurplus,PredHProfitLoss,PredPProfitLoss]