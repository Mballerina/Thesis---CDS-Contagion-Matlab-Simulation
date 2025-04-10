%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   2. SIMULATION PLOT FILE:                              %
%                  FILE FOR CALLING PLOT CASES                            %                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Subroutine Explanation
% Subroutine prompts user for the market scenario that they want to run and
% returns the values to the previous program or returns a string that has
% all the subroutine choices.

function Call=SimulationCallPlot(choice)

    if choice=='N'
        %Determination of Market Liquidity Scenario:
        prompt = 'Which scenario do you want to run?: Stable (S), Decreasing (D), or Crisis (C) ... ';
        scenario = input(prompt, 's');
        if isempty(scenario)
            scenario = 'S';
        end
        Call=[scenario,'X','X','X']; % X will be seen as an empty case, proceed to finish.
    end

    if choice=='Y'
        Call=['S','D','C','X'];
    end
    
end