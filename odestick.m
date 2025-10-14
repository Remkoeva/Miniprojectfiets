function status = odestick(t,state,flag,varargin)

persistent h1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% handleiding voor gebruik van ODESTICK                                          %%
%%                                                                                %%
%% Om eenvoudige animatie te krijgen van een gesimuleerde beweging                %%
%% van een keten segmenten gedefinieerd conform Casius et al. 2005                %%
%% moet je de hieronder beschreven stappen doorlopen                              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% KvS 201709, version 3.0
% Uitgangspunt is dat je een FUNCTIE hebt die de differentiaalvergelijkingen beschrijft
% (zoals 'voordyn2' in de context van practicum 2)
% en dat je een SCRIPTFILE hebt waarmee, gebruikmakend van de matlab integrator ode113,
% een simulatie wordt uitgevoerd
% EN DAT UITVOERING VAN DE SCRIPT HET GEWENSTE RESULTAAT HEEFT
% (dwz dat eea al goed werkt)
%
% Wijzigingen aan te brengen in je SCRIPT
%    - de definitie van de parameters van de integrator moet in ieder geval het volgende bevatten
%      odeparms=('outputfcn',@odestick);
%      voorbeeld: had je een odeparms die er als volgt uitzag:
%      odeparms=('abstol',1e-8,'reltol',1e-8);
%      dan moet dat nu worden
%      odeparms=('abstol',1e-8,'reltol',1e-8,'outputfcn',@odestick);
%    - definieer VOORAFGAAND AAN DE AANROEP VAN ODE113 een substruct genaamd 'stick',
%      binnen de struct parms die al bestond; in deze parms.stick stoppen we 
%      instellingen voor de function die de animatie verzorgt. Deze struct moet er
%      als volgt uitzien: (maw je kunt onderstaande cut-and-pasten naar je SCRIPT en je hoeft
%      dan alleen te uncommenten en de waarden waar nodig aan te passen)

% plotparameters voor odestick
% parms.stick.dostick=1; % 1 for yes, 0 for no
% parms.stick.axisvector=[-1 1 -1 1]; % display area
% parms.stick.timestep=0.05; % animation timestep
% parms.stick.realtime_if_possible=1; % if 1 then wait for real time, if 0 then plot when ready
% parms.stick.fiindex=1; % vector of indices of segment angles in vector state that are to be plotted
% parms.stick.baseindex=NaN; % vector containing index of x and y base
%                            % position in vector state; supply NaN in case of fixed base


% 3. locatie van ODESTICK.m
%    - zorg dat deze mfile staat in een folder waar matlab deze kan vinden!

% preliminary work done at first call (flag = 'init')
if ~isempty(flag) & (flag=='init')
    if nargin>4
        ud=varargin{2};
    else
        ud.stick.dostick=1; % assume we want an animation
        ud.stick.axisvector=[-1 1 -1 1]; % display area
        ud.stick.timestep=0.05; % animation timestep
        ud.stick.realtime_if_possible=1; % if 1 then wait, if 0 then plot when ready
        ud.stick.fiindex=[1]; % index of first segment angle in vector state
        ud.stick.baseindex=NaN; % index of x base position in vector state (y is assumed to follow immediately) 
    end
    if ~isfield(ud.segparms,'L') % inputarg check
        status=1;disp('unable to create stick diagram'),return
    else
        h1=figure;%open new figure for stick
        set(h1,'UserData',ud);
        if ud.stick.dostick==false,status=0;title('parms.stick.dostick == false');drawnow;return,end
        L=ud.segparms.L;
        xcoor(1)=0;
        ycoor(1)=0;
        fi=state(ud.stick.fiindex); 
        for i=1:length(fi)
            xcoor(i+1)=xcoor(i)+L(i)*cos(fi(i));
            ycoor(i+1)=ycoor(i)+L(i)*sin(fi(i));
        end
        if ~isnan(ud.stick.baseindex)
            xcoor=xcoor+state(ud.stick.baseindex(1));
            ycoor=ycoor+state(ud.stick.baseindex(2));
        end
        
        hold off
        ud.line=plot(xcoor,ycoor);%,'EraseMode','Background');
        ud.lasttime=t(1);
        ud.lastplottime=clock;
        axis('equal')
        axis(ud.stick.axisvector);
        
        % The STOP button.
        ud.stop = 0;
        pos = get(0,'DefaultUicontrolPosition');
        pos(1) = pos(1) - 15;
        pos(2) = pos(2) - 15;
        str = 'ud=get(gcf,''UserData''); ud.stop=1; set(gcf,''UserData'',ud)';
        ud.hui=uicontrol( ...
            'Style','push', ...
            'String','Stop', ...
            'Position',pos, ...
            'Callback',str, ...
            'Tag','stop');
        
        set(h1,'UserData',ud);
    end
end

if ~isempty(flag) & (flag=='done')
    ud=get(h1,'userdata');
    if ud.stick.dostick==false,status=0;return,end
    % Change STOP button into ClOSE button.
    str = 'close';
    set(ud.hui,'String','Close','Callback',str,'Tag','close');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% flag = '', normal mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ud=get(h1,'userdata');
if ud.stick.dostick==false,status=0;return,end        

if (t-ud.lasttime)>ud.stick.timestep % time to make a plot
    L=ud.segparms.L;
    xcoor(1)=0;
    ycoor(1)=0;
    fi=state(ud.stick.fiindex);    
    for i=1:length(fi)
        xcoor(i+1)=xcoor(i)+L(i)*cos(fi(i));
        ycoor(i+1)=ycoor(i)+L(i)*sin(fi(i));
    end
    if ~isnan(ud.stick.baseindex)
        xcoor=xcoor+state(ud.stick.baseindex(1));
        ycoor=ycoor+state(ud.stick.baseindex(2));
    end
    
    if ud.stick.realtime_if_possible
        while etime(clock,ud.lastplottime)<ud.stick.timestep,end % wait for real time to catch up
    end    
    set(ud.line,'Xdata',xcoor,'Ydata',ycoor);
    ud.lasttime=t;
    ud.lastplottime=clock;
    set(gcf,'userdata',ud);
end

if ud.stop==0
    status=0;
else
    status=1;
    % Change STOP button into ClOSE button.
    str = 'close';
    set(ud.hui,'String','Close','Callback',str,'Tag','close');
end

drawnow
