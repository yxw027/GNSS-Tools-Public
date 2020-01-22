function plotAgainstTruth(outStruc,filenameTruth,PARAMS)


% parse the truth data
if ~exist('posTruth','var')
    [~,posTruth, epochsTruth,~,~,~,velEnu] = parsePosSolFile(filenameTruth);
    
    posMeas.pos = posTruth;
    posMeas.epochs = epochsTruth;
    
    posMeas.enuPos = nan(size(posTruth));
    
    llh0 = xyz2llh(posTruth(1,:));
    
    utmZone = findUtmZone(llh0(1),llh0(2));
    for idx = 1:size(posTruth,1)
        [posMeas.enuPos(idx,1),posMeas.enuPos(idx,2),posMeas.enuPos(idx,3)] = ...
            cart2utm(posTruth(idx,1),posTruth(idx,2),posTruth(idx,3),utmZone);
    end
end

%% Pull data out of the output saving object
pos = outStruc.posSave;
vel = outStruc.velSave;
att = outStruc.attSave;
attEnu = outStruc.attEnuSave;
biasSave = outStruc.biasSave;
epochSave = outStruc.epochSave;
clockSave = outStruc.clockSave;
covEnu = outStruc.covEnuFull';

%% Compare to truth data
if ~isempty(posMeas)
    % compute the output position to ENU
    posOutEnu = nan(size(pos));
    for idx = 1:size(pos,1)
        [posOutEnu(idx,1),posOutEnu(idx,2),posOutEnu(idx,3)] = ...
            cart2utm(pos(idx,1),pos(idx,2),pos(idx,3),utmZone);
    end
    
    % Convert velocity to ENU velocity
    velOutEnu = nan(size(pos));
    for idx = 1:size(pos,1)
        velOutEnu(idx,:) = XYZ2ENU(vel(idx,:),llh0(1)*pi/180,llh0(2)*pi/180);
    end
    
    % Interpolate truth positions to compare
    truthEnuInterp = nan(size(posOutEnu));
    if isempty(posMeas.epochs)
        truthEnuInterp = repmat(posMeas.enuPos,size(posOutEnu,1),1);
    else
        interpMethod = 'spline';
        truthEnuInterp(:,1) = interp1(posMeas.epochs,posMeas.enuPos(:,1),epochSave,interpMethod);
        truthEnuInterp(:,2) = interp1(posMeas.epochs,posMeas.enuPos(:,2),epochSave,interpMethod);
        truthEnuInterp(:,3) = interp1(posMeas.epochs,posMeas.enuPos(:,3),epochSave,interpMethod);
    end
    
    %     velOutEnu = [nan nan nan; diff(posOutEnu)./diff(epochSave)];
    if ~isempty(velEnu)
        truthVelEnu = nan(size(pos));
        truthVelEnu(:,1) = interp1(posMeas.epochs,velEnu(:,1),epochSave,interpMethod);
        truthVelEnu(:,2) = interp1(posMeas.epochs,velEnu(:,2),epochSave,interpMethod);
        truthVelEnu(:,3) = interp1(posMeas.epochs,velEnu(:,3),epochSave,interpMethod);

    else
        truthVelEnu = [nan nan nan; diff(truthEnuInterp)./diff(epochSave)];
        truthVelEnu = [nan nan nan; (truthEnuInterp(3:end,:)-truthEnuInterp(1:(end-2),:))./(epochSave(3:end)-epochSave(1:(end-2))); ...
            nan nan nan];
       
    end
     truthAttEnu = nan(size(truthVelEnu));
        for idx = 1:size(truthVelEnu,1)
            truthAttEnu(idx,:) = initializeAttitude(truthVelEnu(idx,:)',PARAMS);
        end
    
    errEnu = posOutEnu-truthEnuInterp;
    errVelEnu = velOutEnu -truthVelEnu;
    %
    epochsPlot = epochSave-epochSave(1);
    
    if 0
        indsGnss = find(epochSave == floor(epochSave));
    else
        indsGnss = find(epochSave == epochSave);
    end
    
    fig1 = figure; clf; 
    ha = tight_subplot(3,1,0.05,[0.1 0.1],[0.07 0.05]);
    labels = {'Easting [m]','Northing [m]','Upping [m]'};
    for idx = 1:3
        axes(ha(idx))
        plot(epochsPlot(indsGnss),abs(errEnu(indsGnss,idx)),'.');
                hold on;

%         semilogy(epochsPlot(indsGnss),sqrt(covEnu(indsGnss,idx))*3);
        plot(epochsPlot(indsGnss),outStruc.plFull(idx,indsGnss),'k.','linewidth',2);
        ylabel(labels{idx})
        grid on
%         ylim([0 5])
        if idx == 1
            title('Position Error in ENU') 
        end
    end
    legend('Error','Protection Level')
    xlabel('Time of run [s]')
    
    % Convert position error to body frame
    figure; clf; hold on;
    errBody = nan(size(errEnu));
    for idx = 1:max(size(errBody))
        Rbodyi = euler2dcm123(truthAttEnu(idx,:));
        
        errBody(idx,:) = Rbodyi'*-errEnu(idx,:)';
        
    end
    plot(errBody)
    legend('back','right','up')
    
    if 1
        fig2 = figure; clf; hold on;
        plot(truthEnuInterp(:,1)-truthEnuInterp(1,1),truthEnuInterp(:,2)-truthEnuInterp(1,2),'go')
        plot(posOutEnu(:,1)-truthEnuInterp(1,1),posOutEnu(:,2)-truthEnuInterp(1,2),'.')
        
        distanceTraveledTruth = sum(sqrt(sum(diff(truthEnuInterp).^2,2)));
        distanceTraveledEst = sum(sqrt(sum(diff(posOutEnu).^2,2)));
        
        disp(distanceTraveledTruth)
        disp(distanceTraveledEst)
    end
    
    if 1
        fig3 = figure; clf; hold on;
        plot(epochsPlot,attEnu(:,1)*180/pi,'b')
        plot(epochsPlot,attEnu(:,2)*180/pi,'g')
        plot(epochsPlot,attEnu(:,3)*180/pi,'r')
        
        velMag = sqrt(sum(truthVelEnu.^2,2));
        
        truthAttEnu(repmat(velMag',1,3) < 1) = nan;
        
        plot(epochsPlot,truthAttEnu(:,1)*180/pi,'b--')
        plot(epochsPlot,truthAttEnu(:,2)*180/pi,'g--')
        plot(epochsPlot,truthAttEnu(:,3)*180/pi,'r--')
        plot(epochsPlot,velMag*10,'k','linewidth',2)
    end
    
    fig4 = figure; clf; hold on;
    plot(epochsPlot,biasSave-nanmedian(biasSave)*0)
        
    if 1
        fig5= figure; clf; hold on;
        plot(epochsPlot(indsGnss),errVelEnu(indsGnss,:))
        grid on;
        title('Velocity error in ENU')
    end
end



























end