function draw2(singleLineTrajectory, dataResolution, drawBodyCurve, drawTrajectory, doExport, framerate, clearFigure)
startTime = tic;
config = load('configTim.mat');
config = config.config;
drawTrajectory = drawTrajectory||singleLineTrajectory;
drawBodyCurve = drawBodyCurve||singleLineTrajectory;
doExport = doExport||singleLineTrajectory;
clearFigure = (clearFigure && doExport)||singleLineTrajectory;
if ~(drawBodyCurve || drawTrajectory)
    return
end
maxI = size(config, 1);
framenumber = round(maxI/1000*min([1000 dataResolution]));
indeces = round(linspace(1, maxI, framenumber));
f1 = figure('InnerPosition', [0 0 1280 720]);
colormap([0.3 0.3 0.3])
lightangle(40,15)
xlabel('x');ylabel('y');zlabel('z');
daspect([1 1 1]);
xlim([-60 60]);ylim([-60 60]);zlim([0 200]);
view(40,25);
hold on
if  drawTrajectory
    if ~drawBodyCurve
        xlim([-100 100]);ylim([-100 100]);zlim([140 180]);
        view(60,20)
    end
    headData = zeros(3,framenumber);
end
if drawBodyCurve
    bodyData = cell(1,framenumber);
end
i2 = 1;
tic
for i = indeces
    x = config(i,:);
    if drawTrajectory && ~drawBodyCurve
        headData(:,i2) = line2(x(1), x(2), x(3), x(4), x(5),0);
    elseif ~drawTrajectory && drawBodyCurve
        [~,bodyData{:,i2}] = line2(x(1), x(2), x(3), x(4), x(5),1);
    elseif drawBodyCurve && drawTrajectory
        [headData(:,i2),bodyData{:,i2}] = line2(x(1), x(2), x(3), x(4), x(5),1);
    end
    i2 = i2 + 1;
end
toc
if doExport
    %     v = VideoWriter('out3', 'MPEG-4'); %only use when using win7+ or MacOSX 10.7+
    v = VideoWriter('out3', 'Motion JPEG AVI'); % works on linux
    v.FrameRate = framerate;
    v.Quality = 100;
    open(v);
    meanLoopTime = 0;
    for i = 1:framenumber
        loopStartTime = tic;
        if clearFigure
            if drawTrajectory
                if singleLineTrajectory
                    plot3(headData(1,1:i),headData(2,1:i),headData(3,1:i),'-r');
                else
                    i2 = max(1, i-1);
                    plot3(headData(1,i2:i),headData(2,i2:i),headData(3,i2:i),'-r');
                end
            end
            if drawBodyCurve
                x = bodyData{:,i};
                h = streamtube({[x(1,:);x(2,:);x(3,:)]'}, 3);
                h.FaceLighting = 'gouraud';
                h.AmbientStrength = 0.8;
                h.DiffuseStrength = 0.7;
                h.SpecularStrength = 0.8;
                h.SpecularExponent = 15;
                h.BackFaceLighting = 'unlit';
                lightangle(40,15)
                shading interp
            end
            writeVideo(v, getframe(f1));
            cla(f1);%replace with delete() twice
        else
            if drawTrajectory
                j = plot3(headData(1,1:i),headData(2,1:i),headData(3,1:i),'-r');
            end
            if drawBodyCurve
                x = bodyData{:,i};
                h = streamtube({[x(1,:);x(2,:);x(3,:)]'}, 3);
                h.FaceLighting = 'gouraud';
                h.AmbientStrength = 0.8;
                h.DiffuseStrength = 0.7;
                h.SpecularStrength = 0.8;
                h.SpecularExponent = 15;
                h.BackFaceLighting = 'unlit';
                shading interp
            end
            writeVideo(v, getframe(f1));
            if drawTrajectory
                delete(j)
            end
        end
        loopEndTime = toc(loopStartTime);
        meanLoopTime = (meanLoopTime * (i-1) + loopEndTime)/(i);
        rTime = round(meanLoopTime*(framenumber-i),2);
        fprintf('%d of %d Frames (%.2f%%); estimated time remaining: %.0fh %.0fm %.0fs\n',...
            i, framenumber, round(i/framenumber*100,2), floor(rTime/3600),floor(mod(rTime/60,60)), floor(mod(rTime,60)));
    end
else
    tic
    if drawTrajectory
        plot3(headData(1,:),headData(2,:),headData(3,:),'-r');
    end
    if drawBodyCurve
        for i = 1:framenumber
            x = bodyData{:,i};
            h = streamtube({[x(1,:);x(2,:);x(3,:)]'}, 3);
            h.FaceLighting = 'gouraud';
            h.AmbientStrength = 0.8;
            h.DiffuseStrength = 0.7;
            h.SpecularStrength = 0.8;
            h.SpecularExponent = 15;
            h.BackFaceLighting = 'unlit';
        end
        shading interp
    end
    toc
end
% savefig(f1, 'fig_script_compact.fig', 'compact');
rTime = toc(startTime);
fprintf('Time elapsed: %.0fh %.0fm %.0fs %.0fms\n', floor(rTime/3600),floor(mod(rTime/60,60)), floor(mod(rTime,60)), (rTime-floor(rTime))*1000);