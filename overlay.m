%{
Timer Object: A timer object is used to control the frame updates, ensuring they occur at a rate matching the video's frame rate.
This script provides a basic framework for video playback with play/pause control. 
It should play the video at its natural frame rate, provided the system can handle the processing load. If the video playback is not smooth, it could be due to the limitations of the system's processing power or the MATLAB environment's ability to render video frames in real-time.
%}

function overlay()
    % Load Video and Data
    videoFile = '/Users/haozhema/sign-language-dissertation-tcd-2023/stimuli/V01.mp4';
    vidObj = VideoReader(videoFile);
    %data = load('/Users/haozhema/sign-language-dissertation-tcd-2023/envelope/V01.mat');
    data = load('/Users/haozhema/sign-language-dissertation-tcd-2023/IVC_openpose_envelope/ivc_values_final.mat');
    %frameDifferences = data.feature;
    frameDifferences = data.ivc_values;

    % Create GUI
    hFig = figure;
    hAx1 = subplot(2,1,1);
    hImg = imshow(readFrame(vidObj), 'Parent', hAx1);
    title('Video Playback');
    hAx2 = subplot(2,1,2);
    hPlot = plot(hAx2, frameDifferences(1), 'LineWidth', 2);
    xlim([1 length(frameDifferences)]);
    ylim([min(frameDifferences) max(frameDifferences)]);
    title('Frame Differences');
    hButton = uicontrol('String', 'Play', 'Callback', @togglePlayPause);

    % Initialize Timer
    timerObj = timer('TimerFcn', @updateFrame, 'Period', 1/vidObj.FrameRate, 'ExecutionMode', 'fixedRate');

    % Play/Pause Toggle
    isPlaying = false;
    function togglePlayPause(~, ~)
        if ~isPlaying
            start(timerObj);
            set(hButton, 'String', 'Pause');
        else
            stop(timerObj);
            set(hButton, 'String', 'Play');
        end
        isPlaying = ~isPlaying;
    end

    % Update Frame Function
    function updateFrame(~,~)
        if hasFrame(vidObj) && ishandle(hImg)
            videoFrame = readFrame(vidObj);
            set(hImg, 'CData', videoFrame);
            currentFrame = round(vidObj.CurrentTime * vidObj.FrameRate);
            if currentFrame <= length(frameDifferences)
                set(hPlot, 'YData', frameDifferences(1:currentFrame));
            end
            drawnow;
        else
            stop(timerObj);
            set(hButton, 'String', 'Play');
            isPlaying = false;
        end
    end

    % Clean up on figure close
    set(hFig, 'CloseRequestFcn', @closeFigure);
    function closeFigure(~, ~)
        stop(timerObj);
        delete(timerObj);
        delete(hFig);
    end
end
