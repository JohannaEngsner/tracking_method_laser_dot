clc
clear all
close all

% Remove any old messageboxes
oldMessageBoxes = findall(0, 'Type', 'Figure');
delete(oldMessageBoxes( arrayfun(@(h) contains(h.Tag, 'Msgbox'), oldMessageBoxes)))

% Making windows fullscreen
set(groot, 'defaultFigureWindowState', 'maximized');

%  Instructions for choosing a video
waitfor(msgbox('Start by choosing the video you want to analyze.', 'Choosing video'));
pause(0.2);

% Choosing a video
[videoFileName,videoFilePath] = videoChooser();

%  Instructions for how to crop the video
waitfor(msgbox('Please crop the video so only the paper is shown. Be careful not to accidentally crop the figure. After choosing the correct area, double click to crop.', 'Cropping video'));
pause(0.2);

% Letting the user crop the video
video=VideoReader(videoFileName);
writerObj = VideoWriter('CroppedVideo.avi');
open(writerObj);
firstFrame=readFrame(video);
[croppedImage, rect] = imcrop(firstFrame);

% Creating a processing bar
h = waitbar(0, 'Cropping video...', 'Name', 'Video cropper');
n = video.NumFrames;

% Cropping all frames
for i=1:n
  waitbar(i / n, h);
  im=read(video,i);
  imc=imcrop(im,rect); % Cropping each frame
  scaleFactor = 480 / size(imc, 1); 
  imr=imresize(imc,scaleFactor); % Resizing based on the previous proportions
  writeVideo(writerObj,imr);     
end

pause(0.2);
delete(h) % Removing the processing bar
close(writerObj)

close all;

% Reading Video
vidObj = VideoReader('CroppedVideo.avi');
vidObj.CurrentTime = 0;

% Reading & resizing first frame to make the calculations quicker
vidFrame=readFrame(vidObj);

% Choosing corners
[x, y] = cornerChooser(vidFrame);

close all;

% Creating the different zones based on the chosen corners
[binaryThinLine, binaryThickLine, binaryBlueArea, binaryYellowArea, binaryRedArea] = lineDrawer(x, y, vidFrame);

pause(0.2);
set(groot, 'defaultFigureWindowState', 'maximized');
set(groot,'defaultFigureVisible','on');

% Creating a matrix and making sure the video starts from the beginning
allCenters = [];
vidObj.CurrentTime = 0;
imshow(vidFrame); % To make sure the processing bar overlays the video

% Creating a processing bar
h = waitbar(0, 'Calculating...', 'Name', 'Laser tracker');
a=1;

while(vidObj.hasFrame)
    vidFrame = readFrame(vidObj);
    vidframeGreen= vidFrame(:,:,2);
    binaryFrame = imbinarize(vidframeGreen,0.93); % Making sure only the laser dot is showing

    % Removing unwanted noise from the binary frame
    seSize = 3; % Choose this depending on how far away from the wall the patient is sitting. (If 1m away, choose seSize=1)
    se = strel('disk', seSize);
    binaryFrame=imopen(binaryFrame,se);
  
    stats = regionprops("table",binaryFrame,"Centroid", "MajorAxisLength","MinorAxisLength");

    [~,indexMax]=max(stats.MajorAxisLength);
    centers = stats.Centroid(indexMax,:); % Choosing the largest dot 
     if isempty(stats) == 0    % Checking if the laser dot was found to avoid errors
        allCenters(a,1) = centers(1,1);
        allCenters(a,2) = centers(1,2);
        a=a+1;
     end

    hold on % Visualising the tracking of the laser dot 
    imshow(vidFrame)
    hold on
    viscircles(centers,5,'Color', '#e697c6');
    waitbar(a / n, h);

    pause(1/vidObj.FrameRate)
end

pause(0.2);
delete(h) % Deleting the processing bar

[results] = zoneCalculator(allCenters,binaryThinLine, binaryThickLine, binaryBlueArea, binaryYellowArea, binaryRedArea,vidObj.NumFrames);

% Calculating the total time
totalTime = vidObj.NumFrames / vidObj.FrameRate;

% Creating the result message 
message = sprintf('Thin line: % .2f % % \n', results(1)*100);
message = [message, sprintf('Thick line: % .2f % % \n', results(2)*100)];
message = [message, sprintf('Blue zone: % .2f % % \n', results(3)*100)];
message = [message, sprintf('Yellow zone: % .2f % % \n', results(4)*100)];
message = [message, sprintf('Red zone: % .2f % % \n', results(5)*100)];
message = [message, sprintf('Outside: % .2f % % \n', results(6)*100)];
message = [message, sprintf('Total time: % .2f seconds\n', totalTime)];

% Displaying the results
msgbox(message, 'Results');
%% Function videoChooser
function [videoFileName,videoFilePath] = videoChooser()

[videoFileName,videoFilePath] = uigetfile('*.*', 'Select a file');
if(videoFilePath==0)
    return;
end
end
%% Function cornerChooser
function[x_corners, y_corners] = cornerChooser(frame)

imshow(frame);
waitfor(msgbox('Click on the four corners of the black line in the figure. Please click on the outmost part of the corners, as far out as possible.', 'Click the corners of the black line'));

x_corners = zeros(4, 1);
y_corners = zeros(4, 1);

dotHandles = cell(4, 1);
labels = cell(4, 1);

% Prompt the user to click on each corner and display labels
cornerNames = {'Upper Left', 'Upper Right', 'Lower Left', 'Lower Right'};

for i = 1:4
    waitfor(msgbox(['Click on the ', cornerNames{i}, ' corner.']))
    [x, y] = ginput(1);
    
    % Store corner coordinates
    x_corners(i) = x;
    y_corners(i) = y;

    % Adjust label position based on corner index
    label_x = x;
    if i == 1 || i == 3
        screenSize = get(0, 'ScreenSize');
        label_x = label_x - screenSize(3)/20; % Move left for corners 1 and 3

    elseif i == 2 || i == 4
        label_x = label_x + screenSize(3)/65; % Move right for corners 2 and 4
    end
    
    % Plot the dot and label for the current corner
    hold on;
    dotHandles{i} = plot(x, y, '-o','Color','#e697c6', 'MarkerSize', 8);
    labels{i} = text(label_x, y, cornerNames{i},'BackgroundColor', '#e697c6', 'FontSize', 12);
    hold off;
end

%  Allow the user to adjust corners
while true
    % Ask if the user wants to change the corners
    choice = questdlg('Do you want to change the position of any of the corners?', ...
                      'Corner Adjustment', 'Yes', 'No', 'No');
    
    if strcmp(choice, 'Yes')
        % Prompt the user to click on a corner to adjust
        waitfor(msgbox('Click on the corner you want to adjust.'));
        [x_adjust, y_adjust, button] = ginput(1);

        % Check if the user clicked on a corner
        if button == 1
            % Find the closest corner to the clicked position
            distances = sqrt((x_corners - x_adjust).^2 + (y_corners - y_adjust).^2);
            [~, closestCorner] = min(distances);

            % Remove the dot and label of the selected corner
            delete(dotHandles{closestCorner});
            delete(labels{closestCorner});

            % Prompt the user to click on the new position for the corner
            waitfor(msgbox(['Click where you want to place the ', cornerNames{closestCorner}, ' corner.']));
            [x_new, y_new] = ginput(1);

            % Update corner coordinates
            x_corners(closestCorner) = x_new;
            y_corners(closestCorner) = y_new;

            % Adjust the x-coordinate of the label for the selected corner
            if closestCorner == 1 || closestCorner == 3
                labelOffset = -screenSize(3)/20; %  Move left for corners 1 and 3
            elseif closestCorner == 2 || closestCorner == 4
                labelOffset = screenSize(3)/65; %  Move right for corners 2 and 4
            end
 
            % Plot the updated dot and label for the selected corner
            hold on;
            dotHandles{closestCorner} = plot(x_new, y_new, '-o','Color','#e697c6', 'MarkerSize', 8);
            labels{closestCorner} = text(x_new + labelOffset, y_new, cornerNames{closestCorner}, 'BackgroundColor', '#e697c6', 'FontSize', 12);            hold off;
        end

    elseif strcmp(choice, 'No') || isempty(choice)
        % If the user chooses not to change corners or closes the dialog, exit the loop
        break;
    end
end
end

%% Function zoneCalculator
function [results] = zoneCalculator(allCenters, binaryThinLine, binaryThickLine, binaryBlueArea, binaryYellowArea, binaryRedArea,numFrames)
thinFrames=0;
thickFrames=0;
blueFrames=0;
yellowFrames=0;
redFrames=0;

% Going through all centroid coordinates and checking which zone they match
for i = 1:size(allCenters,1)
    onThin=0; % Keeping track of whether a surrounding pixel is on the thinLine
   for a = 1:9 
        switch a
            case 1
                offset_x = 0;
                offset_y = 0;
            case 2
                offset_x = 1;
                offset_y = 0;
            case 3
                offset_x = -1;
                offset_y = 0;
            case 4
                offset_x = 0;
                offset_y = 1;
            case 5
                offset_x = 0;
                offset_y = -1;
            case 6
                offset_x = 1;
                offset_y = -1;
            case 7
                offset_x = 1;
                offset_y = 1;
            case 8
                offset_x = -1;
                offset_y = -1;
            case 9
                offset_x = -1;
                offset_y = 1;
        end

if allCenters(i,2) >1 && allCenters(i,2)<479 && allCenters(i,1)>1 && allCenters(i,1)<width(binaryThinLine)-1
if binaryThinLine(int64(allCenters(i,2)) + offset_y, int64(allCenters(i,1)) + offset_x)==0
        thinFrames=thinFrames+1; % Updating counter if the coordinates match the zone
        onThin=1;
        break;
        
end
end
    end

if onThin==0    
if binaryThickLine(int64(allCenters(i,2)),int64(allCenters(i,1)))==0
    thickFrames=thickFrames+1; % Updating counter if the coordinates match the zone

elseif binaryBlueArea(int64(allCenters(i,2)),int64(allCenters(i,1)))==0
    blueFrames=blueFrames+1; % Updating counter if the coordinates match the zone

elseif binaryYellowArea(int64(allCenters(i,2)),int64(allCenters(i,1)))==0
    yellowFrames=yellowFrames+1; % Updating counter if the coordinates match the zone

elseif binaryRedArea(int64(allCenters(i,2)),int64(allCenters(i,1)))==0
    redFrames=redFrames+1; % Updating counter if the coordinates match the zone

end
end
end

% Calculating the results
results(1)=thinFrames/numFrames;
results(2)=thickFrames/numFrames;
results(3)=blueFrames/numFrames;
results(4)=yellowFrames/numFrames;
results(5)=redFrames/numFrames;

% All centroids that weren't located in a zone belong to the outside
untracked=numFrames-size(allCenters,1);
results(6)=((size(allCenters,1)-thinFrames-thickFrames-blueFrames-yellowFrames-redFrames)+untracked)/numFrames;

end

%% LineDrawer 
function[binaryThinLine, binaryThickLine, binaryBlueArea, binaryYellowArea, binaryRedArea] = lineDrawer(x_corners, y_corners, frame)
set(groot,'defaultFigureWindowState','normal')
set(groot, 'defaultFigureUnits', 'pixels');
set(groot,'defaultFigureVisible','off');

segment = (x_corners(2)-x_corners(1))*0.048; % Relationship between a segment and the length of thin line. 
% Segment is used to calculate the thickness of each zone.
factor=0.090; % The thickness of green zone compared to its length


thickness(1) = (x_corners(2)-x_corners(1))*0.00427; % Thickness of thin line
thickness(2) = (x_corners(2)-x_corners(1))*(factor/2); % Green  
thickness(3) = thickness(2)+segment; % Blue
thickness(4) = thickness(2)+1.9*segment; % Yellow
thickness(5) = thickness(2)+2.8*segment; % Red

binaryThinLine=drawLine(x_corners,y_corners,thickness(1),frame,0,0);
binaryThickLine=drawLine(x_corners,y_corners,thickness(2),frame,0.8,2);
binaryBlueArea=drawLine(x_corners,y_corners,thickness(3),frame,0.86,1.95);
binaryYellowArea=drawLine(x_corners,y_corners,thickness(4),frame,0.92,1.9);
binaryRedArea=drawLine(x_corners,y_corners,thickness(5),frame,0.98,1.9);

end

%% Function drawLine
function[binaryIm]=drawLine(x_corners,y_corners,thickness,frame,factor,factor2)
whiteImage = ones(size(frame), 'uint8') * 255; % White image used as background for the zones

color = [0 0 1]; % Choosing blue as the color for the plots

% Drawing the zones using lines
whiteImage = insertShape(whiteImage, 'Line', [x_corners(1), y_corners(1), x_corners(2), y_corners(2)], 'Color', color, 'Opacity', 0.25, 'LineWidth', int64(thickness), 'SmoothEdges', true);
whiteImage = insertShape(whiteImage, 'Line', [x_corners(2), y_corners(2), x_corners(3), y_corners(3)], 'Color', color, 'Opacity', 0.25, 'LineWidth', int64(thickness), 'SmoothEdges', true);
whiteImage = insertShape(whiteImage, 'Line', [x_corners(1), y_corners(1), x_corners(4), y_corners(4)], 'Color', color, 'Opacity', 0.25, 'LineWidth', int64(thickness), 'SmoothEdges', true);
whiteImage = insertShape(whiteImage, 'Line', [x_corners(3), y_corners(3), x_corners(4), y_corners(4)], 'Color', color, 'Opacity', 0.25, 'LineWidth', int64(thickness), 'SmoothEdges', true);

% Adding triangles to the corners so they will match the pattern
if factor ~=0
cornerTriangles = [x_corners(1)-(thickness*factor2), y_corners(1)-thickness/2, (x_corners(2)+x_corners(1))/2-(factor*thickness), (y_corners(3) +y_corners(1))/2+(factor2), x_corners(1) + 0.7*thickness, y_corners(1)-0.5*thickness;
                    x_corners(2)+(thickness*factor2), y_corners(2)-thickness/2, (x_corners(2)+x_corners(1))/2+(factor*thickness), (y_corners(3) +y_corners(1))/2+(factor2), x_corners(2) - 0.7*thickness, y_corners(2)-0.5*thickness;
                    x_corners(3)-(thickness*factor2), y_corners(3)+thickness/2, (x_corners(4)+x_corners(3))/2-(factor*thickness), (y_corners(3) +y_corners(1))/2-(factor2), x_corners(3) + 0.7*thickness, y_corners(3)+0.5*thickness;
                    x_corners(4)+(thickness*factor2), y_corners(4)+thickness/2, (x_corners(4)+x_corners(3))/2+(factor*thickness), (y_corners(3) +y_corners(1))/2-(factor2), x_corners(4) - 0.7*thickness, y_corners(4)+0.5*thickness];

for i = 1:size(cornerTriangles, 1)
    whiteImage = insertShape(whiteImage, 'FilledPolygon', cornerTriangles(i,:), 'Color', color, 'Opacity', 0.25);
end
end

grayScaleImage = whiteImage(:,:,1); % Choosing the blue RGB channel
binaryIm = imbinarize(grayScaleImage, 0.9); % Turning it into a binary image
end