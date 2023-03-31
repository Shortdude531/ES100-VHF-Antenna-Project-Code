% Take a thermal RGB image from the FLIR One camera and uses the embedded color bar to determine temperatures from the colors and make a temperature image.
clc;    % Clear the command window.
close all;  % Close all figures (except those of imtool.)
clear;  % Erase all existing variables. Or clearvars if you want.
workspace;  % Make sure the workspace panel is showing.
format long g;
format compact;
fontSize = 15;
%===============================================================================
% Get the name of the image the user wants to use.
baseFileName = 'Lamp.png'; % Base file name with no folder prepended (yet).
% Get the full filename, with path prepended.
folder = pwd; % Change to whatever folder the image lives in.
fullFileName = fullfile(folder, baseFileName);  % Append base filename to folder to get the full file name.
if ~isfile(fullFileName)
	errorMessage = sprintf('Error: file not found:\n%s', fullFileName)
	uiwait(errordlg(errorMessage));
	return;
end
fprintf('Transforming image "%s" to a thermal image.\n', fullFileName);
%===============================================================================
% Read in a demo image.
rgbImage = imread(fullFileName);
% Display the image.
subplot(2, 3, 1);
h1 = imshow(rgbImage, []);
axis on;
caption = sprintf('Original Pseudocolor Image, %s', baseFileName);
title(caption, 'FontSize', fontSize, 'Interpreter', 'None');
xlabel('Column', 'FontSize', fontSize, 'Interpreter', 'None');
ylabel('Row', 'FontSize', fontSize, 'Interpreter', 'None');
drawnow;
%grayImage = min(rgbImage, [], 3); % Useful for finding image and color map regions of image.
%=========================================================================================================
% Need to read in the color bar separately.
% Get the full filename, with path prepended.
baseFileName = 'Screen Shot 2023-03-15 at 12.44.48 PM.png'; % Base file name with no folder prepended (yet).
folder = pwd; % Change to whatever folder the image lives in.
fullFileName = fullfile(folder, baseFileName);  % Append base filename to folder to get the full file name.
if ~isfile(fullFileName)
	errorMessage = sprintf('Error: file not found:\n%s', fullFileName)
	uiwait(errordlg(errorMessage));
	return;
end
colorBarImage = imread(fullFileName);
 %b = colorBarImage(:,:,3);
%=========================================================================================================
% Display the colorbar image.
subplot(2, 3, 2);
h2 = imshow(colorBarImage, []);
axis on;
caption = sprintf('Cropped Colorbar Image');
title(caption, 'FontSize', fontSize, 'Interpreter', 'None');
xlabel('Column', 'FontSize', fontSize, 'Interpreter', 'None');
ylabel('Row', 'FontSize', fontSize, 'Interpreter', 'None');
hp = impixelinfo();
drawnow;
% Set up figure properties:
% Enlarge figure to full screen.
g = gcf;
g.WindowState = 'maximized';
% Get rid of tool bar and pulldown menus that are along top of figure.
% set(gcf, 'Toolbar', 'none', 'Menu', 'none');
% Give a name to the title bar.
set(gcf, 'Name', 'Demo by ImageAnalyst', 'NumberTitle', 'Off')
%=========================================================================================================
% Get the color map from the color bar image.
storedColorMap = colorBarImage(5:445, 19,:);
% Need to call squeeze to get it from a 3D matrix to a 2-D matrix.
% Also need to divide by 255 since colormap values must be between 0 and 1.
storedColorMap = double(squeeze(storedColorMap)) / 255;
% Need to flip up/down because the low rows are the high temperatures, not the low temperatures.
storedColorMap = flipud(storedColorMap);
subplot(2, 3, 3);
plot(storedColorMap(:, 1), 'r-', 'LineWidth', 3);
hold on;
grid on;
plot(storedColorMap(:, 2), 'g-', 'LineWidth', 3);
plot(storedColorMap(:, 3), 'b-', 'LineWidth', 3);
title('Stored Color Map', 'FontSize', fontSize);
% Convert the subject/sample from a pseudocolored RGB image to a grayscale, indexed image.
indexedImage = rgb2ind(rgbImage, storedColorMap);
% Display the indexed image.
subplot(2, 3, 4);
h4 = imshow(indexedImage, []);
axis on;
caption = sprintf('Indexed Image (Gray Scale Thermal Image)');
title(caption, 'FontSize', fontSize, 'Interpreter', 'None');
xlabel('Column', 'FontSize', fontSize, 'Interpreter', 'None');
ylabel('Row', 'FontSize', fontSize, 'Interpreter', 'None');
hp = impixelinfo();
drawnow;
%=========================================================================================================
% Now we need to define the temperatures at the end of the colored temperature scale.
% You can read these off of the image, since we can't figure them out without doing OCR on the image.
% Define the temperature at the top end of the scale.
% This will probably be the high temperature.
highTemp = 432;
% Define the temperature at the dark end of the scale
% This will probably be the low temperature.
lowTemp = 24.2;
% Scale the indexed gray scale image so that it's actual temperatures in degrees C instead of in gray scale indexes.
thermalImage = lowTemp + (highTemp - lowTemp) * mat2gray(indexedImage);
% Display the thermal image.
subplot(2, 3, 5);
h5 = imshow(thermalImage, []);
axis on;
colorbar;
hp = impixelinfo(gcf, [h1, h2, h4, h5]);
title('Floating Point Thermal (Temperature) Image', 'FontSize', fontSize, 'Interpreter', 'None');
xlabel('Column', 'FontSize', fontSize, 'Interpreter', 'None');
ylabel('Row', 'FontSize', fontSize, 'Interpreter', 'None');
% Let user mouse around and see temperatures on the GUI under the temperature image.
hp = impixelinfo();
hp.Units = 'normalized';
hp.Position = [0.45, 0.03, 0.25, 0.05];
%=========================================================================================================
% Get and display the histogram of the thermal image.
subplot(2, 3, 6);
histogram(thermalImage, 'Normalization', 'probability');
axis on;
grid on;
caption = sprintf('Histogram of Thermal Image');
title(caption, 'FontSize', fontSize, 'Interpreter', 'None');
xlabel('Temperature [Degrees]', 'FontSize', fontSize, 'Interpreter', 'None');
ylabel('Frequency [Pixel Count]', 'FontSize', fontSize, 'Interpreter', 'None');
% Get the maximum temperature.
maxTemperature = max(thermalImage(:));
fprintf('The maximum temperature in the image is %.2f\n', maxTemperature);
fprintf('Done!  Thanks Image Analyst!\n');

%% ======================================================
 indices = find(abs(thermalImage)<50);
avgthermal = thermalImage;
fprintf("Temp over area k")
avgthermal(indices) = [];

avgheat = sum(avgthermal+273)/length(avgthermal)
fprintf("T4 k^4")
avgheatt4 = (avgthermal+273).^4;
s1  = sum(avgheatt4)/length(avgheatt4)

fprintf("T4 to the quater k^4") 

s1  = (s1)^(1/4)
set(findall(gcf,'-property','FontSize'),'FontSize',20)