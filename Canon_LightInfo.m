% Canon lux calculation, mask method test
% Hankun Li, University of Kansas
% 04/17/2021
clear all; %#ok<*CLALL>
close all;
clc;
% camera parameters
sy = 14.9; sx = 22.3; f = 4.5; Lfov = 180;

%% Load luminance map
yn = yn_dialog("Get Luminance map from HDR image [YES]? or Load a txt file [NO]?");
if ismember(yn, ['Yes', 'yes'])
    [hdr,lmap,luxmask,sAmask] = preProcess_hdr(f,sy,Lfov);
else
    [lmap,luxmask,sAmask] = preProcess_lmap(f,sy,Lfov);
end

%% display HDR image, optional, ONLY works for HDR input
if hdr
    [hdr_converted] = hdrGammaShow(hdr,0.5); % defautl gamma = 0.5, can be adjusted.
    imshow(hdr_converted);
end

%% visulize luminance map
Luminance_map_show(lmap,'Luminance map');

%% calculate illuminance
lx = sum(sum(luxmask.*lmap));
fprintf('Calculated illuminance is: %.2f lx\n',lx);

%% retrieve luminance, solid angle, illuminance contribution of intersted part
[luminance,shape] = LuminanceRetrieveROI(lmap);
fprintf('[1] Average luminance of selection: %.2f cd/m2 \n',luminance);
fprintf('[2] total solid angle of selection: %.2f sr \n',sum(sum(shape.*sAmask)));
fprintf('[3] illuminance contribution of selection: %.4f percent \n',...
    sum(sum(luxmask.*lmap.*shape))/lx);

% main function end...








%% Functions library
function [lmap,luxmask,sAmask] = preProcess_lmap(f,sy,Lfov)
lmap = circularCrop_img(loadLuminanceMap(),f,sy);
[luxmask, sAmask] = FASTequisolid_lmapVer(lmap,Lfov);
end

function [luxmask,sAmask] = FASTequisolid_lmapVer(lmap,fov)
[yc,xc] = size(lmap); [xf,yf] = meshgrid(1:xc,1:yc); 
xf = (xf - round(xc/2))./round(xc/2);
yf = (yf - round(yc/2))./round(yc/2);
phiS = 2*asind(sqrt(yf.^2+xf.^2)*sind(fov/4));
bound = imbinarize(lmap,0); 
sA = 2*pi/sum(sum(bound));
luxmask = cosd(phiS).*bound.*sA;
if (3.12 <= round(sum(sum(luxmask)),2))&& (round(sum(sum(luxmask)),2)<= 3.16)
    fprintf('solid-angle calculation...OK...error < 0.05 percent\n')
else
    warndlg("solid-angle calculation...error > 0.05 percent\n");
end
luxmask = luxmask.*(pi/sum(sum(luxmask)));
sAmask = sA.*double(bound);
end

function [hdr,lmap,luxmask,sAmask] = preProcess_hdr(f,sy,Lfov)
[hdr, cf] = loadHDR(); hdr = hdr_resize(hdr);
[y,x] = size(hdr(:,:,1));
pr180 = round(sind(45)*2*f/sy*min([x,y])/2)*2;
hdr = circularCrop_hdr(hdr,x,y,pr180);
[luxmask,vcmask,sAmask] = FASTequisolid(hdr,Lfov);
lmap = luminanceMap_gen(hdr,cf).*vcmask;
end

function [luminance,shape] = LuminanceRetrieveROI(lmap)
gm = 0.5; agm = autoGamma(lmap);
fprintf('default gamma= %.2f \n',gm);
fprintf('auto gamma= %.2f \n',agm);
yn = yn_dialog('using auto-calculated gamma? [check Command Window info]');
if ismember(yn, ['Yes', 'yes'])
    gm = agm;
end
%%
lmap(lmap<0) = 0;Limg = (lmap - min(min(lmap)))/(max(max(lmap))-min(min(lmap)));
Limg = uint8((Limg.^gm).*256);
%%
uiwait(msgbox({'select target surface.'},'Notice!'));
figure(1); shape = roipoly(Limg); close(gcf);
if isempty(shape)
    uiwait(msgbox('No ROI selected!','Error','error'));
    error('NO ROI selected\n..');
end
lmap_roi = lmap.*shape; luminance = mean(lmap_roi(lmap_roi>0));
end

% VC correction, require calibration for every Camera + Lens
function [vcf] = canonVC(angle) % vc correction for Canon t2i + sigma f2.8
vcf = 1./(-4.3909e-09.*angle.^(4) - 3.9024e-07.*angle.^(3) +...
    3.3680e-05.*angle.^(2)-0.0018.*angle + 1.0018);
end

function [luxmask,vcmask,sAmask] = FASTequisolid(hdr,fov)
[yc,xc] = size(hdr(:,:,1));
[xf,yf] = meshgrid(1:xc,1:yc); 
xf = (xf - round(xc/2))./round(xc/2);
yf = (yf - round(yc/2))./round(yc/2);
phiS = 2*asind(sqrt(yf.^2+xf.^2)*sind(fov/4));
vcmask = canonVC(phiS); vcmask(vcmask<0) = 0;
bound = imbinarize(rgb2gray(hdr),0); 
vcmask = vcmask.*bound; 
sA = 2*pi/sum(sum(bound));
luxmask = cosd(phiS).*bound.*sA;
if (3.12 <= round(sum(sum(luxmask)),2))&& (round(sum(sum(luxmask)),2)<= 3.16)
    fprintf('solid-angle calculation...OK...error < 0.05 percent\n')
else
    maskShow(luxmask);
    fprintf('solid-angle calculation...error > 0.05 percent\n')
end
luxmask = luxmask.*(pi/sum(sum(luxmask)));
sAmask = sA.*double(bound);
end

function maskShow(luxmask)
immsk = uint8((luxmask - min(min(luxmask)))./(max(max(luxmask))-min(min(luxmask))).*255);
f = figure(1); imshow(immsk);
uiwait(f);
end

function lmap = luminanceMap_gen(hdr,cf)
lmap = (hdr(:,:,1).*0.265 + hdr(:,:,2).*0.670 + hdr(:,:,3).*0.065).*179.*cf;
end

function Luminance_map_show(lmap,name)
cv = std(std(lmap))/mean(mean(lmap));
lmap(lmap<0) = 0;lumimg = (lmap - min(min(lmap)))/(max(max(lmap))-min(min(lmap)));
if  (1.5<cv)&&(cv<10)
    gm = round(1/cv,2);
elseif cv>10
    gm = 0.09;
else
    gm = 1;
end
lumimg = uint8((lumimg.^gm).*256);
rg = max(max(lmap))-min(min(lmap)); crange = jet(256);crange(1,:) = 0;
cb1 = round(rg.*(0.03316.^(1/gm)),7);cb2 = round(rg.*(0.26754.^(1/gm)),2);
cb3 = round(rg.*(0.50191.^(1/gm)),2);cb4 = round(rg.*(0.73629.^(1/gm)),2);
cb5 = round(rg.*(0.97066.^(1/gm)),2);
figure(2);imshow(lumimg,'Colormap',crange); title(name,'FontSize',20);
hcb = colorbar('Ticks',[8,68,128,188,248],'TickLabels',{cb1,cb2,cb3,cb4,cb5});
hcb.FontSize = 18; title(hcb,'luminance (cd/m2)');
end

function fig1 = hdrGammaShow(imHDR,gamma)
fr = imHDR(:,:,1);fg = imHDR(:,:,2);fb = imHDR(:,:,3);
fr = single(fr).^gamma;fg = single(fg).^gamma;fb = single(fb).^gamma;
fig1 = cat(3,fr,fg,fb);
end

function exposure = getExpValue(filename)
fid = fopen(filename); ct = 0;
while ct < 16
    line = fgetl(fid);
    if contains(line, 'EXPOSURE')
        line = erase(line, ' ');
        break
    end
end
fclose(fid); exposure = str2double(erase(line, 'EXPOSURE='));
end

function yn = yn_dialog(ques)
opts.Interpreter = 'tex'; opts.Default = 'No';
yn = questdlg(ques,'Dialog Window',...
    'Yes','No',opts);
end

function Icropped = circularCrop_hdr(I,x,y,r)
xc = round(x/2);yc = round(y/2);
c = zeros(y,x); [L(:,1),L(:,2)] = find(c == 0);
L(:,3) = sqrt((L(:,1) - yc).^2 + (L(:,2) - xc).^2);
L(L(:, 3) > r, :) = [];
for i = 1: size(L,1)
   c(y+1-L(i,1),L(i,2)) = 1;
end
msk = imbinarize(c,0);
ir = double(I(:,:,1)).*msk;
ig = double(I(:,:,2)).*msk;
ib = double(I(:,:,3)).*msk;
Icc = cat(3,ir,ig,ib);
[mski(:,1), mski(:,2)] = find(msk == 1);
Icropped = imcrop(Icc,[min(mski(:,2)),min(mski(:,1)),...
    max(mski(:,2))-min(mski(:,2)),max(mski(:,1))-min(mski(:,1))]);
end

function Icropped = circularCrop_img(I,f,sy)
[y,x] = size(I); r = round(sind(45)*2*f/sy*min([x,y])/2)*2;
xc = round(x/2);yc = round(y/2);
c = zeros(y,x); [L(:,1),L(:,2)] = find(c == 0);
L(:,3) = sqrt((L(:,1) - yc).^2 + (L(:,2) - xc).^2);
L(L(:, 3) > r, :) = [];
for i = 1: size(L,1)
   c(y+1-L(i,1),L(i,2)) = 1;
end
msk = imbinarize(c,0); Icc = double(I).*msk;
[mski(:,1), mski(:,2)] = find(msk == 1);
Icropped = imcrop(Icc,[min(mski(:,2)),min(mski(:,1)),...
    max(mski(:,2))-min(mski(:,2)),max(mski(:,1))-min(mski(:,1))]);
end


function gamma = autoGamma(lmap)
cv = std(std(lmap))/mean(mean(lmap));
if  (1.1<=cv)&&(cv<=10)
    gamma = round(1/cv,2);
elseif cv > 10
    gamma = 0.09;
else
    gamma = 1;
end
end

function [ans1] = ui_cf()
prompt = {'Global luminance calibration factor? [use 1 in the calibration]'};
dlgtitle = 'User Input'; dims = [1 50];definput = {'1.0'};
answer = str2double(inputdlg(prompt,dlgtitle,dims,definput));
if isempty(answer)
    ans1 = 1.0;
else
    ans1 = answer(1);
end
end

function [newHDR] = hdr_resize(HDR)
prompt = {'Specify a scale factor to resize HDR image, [f<1 to compress]'};
dlgtitle = 'User Input'; dims = [1 50];definput = {'1.0'};
answer = str2double(inputdlg(prompt,dlgtitle,dims,definput));
if isempty(answer)
    ans1 = 1.0;
else
    ans1 = answer(1);
end
newHDR = imresize(HDR, ans1);
end

%% load an luminance map, [from .txt file]
function [Lmap] = loadLuminanceMap()
[fn,pn]=uigetfile('*.txt','load 180 FOV luminance map.');
L = table2array(readtable([pn,fn]));
Lmap = transpose(reshape(L(:, 3)*179, max(L(:,1))+1,[]));
Lmap(Lmap < 0) = 0;
end

function [imF, cf] = loadHDR()
[fn,pn]=uigetfile('*.HDR','load 180 FOV hdr image.');
str=[pn,fn]; imF = hdrread(str); cf1 = ui_cf();
yn = yn_dialog("Calibrate EXP VALUE?");
if ismember(yn, ['Yes', 'yes'])
    exp = getExpValue(str);
else
    exp = 1;
end
cf = cf1./exp;
end