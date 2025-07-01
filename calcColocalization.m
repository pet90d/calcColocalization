function [r,m1orig,m2orig,m1,m2,m1area,m2area,numOfPixelPairsinOrMask,label1area,label2area,label12area]=calcColocalization(image1,image2,n,ci,thr1,thr2,whatToCalculate,printReport,showPlots)
% INPUT:
% image1, image2 - you should know..., the two images to be analyzed
% n - number of times scrambling is done
% ci - confidence interval % (0-1)
% thr1, thr2 - thresholds for the two images OR masks for the two images
%              THE CALCULATION IS CARRIED OUT WITH PIXELS THAT ARE IN
%              EITHER ONE OF THE MASKS
% whatToCalculate - 'pearson', 'manders' or 'both'
% printReport - if true, a report will be displayed on the screen
% showPlots - if true, the confidence interval plots will be displayed
% OUTPUT
% r, m1orig, m2orig, m1, m2, m1area, m2area - structures containing the
%    Pearson correlation coefficient, the original Manders coefficients (m1orig,
%    m2orig), the intensity-based Manders coefficients (m1, m2) and the
%    area-based Manders coefficients (m1area, m2area). Each of them is a
%    structure with the following fields:
%       - value: the value calculated for the image pair
%       - array: the values calculated for the scrambled images
%       - confInterval: the confidence interval for no correlation
%                     (confidence interval % provided as an input to the
%                     function)
%       - hist: a histogram of the values present in 'array'
%       M1 IS THE FRACTION OF THE INTERSECTION COMPARED TO THE PIXELS IN
%       IMAGE1 ABOVE THE THRESHOLD, M2 IS THE FRACTION OF THE INTERSECTION
%       TO THE PIXELS IN IMAGE2 ABOVE THE THRESHOLD
%  numOfPixelPairs - the number of pixels pairs present in the mask in which
%    the analysis is performed (determined by the masks and the 'typeOfBoolean'
% label1area, label2area, label12area - area of the mask of image1, image2
%    and the area of the joint mask, respectively
% Written by Peter Nagy, email: peter.v.nagy@gmail.com, web: https://peternagyweb.hu
% V2.03
if or(ci>1,ci<0)
    errordlg('Confidence interval has to be <1 and >0','Oops');
    r=nan;
    m1orig=nan;
    m2orig=nan;
    m1=nan;
    m2=nan;
    m1area=nan;
    m2area=nan;
    numOfPixelPairsinOrMask=nan;
    label1area=nan;
    label2area=nan;
    label12area=nan;
    return;
end
if n<10
    errordlg('Scrambling has to be done at least 10x','Oops');
    r=nan;
    m1orig=nan;
    m2orig=nan;
    m1=nan;
    m2=nan;
    m1area=nan;
    m2area=nan;
    numOfPixelPairsinOrMask=nan;
    label1area=nan;
    label2area=nan;
    label12area=nan;
    return;
end
if ~ismember(whatToCalculate,{'pearson','manders','both'})
    errordlg('''whatToCalculate'' has to be ''pearson'', ''manders'' or ''both''.','Oops');
    r=nan;
    m1orig=nan;
    m2orig=nan;
    m1=nan;
    m2=nan;
    m1area=nan;
    m2area=nan;
    numOfPixelPairsinOrMask=nan;
    label1area=nan;
    label2area=nan;
    label12area=nan;
    return;
end
image1type=whos('image1');
image2type=whos('image2');
if ~strcmp(image1type.class,'dip_image') || ~strcmp(image2type.class,'dip_image')
    errordlg('Type of images has to be ''dip_image''.','Oops');
    r=nan;
    m1orig=nan;
    m2orig=nan;
    m1=nan;
    m2=nan;
    m1area=nan;
    m2area=nan;
    numOfPixelPairsinOrMask=nan;
    label1area=nan;
    label2area=nan;
    label12area=nan;
    return;
end
if isnumeric(thr1)~=isnumeric(thr2)
    errordlg('Either thresholds are provided for both image, or masks for both','Oops');
    r=nan;
    m1orig=nan;
    m2orig=nan;
    m1=nan;
    m2=nan;
    m1area=nan;
    m2area=nan;
    numOfPixelPairsinOrMask=nan;
    label1area=nan;
    label2area=nan;
    label12area=nan;
    return;
end
image1data_2d=double(image1);
image2data_2d=double(image2);
if numel(size(image1))~=numel(size(image2)) || ~all(size(image1)==size(image2))
    errordlg('Image sizes have to be identical','Oops');
    r=nan;
    m1orig=nan;
    m2orig=nan;
    m1=nan;
    m2=nan;
    m1area=nan;
    m2area=nan;
    numOfPixelPairsinOrMask=nan;
    label1area=nan;
    label2area=nan;
    label12area=nan;
    return;
end
if isnumeric(thr1) % threshold is provided
    if numel(thr1)>1 || numel(thr2)>2
        errordlg('Threshold can only contain a single number if it is numeric','Oops');
        r=nan;
        m1orig=nan;
        m2orig=nan;
        m1=nan;
        m2=nan;
        m1area=nan;
        m2area=nan;
        numOfPixelPairsinOrMask=nan;
        label1area=nan;
        label2area=nan;
        label12area=nan;
        return;
    end
    orMask_1d=image1data_2d(:)>thr1 | image2data_2d(:)>thr2;
    bothPixeldata_1d=[image1data_2d(orMask_1d) image2data_2d(orMask_1d)];
    bothMaskdata_1d=[image1data_2d(orMask_1d)>thr1 image2data_2d(orMask_1d)>thr2];
    numOfPixelPairsinOrMask=size(bothPixeldata_1d,1);
    bothAllMaskdata_1d=[image1data_2d(:)>thr1 image2data_2d(:)>thr2];
else % mask is provided
    thr1data_2d=double(thr1)>0;
    thr2data_2d=double(thr2)>0;
    orMask_1d=thr1data_2d(:)>0 | thr2data_2d(:)>0;
    thr1ForOrMask=thr1data_2d(orMask_1d);
    thr2ForOrMask=thr2data_2d(orMask_1d);
    bothPixeldata_1d=[image1data_2d(orMask_1d) image2data_2d(orMask_1d)];
    numOfPixelPairsinOrMask=size(bothPixeldata_1d,1);
    bothMaskdata_1d=false(numOfPixelPairsinOrMask,2);
    bothMaskdata_1d(thr1ForOrMask,1)=true;
    bothMaskdata_1d(thr2ForOrMask,2)=true;
    bothAllMaskdata_1d=[thr1data_2d(:) thr2data_2d(:)];
end
bothAllImagedata_1d=[image1data_2d(:) image2data_2d(:)];
numOfAllPixels=numel(image1data_2d);
r.array=zeros(n,1);
m1orig.array=zeros(n,1);
m2orig.array=zeros(n,1);
m1.array=zeros(n,1);
m2.array=zeros(n,1);
m1area.array=zeros(n,1);
m2area.array=zeros(n,1);
% calculate r and the Manders coefficients for the unscrambled images
if ismember(whatToCalculate,{'pearson','both'})
    rTemp=corr(bothPixeldata_1d);
    r.value=rTemp(1,2);
end
if ismember(whatToCalculate,{'manders','both'})
    [m1orig.value,m2orig.value,m1.value,m2.value,m1area.value,m2area.value,label1area,label2area,label12area]=calculateManders(bothPixeldata_1d(:,1),bothPixeldata_1d(:,2),bothMaskdata_1d(:,1),bothMaskdata_1d(:,2),'mask');
end
% scrambling
for j=1:n
    if ismember(whatToCalculate,{'pearson','both'})
        scrambled=zeros(numOfPixelPairsinOrMask,2);
        permutedIndices1=randperm(numOfPixelPairsinOrMask);
        scrambled(:,2)=bothPixeldata_1d(:,2); % the 2nd column is unscrambled
        scrambled(:,1)=bothPixeldata_1d(permutedIndices1,1); % this is the scrabled column
        rtemp=corr(scrambled);
        r.array(j)=rtemp(1,2);
    end
    % for Manders randomly choose pixels from the whole image
    if ismember(whatToCalculate,{'manders','both'})
        permutedIndices2=randperm(numOfAllPixels);
        permutedIndices3=randperm(numOfAllPixels);
        permutedIndices2=permutedIndices2(1:numOfPixelPairsinOrMask);
        permutedIndices3=permutedIndices3(1:numOfPixelPairsinOrMask);
        [m1orig.array(j),m2orig.array(j),m1.array(j),m2.array(j),m1area.array(j),m2area.array(j)]=calculateManders(bothAllImagedata_1d(permutedIndices2,1),bothAllImagedata_1d(permutedIndices3,2),bothAllMaskdata_1d(permutedIndices2,1),bothAllMaskdata_1d(permutedIndices3,2),'mask');
    end
    % end of Manders block
end
% calculate histograms and confidence intervals, generate plots
if ismember(whatToCalculate,{'pearson','both'}) && sum(isnan(r.array))~=0
    errordlg('NaN generated in the calculation, probably because the threshold is too high.','Oops','modal');
    r.hist=[nan nan];
    r.confInterval=[nan nan];
    numOfPixelPairsinOrMask=nan;
else
    if ismember(whatToCalculate,{'pearson','both'})
        [r.hist,r.confInterval]=calcHistAndConfInt(r.array,n,ci);
    end
    % Manders block
    if ismember(whatToCalculate,{'manders','both'})
        [m1orig.hist,m1orig.confInterval]=calcHistAndConfInt(m1orig.array,n,ci);
        [m2orig.hist,m2orig.confInterval]=calcHistAndConfInt(m2orig.array,n,ci);
        [m1.hist,m1.confInterval]=calcHistAndConfInt(m1.array,n,ci);
        [m2.hist,m2.confInterval]=calcHistAndConfInt(m2.array,n,ci);
        [m1area.hist,m1area.confInterval]=calcHistAndConfInt(m1area.array,n,ci);
        [m2area.hist,m2area.confInterval]=calcHistAndConfInt(m2area.array,n,ci);
    end
    % end of Manders block
    if showPlots>0
        figure;
        switch whatToCalculate
            case 'pearson'
                n1=1;
                n2=1;
            case 'manders'
                n1=2;
                n2=3;
            case 'both'
                n1=2;
                n2=4;
        end
        plotCounter=0;
        if ismember(whatToCalculate,{'pearson','both'})
            plotCounter=plotCounter+1;
            subplot(n1,n2,plotCounter);
            bar(r.hist(:,1),r.hist(:,2));
            title('r');
        end
        if ismember(whatToCalculate,{'manders','both'})
            plotCounter=plotCounter+1;
            subplot(n1,n2,plotCounter);
            bar(m1orig.hist(:,1),m1orig.hist(:,2));
            title('Original M1');
            plotCounter=plotCounter+1;
            subplot(n1,n2,plotCounter);
            bar(m2orig.hist(:,1),m2orig.hist(:,2));
            title('Original M2');
            plotCounter=plotCounter+1;
            subplot(n1,n2,plotCounter);
            bar(m1.hist(:,1),m1.hist(:,2));
            title('Modified M1');
            plotCounter=plotCounter+1;
            subplot(n1,n2,plotCounter);
            bar(m2.hist(:,1),m2.hist(:,2));
            title('Modified M2');
            plotCounter=plotCounter+1;
            subplot(n1,n2,plotCounter);
            bar(m1area.hist(:,1),m1area.hist(:,2));
            title('Area M1');
            plotCounter=plotCounter+1;
            subplot(n1,n2,plotCounter);
            bar(m2area.hist(:,1),m2area.hist(:,2));
            title('Area M2');
        end
    end
end
switch whatToCalculate
    case 'pearson'
        m1orig=nan;
        m2orig=nan;
        m1=nan;
        m2=nan;
        m1area=nan;
        m2area=nan;
        label1area=nan;
        label2area=nan;
        label12area=nan;
    case 'manders'
        r=nan;
end
if printReport>0
    % report
    if ismember(whatToCalculate,{'pearson','both'})
        fprintf('The correlation coefficient (r) is %f.\nThe %1.0f%% confidence range of r, assuming r=0, is between %f and %f.\n',r.value,100*ci, r.confInterval(1),r.confInterval(2));
    end
    if ismember(whatToCalculate,{'manders','both'})
        fprintf('Original Manders coefficients (mean [%1.0f%% confidence range assuming the absence of correlation])\n',100*ci);
        fprintf('M1=ratio of intensity in image 1 for which intensity in image 2 is >thr to total intensity in image 1\n');
        fprintf('M2=ratio of intensity in image 2 for which intensity in image 1 is >thr to total intensity in image 2\n');
        fprintf('M1=\t%f [%f, %f]\n',m1orig.value,m1orig.confInterval(1),m1orig.confInterval(2));
        fprintf('M2=\t%f [%f, %f]\n',m2orig.value,m2orig.confInterval(1),m2orig.confInterval(2));
        fprintf('Modified Manders coefficients (mean [%1.0f%% confidence range assuming the absence of correlation])\n',100*ci);
        fprintf('M1=ratio of intensity in image 1 for which intensity in both images is >thr to intensity in image 1 above the thr\n');
        fprintf('M2=ratio of intensity in image 2 for which intensity in both images is >thr to intensity in image 2 above the thr\n');
        fprintf('M1=\t%f [%f, %f]\n',m1.value,m1.confInterval(1),m1.confInterval(2));
        fprintf('M2=\t%f [%f, %f]\n',m2.value,m2.confInterval(1),m2.confInterval(2));
        fprintf('Manders coefficients for areas (mean [%1.0f%% confidence range assuming the absence of correlation])\n',100*ci);
        fprintf('The same as the modified Manders coefficients, but instead of intensity sums area sums are calculated\n');
        fprintf('M1=\t%f [%f, %f]\n',m1area.value,m1area.confInterval(1),m1area.confInterval(2));
        fprintf('M2=\t%f [%f, %f]\n',m2area.value,m2area.confInterval(1),m2area.confInterval(2));
    end
    switch whatToCalculate
        case 'pearson'
            fprintf('Only the Pearson correlation coefficient was calcualted.\n');
        case 'manders'
            fprintf('Only the Manders coefficient was calcualted.\n');
        case 'both'
            fprintf('Both the Pearson and the Manders coefficients were calcualted.\n');
    end
    if isnumeric(thr1)
        fprintf('Threshold for image1 and image2: %1.0f, %1.0f\n',thr1,thr2);
    else
        fprintf('Masks were provided.\n');
    end
    fprintf('Number of pixels in the thresholded images: %1.0f\n',numOfPixelPairsinOrMask);
    if ismember(whatToCalculate,{'manders','both'})
        fprintf('Number of pixels above threshold in image1: %1.0f\n',label1area);
        fprintf('Number of pixels above threshold in image2: %1.0f\n',label2area);
        fprintf('Number of pixels above threshold in image1 and image2: %1.0f\n',label12area);
    end
    fprintf('Scrambling was done %1.0f-times.\n\n',n);
end

function [hist,confInt]=calcHistAndConfInt(array,n,ci)
k=1+3.322*log10(n);
min_bin=prctile(array,1);
max_bin=prctile(array,99);
bins=min_bin:(max_bin-min_bin)/k:max_bin;
if ~isempty(bins)
    b=histcounts(array,bins);
    hist=[bins(1:end-1)' b'];
else
    hist=nan;
end
confInt(1)=prctile(array,100*(1-ci)/2);
confInt(2)=prctile(array,100*(1-(1-ci)/2));

function [m1orig,m2orig,m1,m2,m1area,m2area,label1area,label2area,label12area]=calculateManders(image1data,image2data,thr1,thr2,thresholdOrMask)
% image1data,image2data - 1D column arrays of image data
% thr1,thr2 - threshold values OR 1D column arrays of the mask
switch thresholdOrMask
    case 'threshold'
        m1orig=sum(image1data(image2data>thr2))/sum(image1data(:));
        m2orig=sum(image2data(image1data>thr1))/sum(image2data(:));
        m1=sum(image1data(image1data>thr1 & image2data>thr2))/sum(image1data(image1data>thr1));
        m2=sum(image2data(image1data>thr1 & image2data>thr2))/sum(image2data(image2data>thr2));
        m1area=sum(image1data>thr1 & image2data>thr2)/sum(image1data>thr1);
        m2area=sum(image1data>thr1 & image2data>thr2)/sum(image2data>thr2);
        label1area=sum(image1data>thr1);
        label2area=sum(image2data>thr2);
        label12area=sum(image1data>thr1 & image2data>thr2);
    case 'mask'
        m1orig=sum(image1data(thr2>0))/sum(image1data(:));
        m2orig=sum(image2data(thr1>0))/sum(image2data(:));
        m1=sum(image1data(thr1>0 & thr2>0))/sum(image1data(thr1>0));
        m2=sum(image2data(thr1>0 & thr2>0))/sum(image2data(thr2>0));
        m1area=sum(thr1>0 & thr2>0)/sum(thr1>0);
        m2area=sum(thr1>0 & thr2>0)/sum(thr2>0);
        label1area=sum(thr1>0);
        label2area=sum(thr2>0);
        label12area=sum(thr1>0 & thr2>0);
end
