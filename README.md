[r,m1orig,m2orig,m1,m2,m1area,m2area,numOfPixelPairsinOrMask,label1area,label2area,label12area]=calcColocalization(image1,image2,n,ci,thr1,thr2,whatToCalculate,printReport,showPlots)
The program evaluations the colocalization in images image1 and image2.
n - The program provides confidence intervals for the case of no correlation by calculating the coefficients between one of the original image and a scrambled version of the other image. n determines how many times this scrambling is performed and therefore from how many random colocalizations the confidence interval of no-colocalization is calculated
ci - The confidence interval to be determined (0-1 corresponding to confidence interval of 0% and 100%)
thr1, thr2 - thresholds for the two images OR masks for the two images. THE CALCULATION IS CARRIED OUT WITH PIXELS THAT ARE IN EITHER ONE OF THE MASKS
whatToCalculate - 'pearson', 'manders' or 'both'
printReport - if true, a report will be displayed on the screen
showPlots - if true, the confidence interval plots will be displayed

OUTPUT:
r, m1orig, m2orig, m1, m2, m1area, m2area - structures containing the Pearson correlation coefficient, the original Manders coefficients (m1orig, m2orig), the intensity-based Manders coefficients (m1, m2) and the area-based Manders coefficients (m1area, m2area). Each of them is a structure with the following fields:
 - value: the value calculated for the image pair
 - array: the values calculated for the scrambled images
 - confInterval: the confidence interval for no correlation (confidence interval % provided as an input to the function)
 - hist: a histogram of the values present in 'array'
  M1 IS THE FRACTION OF THE INTERSECTION COMPARED TO THE PIXELS IN IMAGE1 ABOVE THE THRESHOLD, M2 IS THE FRACTION OF THE INTERSECTION TO THE PIXELS IN IMAGE2 ABOVE THE THRESHOLD
numOfPixelPairs - the number of pixels pairs present in the mask in which the analysis is performed (determined by the masks and the 'typeOfBoolean'
label1area, label2area, label12area - area of the mask of image1, image2 and the area of the joint mask, respectively
