// Author: Agnieszka Jurkiewicz

#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include <iostream>
#include <numeric>
#include <algorithm>
#include "Pixel.h"

void makeNewPixelAfterFiltering(int coordinateX, int coordinateY, cv::Mat_<cv::Vec3b> _I, float *blurCoefficients, float* colors, int sum) {
    colors[0]=0.0;
    colors[1]=0.0;
    colors[2]=0.0;
    float B = 0.0;
    float G = 0.0;
    float R = 0.0;
    int tabNum = 0;
    for (int k = coordinateY - 2; k < coordinateY + 3; k++) {
        for (int l = coordinateX - 2; l < coordinateX + 3; l++) {
            B = B + _I(l,k)[0] * blurCoefficients[tabNum];
            G = G + _I(l,k)[1] * blurCoefficients[tabNum];
            R = R + _I(l,k)[2] * blurCoefficients[tabNum];
            //std::cout << blurCoefficients[tabNum] <<std::endl;
            tabNum++;

        }
    }

    B = B/sum;
    G = G/sum;
    R = R/sum;

    if (B < 0) {B = 0;}
    if (G < 0) {G = 0;}
    if (R < 0) {R = 0;}
    if (B > 255) {B = 255;}
    if (G > 255) {G = 255;}
    if (R > 255) {R = 255;}

    colors[0] = B;
    colors[1] = G;
    colors[2] = R;

}

cv::Mat blurAndSharpFiltering(cv::Mat& I, float *blurCoefficients, int filterLength){
    CV_Assert(I.depth() != sizeof(uchar));
    cv::Mat  res(I.rows,I.cols, CV_8UC3);

    int sum = 0;
    sum = std::accumulate(blurCoefficients, blurCoefficients+filterLength , sum);
    if (sum == 0) {
        sum = 1;
    }

    switch(I.channels())  {
        case 3:
            cv::Mat_<cv::Vec3b> _I = I;
            cv::Mat_<cv::Vec3b> _R = res;
            float colorsRGB[3];
            for( int i = 1; i < I.rows - 2; i++)
                for( int j = 1; j < I.cols - 2; j++){
                    makeNewPixelAfterFiltering(i, j, _I, blurCoefficients, colorsRGB, sum);
                    _R(i,j)[0] = colorsRGB[0];
                    _R(i,j)[1] = colorsRGB[1];
                    _R(i,j)[2] = colorsRGB[2];
                }
            res = _R;
            break;
    }
    return res;
}

void makeFilters(int filterLength, float mainElement, float centerElement, float *filter) {
    for (int i = 0; i < filterLength; i++) {
        filter[i] = mainElement;
        if (i == filterLength/2) {
            filter[i] = centerElement;
        }
    }
}




cv::Mat makeRangeFilters(cv::Mat& I, int size, int kindOfFilter) {
    CV_Assert(I.depth() != sizeof(uchar));
    cv::Mat res(I.rows,I.cols, CV_8UC3);

    int pixelNumberAfterFiltering = 0;
    if (kindOfFilter == 3) {
        pixelNumberAfterFiltering = size * size - 1;
    } else if (kindOfFilter == 2) {
        pixelNumberAfterFiltering = (size * size)/2;
    } else {
        pixelNumberAfterFiltering = 0;
    }

    switch(I.channels())  {
        case 3:
            cv::Mat_<cv::Vec3b> _I = I;
            cv::Mat_<cv::Vec3b> _R = res;

            int filterSize = size * size;

            Pixel lightness[filterSize];

            // iterate after all pixels
            for( int i = size/2; i < I.rows - size/2; i++) {
                for (int j = size / 2; j < I.cols - size / 2; j++) {
                    // iterate only through the surroundings of this one pixel
                    int lightnessIter = 0;
                    for (int k = i - size / 2; k < i + size / 2 + 1; k++) {
                        for (int l = j - size / 2; l < j + size / 2 + 1; l++) {
                            float lightnessValue = _I(k,l)[0] + _I(k,l)[1] + _I(k,l)[2];
                            Pixel pixel(k, l, lightnessValue);
                            lightness[lightnessIter] = pixel;
                            lightnessIter++;
                        }
                    }

                    std::sort(lightness, lightness + filterSize,
                              [](Pixel const & a, Pixel const & b) -> bool
                              { return a.getLightness() < b.getLightness(); } );

                    _R(i, j)[0] = _I(lightness[pixelNumberAfterFiltering].getCoordinateX(),
                                     lightness[pixelNumberAfterFiltering].getCoordinateY())[0];
                    _R(i, j)[1] = _I(lightness[pixelNumberAfterFiltering].getCoordinateX(),
                                     lightness[pixelNumberAfterFiltering].getCoordinateY())[1];
                    _R(i, j)[2] = _I(lightness[pixelNumberAfterFiltering].getCoordinateX(),
                                     lightness[pixelNumberAfterFiltering].getCoordinateY())[2];
                }
            }
            res = _R;
            break;
    }
    return res;
}


int main(int, char *[]) {
    std::cout << "Start ..." << std::endl;
    cv::Mat image = cv::imread("Lena.png");

    int filterLength = 25;
    float blurCoefficients[filterLength];
    float findEdges[filterLength];
    makeFilters(filterLength, 1.0, 1.0, blurCoefficients);
    makeFilters(filterLength, -1.0, 24.0, findEdges);

    cv::Mat blur = blurAndSharpFiltering(image, blurCoefficients, filterLength);
    cv::Mat edges = blurAndSharpFiltering(image, findEdges, filterLength);

    int sizeR = 4;
    int kindOfFilterR = 5;

    while (sizeR%2 == 0) {
        std::cout << "Wpisz rozmiar filtru (wartość dodatnią nieparzystą) i naciśnij Enter: ";
        std::cin >> sizeR;
    }
    std::cout << std::endl;
    while (!(kindOfFilterR == 1 || kindOfFilterR == 2 || kindOfFilterR == 3)) {
        std::cout
                << "Wpisz symbol filtru rankingowego: 1 dla minimum, 2 dla mediany, 3 dla maximum i naciśnij Enter: ";
        std::cin >> kindOfFilterR;
    }
    cv::Mat range = makeRangeFilters(image, sizeR, kindOfFilterR);

    cv::imshow("Lena",image);
    cv::imshow("Blur", blur);
    cv::imshow("Edges", edges);
    cv::imshow("Range", range);
    //std::cout << image2.isContinuous() << blur.isContinuous() << std::endl;
    cv::imwrite("Blur.png", blur);
    cv::imwrite("Edges.png", edges);
    cv::imwrite("Range.png", range);
    cv::waitKey(-1);
    return 0;
}


