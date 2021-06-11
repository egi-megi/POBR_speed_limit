#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include <iostream>
#include <numeric>

int amountOfPixelsArrayO[8] = {0};

cv::Mat &perform(cv::Mat &I) {
    CV_Assert(I.depth() != sizeof(uchar));
    switch (I.channels()) {
        case 1:
            for (int i = 0; i < I.rows; ++i)
                for (int j = 0; j < I.cols; ++j)
                    I.at<uchar>(i, j) = (I.at<uchar>(i, j) / 32) * 32;
            break;
        case 3:
            cv::Mat_<cv::Vec3b> _I = I;
            int brightness = 0;
            for (int i = 0; i < I.rows; ++i) {
                for (int j = 0; j < I.cols; ++j) {
                    // contrast at the top
                    if (j > i && j < I.rows - i) {
                        float B = _I(i, j)[0] * 1.15;
                        float G = _I(i, j)[1] * 1.15;
                        float R = _I(i, j)[2] * 1.15;
                        if (B < 0) {B = 0;}
                        if (G < 0) {G = 0;}
                        if (R < 0) {R = 0;}
                        if (B > 255) {B = 255;}
                        if (G > 255) {G = 255;}
                        if (R > 255) {R = 255;}
                        _I(i, j)[0] = B;
                        _I(i, j)[1] = G;
                        _I(i, j)[2] = R;
                        brightness = (_I(i, j)[0] + _I(i, j)[1] + _I(i, j)[2]) / 3;
                        amountOfPixelsArrayO[brightness / 32] = amountOfPixelsArrayO[brightness / 32] + 1;
                    } else if (j < i && j > I.rows - i) {
                        float B = _I(i, j)[0] + 50;
                        float G = _I(i, j)[1] + 50;
                        float R = _I(i, j)[2] + 50;
                        if (B < 0) {B = 0;}
                        if (G < 0) {G = 0;}
                        if (R < 0) {R = 0;}
                        if (B > 255) {B = 255;}
                        if (G > 255) {G = 255;}
                        if (R > 255) {R = 255;}
                        _I(i, j)[0] = B;
                        _I(i, j)[1] = G;
                        _I(i, j)[2] = R;
                        brightness = (_I(i, j)[0] + _I(i, j)[1] + _I(i, j)[2]) / 3;
                        amountOfPixelsArrayO[brightness / 32] = amountOfPixelsArrayO[brightness / 32] + 1;
                    } else if (j < i && j < I.rows - i) {
                        float grey = (_I(i, j)[0] + _I(i, j)[1] + _I(i, j)[2]) / 3;
                        float gr = _I(i, j)[2] - grey;
                        if (gr < 0) {gr = 0;}
                        if (gr > 255) {gr = 255;}
                        _I(i, j)[0] = gr;
                        _I(i, j)[1] = gr;
                        _I(i, j)[2] = gr;
                        brightness = (_I(i, j)[0] + _I(i, j)[1] + _I(i, j)[2]) / 3;
                        amountOfPixelsArrayO[brightness / 32] = amountOfPixelsArrayO[brightness / 32] + 1;
                    } else {
                        _I(i, j)[0] = _I(i, j)[0];
                        _I(i, j)[1] = _I(i, j)[1];
                        _I(i, j)[2] = _I(i, j)[2];
                        brightness = (_I(i, j)[0] + _I(i, j)[1] + _I(i, j)[2]) / 3;
                        amountOfPixelsArrayO[brightness / 32] = amountOfPixelsArrayO[brightness / 32] + 1;
                    }
                }
            }
            I = _I;
            break;
    }
    return I;
}


void displayHistogram() {

    int numberOfLevels = 8;
    int maxValueOfPixel = 256;
    int numberOfPixelsInLevel = maxValueOfPixel / numberOfLevels;
    for (int i = 0; i < (sizeof(amountOfPixelsArrayO) / sizeof(*amountOfPixelsArrayO)); i++) {
        std::cout << "Pikseli o jaskości " << i * numberOfPixelsInLevel << "-" << (i + 1) * numberOfPixelsInLevel - 1 << " jest "
                  << amountOfPixelsArrayO[i] << std::endl;
    }

    long sumOfPixels = std::accumulate(amountOfPixelsArrayO, amountOfPixelsArrayO + 8, 0);
    std::cout << "Suma pikseli " << sumOfPixels << std::endl;
}


int main(int, char *[]) {
    std::cout << "Start ..." << std::endl;
    cv::Mat image = cv::imread("Lena.png");
    perform(image);
    displayHistogram();
    cv::imshow("Max", image);
    cv::imshow("Lena", image);
    cv::imwrite("Max.png", image);
    cv::waitKey(-1);
    return 0;
}
