#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include <iostream>
#include <math.h>
#include "FigureCoefficient.h"

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
            for (int i = 0; i < I.rows; ++i)
                for (int j = 0; j < I.cols; ++j) {
                    _I(i, j)[0] = (_I(i, j)[0] / 32) * 32;
                    _I(i, j)[1] = (_I(i, j)[1] / 32) * 32;
                    _I(i, j)[2] = (_I(i, j)[2] / 32) * 32;
                }
            I = _I;
            break;
    }
    return I;
}

cv::Mat selectMax(cv::Mat &I) {
    CV_Assert(I.depth() != sizeof(uchar));
    cv::Mat res(I.rows, I.cols, CV_8UC3);
    switch (I.channels()) {
        case 3:
            cv::Mat_<cv::Vec3b> _I = I;
            cv::Mat_<cv::Vec3b> _R = res;
            for (int i = 0; i < I.rows; ++i)
                for (int j = 0; j < I.cols; ++j) {
                    int sel = (_I(i, j)[0] < _I(i, j)[1]) ? 1 : 0;
                    sel = _I(i, j)[sel] < _I(i, j)[2] ? 2 : sel;
                    _R(i, j)[0] = sel == 0 ? 255 : 0;
                    _R(i, j)[1] = sel == 1 ? 255 : 0;
                    _R(i, j)[2] = sel == 2 ? 255 : 0;
                }
            res = _R;
            break;
    }
    return res;
}

bool computeBordersOfColor(cv::Mat_<cv::Vec3b> &_I, FigureCoefficient *figureCoefficient, int i, int j) {

            int b_min = figureCoefficient->getColors()[0];
            int g_min = figureCoefficient->getColors()[1];
            int r_min = figureCoefficient->getColors()[2];
            int b_max = figureCoefficient->getColors()[3];
            int g_max = figureCoefficient->getColors()[4];
            int r_max = figureCoefficient->getColors()[5];

            bool condition = _I(i, j)[0] != b_min  && _I(i, j)[1] != g_min && _I(i, j)[2] != r_min;
            bool condition1 = (_I(i, j)[0] > b_min &&  _I(i, j)[0] > b_max) &&
                    (_I(i, j)[1] < g_min && _I(i, j)[1] > g_max) && (_I(i, j)[2] < r_min && _I(i, j)[2] > r_max);
    return condition;
}

void computeBox(cv::Mat &I, FigureCoefficient *figureCoefficient) {
    CV_Assert(I.depth() != sizeof(uchar));
    cv::Mat res(I.rows, I.cols, CV_8UC3);
    switch (I.channels()) {
        case 3:
            cv::Mat_<cv::Vec3b> _I = I;

            int coorUpX = 0;
            int coorUpY = 0;
            int coorDownX = I.cols;
            int coorDownY = I.rows;
            int coorLeftX = 0;
            int coorLeftY = I.rows;
            int coorRightX = I.cols;
            int coorRightY = 0;

            int i = 0;
            int j = 0;
            while (computeBordersOfColor(_I, figureCoefficient, i, j) && i < I.rows) {
                j = 0;
                while (computeBordersOfColor(_I, figureCoefficient, i, j) && j < I.cols) {
                    coorUpX = j;
                    coorUpY = i;
                    ++j;
                }
                ++i;
            }

            j = I.cols - 1;
            i = 0;
            while (computeBordersOfColor(_I, figureCoefficient, i, j) && j >= 0) {
                --j;
                i = 0;
                while (computeBordersOfColor(_I, figureCoefficient, i, j) && i < I.rows) {
                    coorRightX = j;
                    coorRightY = i;
                    ++i;
                }
                if (i > 255)
                    i = 255;
            }

            j = 0;
            i = I.rows - 1;
            while (computeBordersOfColor(_I, figureCoefficient, i, j) && j < I.cols) {
                ++j;
                i = I.rows - 1;
                while (computeBordersOfColor(_I, figureCoefficient, i, j) && i >= 0) {
                    coorLeftX = j;
                    coorLeftY = i;
                    --i;
                }
                if (i < 0)
                    i = 0;
            }

            i = I.rows - 1;
            j = I.cols - 1;
            while (computeBordersOfColor(_I, figureCoefficient, i, j) && i >= 0) {
                --i;
                j = I.cols - 1;
                while (computeBordersOfColor(_I, figureCoefficient, i, j) && j >= 0 ) {
                    coorDownX = j;
                    coorDownY = i;
                    j--;
                }
            }

            figureCoefficient->addPixel(coorLeftX - 2, 0);
            figureCoefficient->addPixel(coorUpY - 2, 1);
            figureCoefficient->addPixel(coorRightX + 4, 2);
            figureCoefficient->addPixel(coorDownY + 4, 3);
    }
}

void computeField(FigureCoefficient *figureCoefficient, cv::Mat &I) {
    int field = 0;
    int b = figureCoefficient->getColors()[0];
    int g = figureCoefficient->getColors()[1];
    int r = figureCoefficient->getColors()[2];
    CV_Assert(I.depth() != sizeof(uchar));
    cv::Mat res(I.rows, I.cols, CV_8UC3);
    switch (I.channels()) {
        case 3:
            cv::Mat_<cv::Vec3b> _I = I;
            for (int j = figureCoefficient->getPixels()[0]; j <= figureCoefficient->getPixels()[2]; ++j) {
                for (int i = figureCoefficient->getPixels()[1]; i <= figureCoefficient->getPixels()[3]; ++i) {
                    if (_I(i, j)[0] <= b && _I(i, j)[1] <= g && _I(i, j)[2] <= r) {
                        field = field + 1;
                        //_I(i, j)[0] = 255;
                    }
                }
            }
    }
    figureCoefficient->setField(field);
    std::cout <<  ": S = " << figureCoefficient->getField();
}

void computeCircumference(FigureCoefficient *figureCoefficient, cv::Mat &I) {
    int circumference = 0;
    int b = figureCoefficient->getColors()[0];
    int g = figureCoefficient->getColors()[1];
    int r = figureCoefficient->getColors()[2];
    CV_Assert(I.depth() != sizeof(uchar));
    cv::Mat res(I.rows, I.cols, CV_8UC3);
    switch (I.channels()) {
        case 3:
            cv::Mat_<cv::Vec3b> _I = I;
            for (int j = figureCoefficient->getPixels()[0]; j <= figureCoefficient->getPixels()[2]; ++j) {
                for (int i = figureCoefficient->getPixels()[1]; i <= figureCoefficient->getPixels()[3]; ++i) {
                    if (_I(i, j)[0] <= b && _I(i, j)[1] <= g && _I(i, j)[2] <= r) {
                        if ((_I(i - 1, j)[0] > b && _I(i - 1, j)[1] > g && _I(i - 1, j)[2] > r) ||
                            (_I(i, j - 1)[0] > b && _I(i, j - 1)[1] > g && _I(i, j - 1)[2] > r) ||
                            (_I(i + 1, j)[0] > b && _I(i + 1, j)[1] > g && _I(i + 1, j)[2] > r) ||
                            (_I(i, j + 1)[0] > b && _I(i, j + 1)[1] > g && _I(i, j + 1)[2] > r)) {
                            circumference = circumference + 1;
                            //_I(i, j)[1] = 255;
                        }
                    }
                }
            }
    }
    figureCoefficient->setCircumference(circumference);
    std::cout << ", L = " << figureCoefficient->getCircumference();
}

void computeCoefficientOfMalinowska (FigureCoefficient *figureCoefficient) {
    float coefficientOfMalinowska = (figureCoefficient->getCircumference() / (2 * std::sqrtf(M_PI * figureCoefficient->getField()))) -1;
    figureCoefficient->setW3(coefficientOfMalinowska);
    std::cout << ", W3 = " << figureCoefficient->getW3();
}

void computeAngle(FigureCoefficient *figureCoefficient, double m00, double m01, double m10){
    float middleX = (figureCoefficient->getPixels()[2] - figureCoefficient->getPixels()[0])/ 2 + figureCoefficient->getPixels()[0];
    float middleY = (figureCoefficient->getPixels()[3] - figureCoefficient->getPixels()[1])/ 2 + figureCoefficient->getPixels()[1];
    float middleI = m10/m00;
    float middleJ = m01/m00;

    float value = (middleI - middleX)/(middleJ-middleY);
    float angle =  atan(value);
    float angleRad = angle * 180 / M_PI;
    figureCoefficient->setAngle(angleRad);

}

void computeMoments(FigureCoefficient *figureCoefficient, cv::Mat &I) {

    double m00 = 0;
    double m01 = 0;
    double m10 = 0;
    double m11 = 0;
    double m02 = 0;
    double m20 = 0;

    int b = figureCoefficient->getColors()[0];
    int g = figureCoefficient->getColors()[1];
    int r = figureCoefficient->getColors()[2];

    CV_Assert(I.depth() != sizeof(uchar));
    cv::Mat res(I.rows, I.cols, CV_8UC3);
    switch (I.channels()) {
        case 3:
            cv::Mat_<cv::Vec3b> _I = I;
            for (int j = figureCoefficient->getPixels()[0]; j <= figureCoefficient->getPixels()[2]; ++j) {
                for (int i = figureCoefficient->getPixels()[1]; i <= figureCoefficient->getPixels()[3]; ++i) {
                    if (_I(i, j)[0] == b && _I(i, j)[1] == g && _I(i, j)[2] == r) {
                        m00 = m00 + pow(i, 0) * pow(j, 0);
                        m01 = m01 + pow(i, 0) * pow(j, 1);
                        m10 = m10 + pow(i, 1) * pow(j, 0);
                        m11 = m11 + pow(i, 1) * pow(j, 1);
                        m02 = m02 + pow(i, 0) * pow(j, 2);
                        m20 = m20 + pow(i, 2) * pow(j, 0);
                    }
                }
            }
    }

    float M20 = m20 - (pow(m10, 2))/m00;
    float M02 = m02 - (pow(m01, 2))/m00;
    float M11 = m11 - m10*m01/m00;

    float M1 = (M20 + M02)/pow(m00, 2);
    float M7 = (M20 * M02 - pow(M11, 2))/pow(m00, 4);

    figureCoefficient->setM1(M1);
    figureCoefficient->setM7(M7);

    computeAngle(figureCoefficient, m00, m01, m10);

}




int main(int, char *[]) {
    std::cout << "Start ..." << std::endl;

    // zadanie 1
    /*std::string fileNames[5] = {"elipsa", "elipsa1", "kolo", "prost", "troj"};
    for (int i = 0; i < 5; i++) {
        FigureCoefficient figureCoefficient;
        figureCoefficient.addColor(0,0);
        figureCoefficient.addColor(0,1);
        figureCoefficient.addColor(0,2);
        figureCoefficient.setFileName(fileNames[i]);
        cv::Mat image = cv::imread("SourceImages/" +  fileNames[i] + ".dib");

        computeBox(image, &figureCoefficient);
        cv::Mat image2 = image(cv::Rect(figureCoefficient.getPixels()[0], figureCoefficient.getPixels()[1],
                figureCoefficient.getPixels()[2] - figureCoefficient.getPixels()[0], figureCoefficient.getPixels()[3]
                - figureCoefficient.getPixels()[1]));
        std::cout << i + 1 << ". Plik " + fileNames[i];
        computeField(&figureCoefficient, image);
        computeCircumference(&figureCoefficient, image);
        computeCoefficientOfMalinowska(&figureCoefficient);
        computeMoments(&figureCoefficient, image);
        std::cout << ", M1 = " << figureCoefficient.getM1() << ", M7 = " << figureCoefficient.getM7() << std::endl;

        cv::imshow("Shape", image);
        cv::waitKey(-1);
    }*/

    //zadanie 2
    cv::Mat imageArrow1 = cv::imread("SourceImages/strzalki_2.dib");
    for (int i = 0; i <= 180; i = i + 45) {
        std::string figureName = "Strzalka R = " + std::to_string(i);

        FigureCoefficient figureCoefficient;
        figureCoefficient.addColor(0,0);
        figureCoefficient.addColor(180-i,1);
        figureCoefficient.addColor(i,2);
        figureCoefficient.addColor(0,0);
        figureCoefficient.addColor(180-i,1);
        figureCoefficient.addColor(i,2);


        computeBox(imageArrow1, &figureCoefficient);
        computeMoments(&figureCoefficient, imageArrow1);


        std::cout << figureName + ", nachylenie " << figureCoefficient.getAngle();
        figureCoefficient.setFileName(figureName);

        computeField(&figureCoefficient, imageArrow1);
        computeCircumference(&figureCoefficient, imageArrow1);
        computeCoefficientOfMalinowska(&figureCoefficient);

        std::cout << ", M1 = " << figureCoefficient.getM1() << ", M7 = " << figureCoefficient.getM7() << std::endl;
        cv::Mat image2 = imageArrow1(cv::Rect(figureCoefficient.getPixels()[0], figureCoefficient.getPixels()[1],
                                              figureCoefficient.getPixels()[2] - figureCoefficient.getPixels()[0], figureCoefficient.getPixels()[3]
                                                                                                                   - figureCoefficient.getPixels()[1]));
        cv::imshow("Shape", imageArrow1);
        cv::waitKey(-1);
    }


    return 0;
}
