#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include <opencv2/imgproc.hpp>
#include <iostream>
#include <math.h>
#include "FigureCoefficient.h"
#include <vector>

cv::Mat makeWhiteBoard(cv::Mat &I) {
    CV_Assert(I.depth() != sizeof(uchar));
    cv::Mat res(I.rows, I.cols, CV_8UC3);
    switch (I.channels()) {
        case 3:
            cv::Mat_<cv::Vec3b> _R = res;
            for (int x = 0; x < I.rows; ++x)
                for (int y = 0; y < I.cols; ++y) {
                    _R(x, y)[0] = 255;
                    _R(x, y)[1] = 255;
                    _R(x, y)[2] = 255;
                }
            res = _R;
            break;
    }
    return res;
}

bool computeConditionForColor(cv::Mat_<cv::Vec3b> &_I, FigureCoefficient *figureCoefficient, int x, int y) {

    int b_min = figureCoefficient->getColors()[0];
    int g_min = figureCoefficient->getColors()[1];
    int r_min = figureCoefficient->getColors()[2];
    int b_max = figureCoefficient->getColors()[3];
    int g_max = figureCoefficient->getColors()[4];
    int r_max = figureCoefficient->getColors()[5];

    bool condition = (_I(x, y)[0] >= b_min && _I(x, y)[0] <= b_max) &&
                     (_I(x, y)[1] >= g_min && _I(x, y)[1] <= g_max) && (_I(x, y)[2] >= r_min && _I(x, y)[2] <= r_max);
    return condition;
}

void computeBox(cv::Mat &I, FigureCoefficient *figureCoefficient) {
    CV_Assert(I.depth() != sizeof(uchar));
    switch (I.channels()) {
        case 3:
            cv::Mat_<cv::Vec3b> _I = I;

            int x = 0;
            int y = 0;
            int coorUpX = x;
            int coorUpY = y;
            while (!computeConditionForColor(_I, figureCoefficient, x, y) && x < I.rows) {
                y = 0;
                ++x;
                while (!computeConditionForColor(_I, figureCoefficient, x, y) && y < I.cols) {
                    coorUpX = x;
                    coorUpY = y;
                    ++y;
                }
                if (y >= I.cols) {
                    y = I.cols - 1;
                }

            }

            x = 0;
            y = I.cols - 1;
            int coorRightX = x;
            int coorRightY = y;
            while (!computeConditionForColor(_I, figureCoefficient, x, y) && y >= 0) {
                --y;
                x = 0;
                while (!computeConditionForColor(_I, figureCoefficient, x, y) && x < I.rows) {
                    coorRightX = x;
                    coorRightY = y;
                    ++x;
                }
                if (x >= I.rows)
                    x = I.rows - 1;
            }

            x = 0;
            y = 0;
            int coorLeftX = x;
            int coorLeftY = y;
            while (!computeConditionForColor(_I, figureCoefficient, x, y) && y < I.cols) {
                ++y;
                x = 0;
                while (!computeConditionForColor(_I, figureCoefficient, x, y) && x < I.rows) {
                    coorLeftX = x;
                    coorLeftY = y;
                    ++x;
                }
                if (x == I.rows) {
                    x = I.rows - 1;
                }
            }

            x = I.rows - 1;
            y = I.cols - 1;
            int coorDownX = x;
            int coorDownY = y;
            while (!computeConditionForColor(_I, figureCoefficient, x, y) && x >= 0) {
                --x;
                y = I.cols - 1;
                while (!computeConditionForColor(_I, figureCoefficient, x, y) && y >= 0) {
                    coorDownX = x;
                    coorDownY = y;
                    y--;
                }
                if (y < 0) {
                    y = 0;
                }
            }

            figureCoefficient->setCoorMinX(coorUpX);
            figureCoefficient->setCoorMinY(coorLeftY);
            figureCoefficient->setCoorMaxX(coorDownX);
            figureCoefficient->setCoorMaxY(coorRightY);
    }
}


void edgeDetect(FigureCoefficient *figureCoefficientMain, FigureCoefficient *figureCoefficientNeighbors, cv::Mat &I,
                cv::Mat whiteBoard, int sizeOfNeighborhood) {
    FigureCoefficient figureCoefficientBlackNumbers;
    int colorsBlack[6] = {0, 0, 0, 20, 20, 20};
    figureCoefficientBlackNumbers.setColors(colorsBlack);
    CV_Assert(I.depth() != sizeof(uchar));
    switch (I.channels()) {
        case 3:
            cv::Mat_<cv::Vec3b> _I = I;
            cv::Mat_<cv::Vec3b> _R = whiteBoard;
            for (int x = figureCoefficientMain->getCoorMinX() + sizeOfNeighborhood;
                 x <= figureCoefficientMain->getCoorMaxX() - sizeOfNeighborhood; ++x) {
                for (int y = figureCoefficientMain->getCoorMinY() + sizeOfNeighborhood;
                     y <= figureCoefficientMain->getCoorMaxY() - sizeOfNeighborhood; ++y) {
                    int cwhite = 0;
                    int cred = 0;
                    int cblack = 0;
                    for (int k = x - sizeOfNeighborhood; k <= x + sizeOfNeighborhood; k++) {
                        for (int l = y - sizeOfNeighborhood; l <= y + sizeOfNeighborhood; l++) {
                            if (computeConditionForColor(_I, figureCoefficientMain, k, l)) {
                                cred++;
                            }
                            if (computeConditionForColor(_I, figureCoefficientNeighbors, k, l)) {
                                cwhite++;
                            }
                            if (computeConditionForColor(_I, &figureCoefficientBlackNumbers, k, l)) {
                                cblack++;
                            }
                        }
                    }
                        if (cred > 25 and (cwhite + cblack > 20)) {
                            _R(x, y)[0] = 0;
                            _R(x, y)[1] = 0;
                            _R(x, y)[2] = 0;

                        }

                }
            }
    }
}

void frameForSpeedValue(cv::Mat whiteBoard, FigureCoefficient *figureCoefficientMain) {
    FigureCoefficient figureCoefficientFrameForSpeedValue=*figureCoefficientMain;
    int colorsBlack[6] = {0, 0, 0, 10, 10, 10};

    figureCoefficientFrameForSpeedValue.setColors(colorsBlack);

    int width = figureCoefficientFrameForSpeedValue.getCoorMaxY() - figureCoefficientFrameForSpeedValue.getCoorMinY();
    int hight = figureCoefficientFrameForSpeedValue.getCoorMaxX() - figureCoefficientFrameForSpeedValue.getCoorMinX();

    figureCoefficientMain->setCoorMinY(figureCoefficientFrameForSpeedValue.getCoorMinY() + (int) width / 10);
    figureCoefficientMain->setCoorMaxY(figureCoefficientFrameForSpeedValue.getCoorMaxY() - (int) width / 10);
    figureCoefficientMain->setCoorMinX(figureCoefficientFrameForSpeedValue.getCoorMinX() + (int) hight / 6);
    figureCoefficientMain->setCoorMaxX(figureCoefficientFrameForSpeedValue.getCoorMaxX() - (int) hight / 6);
}

void findPixelsOfSpeedValues(cv::Mat whiteBoard, cv::Mat image, cv::Mat whiteBoardForSpeedValues,
                             FigureCoefficient figureCoefficientMain) {

    frameForSpeedValue(whiteBoard, &figureCoefficientMain);

    CV_Assert(image.depth() != sizeof(uchar));
    cv::Mat res(image.rows, image.cols, CV_8UC3);
    switch (image.channels()) {
        case 3:
            cv::Mat_<cv::Vec3b> _I = image;
            cv::Mat_<cv::Vec3b> _R = whiteBoardForSpeedValues;
            for (int x = figureCoefficientMain.getCoorMinX(); x <= figureCoefficientMain.getCoorMaxX(); ++x) {
                for (int y = figureCoefficientMain.getCoorMinY(); y <= figureCoefficientMain.getCoorMaxY(); ++y) {
                    if (computeConditionForColor(_I, &figureCoefficientMain, x, y)) {
                        _R(x, y)[0] = 0;
                        _R(x, y)[1] = 0;
                        _R(x, y)[2] = 0;
                    }
                }
            }
    }
}

void findFirstBlackPixel(cv::Mat whiteBoardForSpeedValues, int *firstX, int *firstY,
                         FigureCoefficient figureCoefficientMain) {
    bool isBlack = false;
    CV_Assert(whiteBoardForSpeedValues.depth() != sizeof(uchar));
    switch (whiteBoardForSpeedValues.channels()) {
        case 3:
            cv::Mat_<cv::Vec3b> _R = whiteBoardForSpeedValues;
            int x = figureCoefficientMain.getCoorMinX();
            while (x <= figureCoefficientMain.getCoorMaxX() && !isBlack) {
                int y = figureCoefficientMain.getCoorMinY();
                while (y <= figureCoefficientMain.getCoorMaxY() && !isBlack) {
                    if (computeConditionForColor(_R, &figureCoefficientMain, x, y)) {
                        *firstX = x;
                        *firstY = y;
                        isBlack = true;
                    }
                    y++;
                }
                x++;
            }
    }
}


void findStartingOfNumbers(FigureCoefficient *figureCoefficientMain, cv::Mat whiteBoard, cv::Mat image,
                           cv::Mat whiteBoardForSpeedValues, std::vector<int> &vecXsr, std::vector<int> &vecYsr) {

    findPixelsOfSpeedValues(whiteBoard, image, whiteBoardForSpeedValues, *figureCoefficientMain);

    int firstX=-1;
    int firstY=-1;
    findFirstBlackPixel(whiteBoardForSpeedValues, &firstX, &firstY, *figureCoefficientMain);
    if (firstX<0) {
        return;
    }
    std::vector<std::vector<int>> vecX;
    std::vector<std::vector<int>> vecY;
    std::vector<int> vecXs;
    std::vector<int> vecYs;
    vecXs.push_back(firstX);
    vecYs.push_back(firstY);

    bool isBlack = true;

    CV_Assert(whiteBoardForSpeedValues.depth() != sizeof(uchar));
    switch (whiteBoardForSpeedValues.channels()) {
        case 3:
            cv::Mat_<cv::Vec3b> _wB = whiteBoardForSpeedValues;
            for (int x = firstX; x < firstX + 14; x++) {

                for (int y = figureCoefficientMain->getCoorMinY(); y <= figureCoefficientMain->getCoorMaxY(); y++) {
                    if (computeConditionForColor(_wB, figureCoefficientMain, x, y) &&
                        !isBlack) {
                        if (vecXs.size() == 0) {
                            vecXs.push_back(x);
                            vecYs.push_back(y);
                            isBlack = true;
                        }
                        if ((vecXs.size() > 0) && (y - vecYs.back() > 3)) {
                            vecXs.push_back(x);
                            vecYs.push_back(y);
                            isBlack = true;
                        }
                    }
                    if (!computeConditionForColor(_wB, figureCoefficientMain, x, y)) {
                        isBlack = false;
                    }
                }
                vecX.push_back(vecXs);
                vecY.push_back(vecYs);
                vecXs.clear();
                vecYs.clear();
            }
    }
    vecXsr.insert(vecXsr.begin(), vecX[0].begin(), vecX[0].end());
    vecYsr.insert(vecYsr.begin(), vecY[0].begin(), vecY[0].end());

    std::cout << "vec 1: " << vecY[0].size() << std::endl;
    std::cout << "vec 2: " << vecY[1].size() << std::endl;
    std::cout << "vec 3: " << vecY[2].size() << std::endl;
    std::cout << "vec 4: " << vecY[3].size() << std::endl;
    for (int k = 1; k < 14; k++) {
        if (vecX[k].size() > vecXsr.size()) {
            vecXsr.clear();
            vecYsr.clear();
            vecXsr.insert(vecXsr.begin(), vecX[k].begin(), vecX[k].end());
            vecYsr.insert(vecYsr.begin(), vecY[k].begin(), vecY[k].end());
        }
    }
}



void findNeighborhood(cv::Mat whiteBoardForSpeedValues, int firstX, int firstY,
                      FigureCoefficient figureCoefficientMain, std::vector<cv::Point> &checkedPixels, int distance) {

    std::vector<cv::Point> pixelsToCheck;
    cv::Point pixelFirst(firstX, firstY);
    pixelsToCheck.push_back(pixelFirst);

    CV_Assert(whiteBoardForSpeedValues.depth() != sizeof(uchar));
    switch (whiteBoardForSpeedValues.channels()) {
        case 3:
            cv::Mat_<cv::Vec3b> _R = whiteBoardForSpeedValues;
            while (!pixelsToCheck.empty()) {
                cv::Point pixel = pixelsToCheck.back();
                if (!std::count(checkedPixels.begin(), checkedPixels.end(), pixel)) {
                    checkedPixels.push_back(pixel);
                }
                pixelsToCheck.pop_back();
                for (int x = pixel.x - distance; x <= pixel.x + distance; x++) {
                    for (int y = pixel.y - distance; y <= pixel.y + distance; y++) {
                        if (computeConditionForColor(_R, &figureCoefficientMain, x, y)) {
                            cv::Point pixelNext(x, y);
                            if (!std::count(checkedPixels.begin(), checkedPixels.end(), pixelNext)) {
                                if (!std::count(pixelsToCheck.begin(), pixelsToCheck.end(), pixelNext)) {
                                    pixelsToCheck.push_back(pixelNext);
                                }
                            }
                        }
                    }
                }
            }
    }
}


std::vector<std::vector<cv::Point>> findCircles(FigureCoefficient *figureCoefficientMain,
                                                cv::Mat whiteBoardForSpeedValues) {
    std::vector<std::vector<cv::Point>> ret;
    int firstX=-1;
    int firstY=-1;
    int distance=9;
    figureCoefficientMain->setCoorMaxY(whiteBoardForSpeedValues.cols-(distance+1));
    figureCoefficientMain->setCoorMaxX(whiteBoardForSpeedValues.rows-(distance+1));
    figureCoefficientMain->setCoorMinY((distance+1));
    figureCoefficientMain->setCoorMinX((distance+1));
    findFirstBlackPixel(whiteBoardForSpeedValues, &firstX, &firstY, *figureCoefficientMain);

    while (firstX>=0) {
        std::vector<cv::Point> points;
        findNeighborhood(whiteBoardForSpeedValues, firstX, firstY, *figureCoefficientMain,
                         points,distance);
        cv::Mat_<cv::Vec3b> _I = whiteBoardForSpeedValues;
        for (auto p:points) {
            _I(p.x,p.y)[0]=255;
            _I(p.x,p.y)[1]=255;
            _I(p.x,p.y)[2]=255;
        }
        if (points.size()>20) {
            ret.push_back(points);
        }
        firstX=-1;
        firstY=-1;
        findFirstBlackPixel(whiteBoardForSpeedValues, &firstX, &firstY, *figureCoefficientMain);
    }
    return ret;
}



void computeField(FigureCoefficient *figureCoefficient, cv::Mat &I) {
    int field = 0;
    CV_Assert(I.depth() != sizeof(uchar));
    switch (I.channels()) {
        case 3:
            cv::Mat_<cv::Vec3b> _I = I;
            for (int x = figureCoefficient->getCoorMinX(); x <= figureCoefficient->getCoorMaxX(); ++x) {
                for (int y = figureCoefficient->getCoorMinY(); y <= figureCoefficient->getCoorMaxY(); ++y) {
                    if (computeConditionForColor(_I, figureCoefficient, x, y)) {
                        field = field + 1;
                    }
                }
            }
    }
    figureCoefficient->setField(field);
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
            for (int x = figureCoefficient->getCoorMinX(); x <= figureCoefficient->getCoorMaxX(); ++x) {
                for (int y = figureCoefficient->getCoorMinY(); y <= figureCoefficient->getCoorMaxY(); ++y) {
                    if (_I(x, y)[0] <= b && _I(x, y)[1] <= g && _I(x, y)[2] <= r) {
                        if ((_I(x - 1, y)[0] > b && _I(x - 1, y)[1] > g && _I(x - 1, y)[2] > r) ||
                            (_I(x, y - 1)[0] > b && _I(x, y - 1)[1] > g && _I(x, y - 1)[2] > r) ||
                            (_I(x + 1, y)[0] > b && _I(x + 1, y)[1] > g && _I(x + 1, y)[2] > r) ||
                            (_I(x, y + 1)[0] > b && _I(x, y + 1)[1] > g && _I(x, y + 1)[2] > r)) {
                            circumference = circumference + 1;
                        }
                    }
                }
            }
    }
    figureCoefficient->setCircumference(circumference);
}

void computeCoefficientOfMalinowska(FigureCoefficient *figureCoefficient) {
    float coefficientOfMalinowska =
            (figureCoefficient->getCircumference() / (2 * std::sqrtf(M_PI * figureCoefficient->getField()))) - 1;
    figureCoefficient->setW3(coefficientOfMalinowska);
    std::cout << figureCoefficient->getW3();
}

void computeMomentsToCheckNumbers(FigureCoefficient *figureCoefficient, cv::Mat &I) {

    double m00 = 0;
    double m01 = 0;
    double m10 = 0;
    double m11 = 0;
    double m02 = 0;
    double m20 = 0;
    double m12 = 0;
    double m21 = 0;
    double m03 = 0;
    double m30 = 0;

    int b = figureCoefficient->getColors()[0];
    int g = figureCoefficient->getColors()[1];
    int r = figureCoefficient->getColors()[2];

    CV_Assert(I.depth() != sizeof(uchar));
    cv::Mat res(I.rows, I.cols, CV_8UC3);
    switch (I.channels()) {
        case 3:
            cv::Mat_<cv::Vec3b> _I = I;
            for (int x = figureCoefficient->getCoorMinX(); x <= figureCoefficient->getCoorMaxX(); ++x) {
                for (int y = figureCoefficient->getCoorMinY(); y <= figureCoefficient->getCoorMaxY(); ++y) {
                    if (_I(x, y)[0] == b && _I(x, y)[1] == g && _I(x, y)[2] == r) {
                        m00 = m00 + pow(x, 0) * pow(y, 0);
                        m01 = m01 + pow(x, 0) * pow(y, 1);
                        m10 = m10 + pow(x, 1) * pow(y, 0);
                        m11 = m11 + pow(x, 1) * pow(y, 1);
                        m02 = m02 + pow(x, 0) * pow(y, 2);
                        m20 = m20 + pow(x, 2) * pow(y, 0);
                        m12 = m12 + pow(x, 1) * pow(y, 2);
                        m21 = m21 + pow(x, 2) * pow(y, 1);
                        m03 = m03 + pow(x, 0) * pow(y, 3);
                        m30 = m30 + pow(x, 3) * pow(y, 0);
                    }
                }
            }
    }

    float i = m10 / m00;
    float j = m01 / m00;
    float M20 = m20 - (pow(m10, 2)) / m00;
    float M02 = m02 - (pow(m01, 2)) / m00;
    float M11 = m11 - m10 * m01 / m00;
    float M12 = m12 - 2 * m11 * j - m02 * i + 2 * m10 * pow(j, 2);
    float M21 = m21 - 2 * m11 * i - m20 * j + 2 * m01 * pow(i, 2);
    float M30 = m30 - 3 * m20 * i + 2 * m10 * pow(i, 2);
    float M03 = m03 - 3 * m02 * j + 2 * m01 * pow(j, 2);

    float M1 = (M20 + M02) / pow(m00, 2);
    float M2 = (pow((M20 - M02), 2) + 4 * pow(M11, 2)) / pow(m00, 4);
    float M3 = (pow((M30 - 3 * M12), 2) + pow((3 * M21 - M03), 2)) / pow(m00, 5);
    float M4 = (pow((M30 + M12), 2) + pow((M21 + M03), 2)) / pow(m00, 5);
    float M7 = (M20 * M02 - pow(M11, 2)) / pow(m00, 4);
    float M8 = (M30 * M12 + M21 * M03 - pow(M12, 2) - pow(M21, 2)) / pow(m00, 5);
    float M10 = ((pow((M30 * M03 - M12 * M21), 2)) - 4 * (M30 * M12 - pow(M21, 2)) * (M03 * M21 - M12)) / pow(m00, 10);

    figureCoefficient->setM1(M1);
    figureCoefficient->setM2(M2);
    figureCoefficient->setM3(M3);
    figureCoefficient->setM4(M4);
    figureCoefficient->setM7(M7);
    figureCoefficient->setM8(M8);
    figureCoefficient->setM10(M10);
}

void computeMoments(FigureCoefficient *figureCoefficient, std::vector<cv::Point> checkedPixels) {

    double m00 = 0;
    double m01 = 0;
    double m10 = 0;
    double m11 = 0;
    double m02 = 0;
    double m20 = 0;
    double m12 = 0;
    double m21 = 0;
    double m03 = 0;
    double m30 = 0;

    int b = figureCoefficient->getColors()[0];
    int g = figureCoefficient->getColors()[1];
    int r = figureCoefficient->getColors()[2];

    for (auto i = checkedPixels.cbegin(); i != checkedPixels.cend(); ++i) {
        int x = i->x;
        int y = i->y;

        m00 = m00 + pow(x, 0) * pow(y, 0);
        m01 = m01 + pow(x, 0) * pow(y, 1);
        m10 = m10 + pow(x, 1) * pow(y, 0);
        m11 = m11 + pow(x, 1) * pow(y, 1);
        m02 = m02 + pow(x, 0) * pow(y, 2);
        m20 = m20 + pow(x, 2) * pow(y, 0);
        m12 = m12 + pow(x, 1) * pow(y, 2);
        m21 = m21 + pow(x, 2) * pow(y, 1);
        m03 = m03 + pow(x, 0) * pow(y, 3);
        m30 = m30 + pow(x, 3) * pow(y, 0);
    }


    float i = m10 / m00;
    float j = m01 / m00;
    float M20 = m20 - (pow(m10, 2)) / m00;
    float M02 = m02 - (pow(m01, 2)) / m00;
    float M11 = m11 - m10 * m01 / m00;
    float M12 = m12 - 2 * m11 * j - m02 * i + 2 * m10 * pow(j, 2);
    float M21 = m21 - 2 * m11 * i - m20 * j + 2 * m01 * pow(i, 2);
    float M30 = m30 - 3 * m20 * i + 2 * m10 * pow(i, 2);
    float M03 = m03 - 3 * m02 * j + 2 * m01 * pow(j, 2);

    float M1 = (M20 + M02) / pow(m00, 2);
    float M2 = (pow((M20 - M02), 2) + 4 * pow(M11, 2)) / pow(m00, 4);
    float M3 = (pow((M30 - 3 * M12), 2) + pow((3 * M21 - M03), 2)) / pow(m00, 5);
    float M4 = (pow((M30 + M12), 2) + pow((M21 + M03), 2)) / pow(m00, 5);
    float M7 = (M20 * M02 - pow(M11, 2)) / pow(m00, 4);
    float M8 = (M30 * M12 + M21 * M03 - pow(M12, 2) - pow(M21, 2)) / pow(m00, 5);
    float M10 = ((pow((M30 * M03 - M12 * M21), 2)) - 4 * (M30 * M12 - pow(M21, 2)) * (M03 * M21 - M12)) / pow(m00, 10);

    figureCoefficient->setM1(M1);
    figureCoefficient->setM2(M2);
    figureCoefficient->setM3(M3);
    figureCoefficient->setM4(M4);
    figureCoefficient->setM7(M7);
    figureCoefficient->setM8(M8);
    figureCoefficient->setM10(M10);
}


cv::Mat changeColorInRedCircle(cv::Mat image) {
    CV_Assert(image.depth() != sizeof(uchar));
    cv::Mat res(image.rows, image.cols, CV_8UC3);
    switch (image.channels()) {
        case 3:
            cv::Mat_<cv::Vec3b> _R = res;
            cv::Mat_<cv::Vec3b> _I = image;
            for (int x = 0; x < image.rows; ++x)
                for (int y = 0; y < image.cols; ++y) {
                    if (_I(x, y)[0] < 105 && _I(x, y)[0] > 85 && _I(x, y)[1] < 130 && _I(x, y)[1] > 80 && _I(x, y)[2] > 215) {
                        _R(x, y)[0] = 92;
                        _R(x, y)[1] = 81;
                        _R(x, y)[2] = _I(x, y)[2];
                    } else {
                        _R(x, y)[0] = _I(x, y)[0];
                        _R(x, y)[1] = _I(x, y)[1];
                        _R(x, y)[2] = _I(x, y)[2];
                    }
                }
            res = _R;
            break;
    }
    return res;
}




void displayResult(FigureCoefficient figureCoefficientMain, std::vector<std::vector<cv::Point>> points, cv::Mat &I, int colorChange) {
    std::cout << "Groups num: " << points.size() << "\n";

    if (points.size() == 3) {
        computeMoments(&figureCoefficientMain, points[1]);
    }
    else if (points.size()>0) {
        computeMoments(&figureCoefficientMain, points[0]);
    }
    bool one = figureCoefficientMain.getM7() < 0.019;
    bool two = figureCoefficientMain.getM2() >= 0.14 && figureCoefficientMain.getM2() < 0.28;
    bool four = figureCoefficientMain.getM10() < 0 && figureCoefficientMain.getM8() > -0.005 && figureCoefficientMain.getM3() > 0.01;
    bool eight = figureCoefficientMain.getM8() > -0.001 && figureCoefficientMain.getM1() > 0.3 && figureCoefficientMain.getM10() > -0.001
                 && figureCoefficientMain.getM3() < 0.01 && figureCoefficientMain.getM2() > 0.05 && figureCoefficientMain.getM7() > 0.025;
    bool zero = figureCoefficientMain.getM3() < 0.003 && figureCoefficientMain.getM2() < 0.05;
    if (points.size() == 3) {
        std::cout << "1";
        if (one) {
            std::cout << "1";
        }
        if (two) {
            std::cout << "2";
        }
        if (four) {
            std::cout << "4";
        }
        if (eight) {
            std::cout << "8";
        }
        if (zero) {
            std::cout << "0";
        }
        std::cout << "0";
        cv::Point p1(figureCoefficientMain.getCoorMinY(), figureCoefficientMain.getCoorMinX());
        cv::Point p2(figureCoefficientMain.getCoorMaxY(), figureCoefficientMain.getCoorMaxX());
        cv::rectangle(I, p1, p2,
                  cv::Scalar(255, 0, 0 + colorChange),
                  2, cv::LINE_8);
    } else if (points.size() == 2) {
        if (one) {
            std::cout << "1";
        }
        if (two) {
            std::cout << "2";
        }
        if (four) {
            std::cout << "4";
        }
        if (eight) {
            std::cout << "8";
        }
        std::cout << "0";
        cv::Point p1(figureCoefficientMain.getCoorMinY(), figureCoefficientMain.getCoorMinX());
        cv::Point p2(figureCoefficientMain.getCoorMaxY(), figureCoefficientMain.getCoorMaxX());
        cv::rectangle(I, p1, p2,
                      cv::Scalar(255, 0, 0 + colorChange),
                      2, cv::LINE_8);
    } else if (points.size() == 1) {
        if (one) {
            std::cout << "1";
        }
        if (two) {
            std::cout << "2";
        }
        if (four) {
            std::cout << "4";
        }
        if (eight) {
            std::cout << "8";
        }
        cv::Point p1(figureCoefficientMain.getCoorMinY(), figureCoefficientMain.getCoorMinX());
        cv::Point p2(figureCoefficientMain.getCoorMaxY(), figureCoefficientMain.getCoorMaxX());
        cv::rectangle(I, p1, p2,
                      cv::Scalar(255, 0, 0 + colorChange),
                      2, cv::LINE_8);
    } else {
        std::cout << "Tu nie ma znaku ograniczenia prędkości. ";
    }

}

void mergingAndCutToSmallGroups( std::vector<std::vector<cv::Point>> &points, cv::Mat whiteBoardForSpeedValues,
                                 FigureCoefficient figureCoefficientMain, std::vector<int> &vecXs, std::vector<int> &vecYs, int minGroupNumber) {
    std::vector<std::vector<long>> pointsHash;
    for (int j = 0; j < vecXs.size(); ++j) {
        int firstX = vecXs[j];
        int firstY = vecYs[j];
        std::vector<cv::Point> checkedPixels;
        std::vector<long> hash;
        findNeighborhood(whiteBoardForSpeedValues, firstX, firstY,
                         figureCoefficientMain, checkedPixels,1);
        for (auto p:checkedPixels) {
            hash.push_back(p.x * 512 + p.y);
        }

        std::sort(hash.begin(), hash.end());
        if (checkedPixels.size() > minGroupNumber) {
            points.push_back(checkedPixels);
            pointsHash.push_back(hash);
        }
    }
    if (points.size()<1){
        return;
    }
    int k = 0;
    while (k < points.size() - 1) {
        int l = k + 1;
        while (l < points.size()) {
            std::vector<long> intersect;
            std::set_intersection(pointsHash[k].begin(), pointsHash[k].end(), pointsHash[l].begin(),
                                  pointsHash[l].end(), std::back_inserter(intersect));
            if (intersect.size() > 0) {
                points.erase(points.begin() + l);
                pointsHash.erase(pointsHash.begin() + l);
            } else {
                l = l + 1;
            }
        }
        k = k + 1;
    }
}


void findSpeedLimitSign() {
    cv::Mat speedLimitSign1 = cv::imread("SourceImages/road112_AgnieszkaJurkiewicz.png");
    cv::Mat speedLimitSign2 = cv::imread("SourceImages/road113_AgnieszkaJurkiewicz.png");
    cv::Mat speedLimitSign3 = cv::imread("SourceImages/road119_AgnieszkaJurkiewicz.png");
    cv::Mat speedLimitSign4 = cv::imread("SourceImages/road734_AgnieszkaJurkiewicz_dwaZnaki.png");
    cv::Mat imagesList[4] = {speedLimitSign1, speedLimitSign2, speedLimitSign3, speedLimitSign4};

    int i = 1;

    for (cv::Mat image : imagesList) {
        int colorChange = 10;
        cv::Mat whiteBoardWithCircle = makeWhiteBoard(image);
        cv::Mat clearWhiteBoard = makeWhiteBoard(image);

        FigureCoefficient figureCoefficientRedCircle;
        int colorsRed[6] = {0, 0, 70, 93, 82, 255};
        figureCoefficientRedCircle.setColors(colorsRed);

        FigureCoefficient figureCoefficientWhiteInnerCircle;
        int colorsWhite[6] = {155, 155, 155, 255, 255, 255};
        figureCoefficientWhiteInnerCircle.setColors(colorsWhite);

        FigureCoefficient figureCoefficientMain;
        int colorsBlack[6] = {0, 0, 0, 130, 130, 130};
        figureCoefficientMain.setColors(colorsBlack);

        cv::Mat imageMoreRed = changeColorInRedCircle(image);

        computeBox(imageMoreRed, &figureCoefficientRedCircle);

        edgeDetect(&figureCoefficientRedCircle, &figureCoefficientWhiteInnerCircle, imageMoreRed, whiteBoardWithCircle, 5);

        std::vector<std::vector<cv::Point>> circles=findCircles(&figureCoefficientMain, whiteBoardWithCircle);
        for (std::vector<cv::Point> circle: circles) {
            int minx=10000,miny=10000;
            int maxx=-11,maxy=-1;
            for (cv::Point p:circle) {
                if (minx>p.x) {
                    minx=p.x;
                }
                if (miny>p.y) {
                    miny=p.y;
                }
                if (maxx<p.x) {
                    maxx=p.x;
                }
                if (maxy<p.y) {
                    maxy=p.y;
                }
            }

            std::vector<int> vecXsCircle;
            std::vector<int> vecYsCircle;
            figureCoefficientMain.setCoorMaxY(maxy);
            figureCoefficientMain.setCoorMinY(miny);
            figureCoefficientMain.setCoorMaxX(maxx);
            figureCoefficientMain.setCoorMinX(minx);

            cv::Mat whiteBoardForSpeedValues = makeWhiteBoard(image);
            findStartingOfNumbers(&figureCoefficientMain, whiteBoardWithCircle, image, whiteBoardForSpeedValues, vecXsCircle,
                                  vecYsCircle);
            std::vector<std::vector<cv::Point>> groupOfPointsInCircle;
            int minGroupNumber = 20;
            mergingAndCutToSmallGroups(groupOfPointsInCircle, whiteBoardForSpeedValues,
                                       figureCoefficientMain, vecXsCircle, vecYsCircle, minGroupNumber);

            displayResult(figureCoefficientMain, groupOfPointsInCircle, image, colorChange);
            colorChange = colorChange * 20;

                std::cout << "\ncircle"<<i<<"\n";

                //cv::Mat image2 = speedLimitSign1(cv::Rect(figureCoefficientRedCircle.getPixels()[0], figureCoefficientRedCircle.getPixels()[1],
                //                                          figureCoefficientRedCircle.getPixels()[2] - figureCoefficientRedCircle.getPixels()[0], figureCoefficientRedCircle.getPixels()[3]
                //                                                                                                                                 - figureCoefficientRedCircle.getPixels()[1]));
                cv::imshow("Image" + std::to_string(i), image);
                //cv::imshow("Shape" + std::to_string(i), whiteBoardWithCircle);
                //cv::imshow("Values" + std::to_string(i), whiteBoardForSpeedValues);
                cv::waitKey(-1);
        }
    }
}


void countCoefficientFoNumbers() {
    //int amountOfNumbers = 12;
    //std::string fileNames[12] = {"0a", "0b", "0c", "0d", "0e", "0f", "0g", "0h", "0i", "0j", "0k", "0l"};
    //int amountOfNumbers = 3;
    //std::string fileNames[3] = {"1a", "1b", "1c"};
    //int amountOfNumbers = 8;
    //std::string fileNames[8] = {"2a", "2b", "2c", "2d", "2e", "2f", "2g", "2h"};
    int amountOfNumbers = 4;
    std::string fileNames[4] = {"4a", "4b", "4c", "4d"};
    //int amountOfNumbers = 4;
    //std::string fileNames[4] = {"8a", "8b", "8c", "8d"};

    for (int i = 0; i < amountOfNumbers; i++) {
        FigureCoefficient figureCoefficient;
        int colorsBlack[6] = {0, 0, 0, 100, 100, 100};
        figureCoefficient.setColors(colorsBlack);
        figureCoefficient.setFileName(fileNames[i]);
        cv::Mat image = cv::imread("SourceImages/Numbers/" + fileNames[i] + ".png");

        computeBox(image, &figureCoefficient);
        computeField(&figureCoefficient, image);
        computeCircumference(&figureCoefficient, image);
        computeCoefficientOfMalinowska(&figureCoefficient);

        computeMomentsToCheckNumbers(&figureCoefficient, image);
        std::cout << ", " << figureCoefficient.getM1() << ", " << figureCoefficient.getM2() << ", "
                  << figureCoefficient.getM3()
                  << ", " << figureCoefficient.getM4() << ", " << figureCoefficient.getM7() << ", "
                  << figureCoefficient.getM8() << ", " << figureCoefficient.getM10()
                  << std::endl;
        cv::imshow("Shape", image);
        cv::waitKey(-1);
    }
}


int main(int, char *[]) {
    std::cout << "Start ..." << std::endl;

    findSpeedLimitSign();


    cv::waitKey(-1);


    //countCoefficientFoNumbers();
    return 0;
}