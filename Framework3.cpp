#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include <iostream>
#include <math.h>
#include "FigureCoefficient.h"
#include <vector>
#include "Pixel.h"

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
    cv::Mat res(I.rows, I.cols, CV_8UC3);
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
            //std::cout << "coorDownX = " << figureCoefficient->getCoorMinX() << ", coorUpX = " << figureCoefficient->getCoorMaxX()
            //        << ", coorLeftY = " << figureCoefficient->getCoorMinY() << ", coorRightY = " << figureCoefficient->getCoorMaxY()<< std::endl;
            //std::cout << "Width = " << figureCoefficient->getCoorMaxY() - figureCoefficient->getCoorMinY() << "Hight = " << figureCoefficient->getCoorMaxX()
            //                                                                     - figureCoefficient->getCoorMinX() << std::endl;
    }
}


void findNeighborhoodNotDiagonally(FigureCoefficient *figureCoefficientWhiteBoard, cv::Mat whiteBoard,
                                   int sizeOfNeighborhood, int r, int g, int b) {
    CV_Assert(whiteBoard.depth() != sizeof(uchar));
    cv::Mat res(whiteBoard.rows, whiteBoard.cols, CV_8UC3);
    switch (whiteBoard.channels()) {
        case 3:
            cv::Mat_<cv::Vec3b> _R = whiteBoard;
            for (int x = figureCoefficientWhiteBoard->getCoorMinX() + sizeOfNeighborhood;
                 x <= figureCoefficientWhiteBoard->getCoorMaxX() - sizeOfNeighborhood; x++) {
                for (int y = figureCoefficientWhiteBoard->getCoorMinY() + sizeOfNeighborhood;
                     y <= figureCoefficientWhiteBoard->getCoorMaxY() - sizeOfNeighborhood; y++) {
                    if (computeConditionForColor(_R, figureCoefficientWhiteBoard, x, y)) {
                        int numberOfNeighbors = 0;
                        if (computeConditionForColor(_R, figureCoefficientWhiteBoard, x - 1, y)) {
                            numberOfNeighbors++;
                        }
                        if (computeConditionForColor(_R, figureCoefficientWhiteBoard, x, y - 1)) {
                            numberOfNeighbors++;
                        }
                        if (computeConditionForColor(_R, figureCoefficientWhiteBoard, x + 1, y)) {
                            numberOfNeighbors++;
                        }
                        if (computeConditionForColor(_R, figureCoefficientWhiteBoard, x, y + 1)) {
                            numberOfNeighbors++;
                        }
                        if (numberOfNeighbors <= 2) {
                            _R(y, x)[0] = b;
                            _R(y, x)[1] = g;
                            _R(y, x)[2] = r;
                        }
                    }
                }
            }
    }
}


void edgeDetect(FigureCoefficient *figureCoefficientMain, FigureCoefficient *figureCoefficientNeighbors, cv::Mat &I,
                cv::Mat whiteBoard, int sizeOfNeighborhood) {
    FigureCoefficient figureCoefficientBlackNumbers;
    int colorsBlack[6] = {0, 0, 0, 20, 20, 20};
    figureCoefficientBlackNumbers.setColors(colorsBlack);
    CV_Assert(I.depth() != sizeof(uchar));
    cv::Mat res(I.rows, I.cols, CV_8UC3);
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
                        // if (cwhite>sizeOfNeighborhood and cred>sizeOfNeighborhood*sizeOfNeighborhood*2) {
                        if (cred > 25 and (cwhite > 20 or cwhite + cblack > 20)) {
                            _R(x, y)[0] = 0;
                            _R(x, y)[1] = 0;
                            _R(x, y)[2] = 0;
                        }
                    }

                }
            }
    }
}

void frameForSpeedValue(cv::Mat whiteBoard, FigureCoefficient *figureCoefficientMain) {
    FigureCoefficient figureCoefficientFrameForSpeedValue;
    int colorsBlack[6] = {0, 0, 0, 10, 10, 10};
    figureCoefficientFrameForSpeedValue.setColors(colorsBlack);

    computeBox(whiteBoard, &figureCoefficientFrameForSpeedValue);

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
    //cv::Mat res(whiteBoardForSpeedValues.rows, whiteBoardForSpeedValues.cols, CV_8UC3);
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

void findStartingPointsForAllNumbers(cv::Mat whiteBoardForSpeedValues, int firstX, int firstY,
                                     FigureCoefficient figureCoefficientMain,
                                     std::vector<std::vector<int>> vecX, std::vector<std::vector<int>> vecY) {
    bool isBlack = true;
    std::vector<int> vecXs;
    std::vector<int> vecYs;
    CV_Assert(whiteBoardForSpeedValues.depth() != sizeof(uchar));
    //cv::Mat res(I.rows, I.cols, CV_8UC3);
    switch (whiteBoardForSpeedValues.channels()) {
        case 3:
            cv::Mat_<cv::Vec3b> _wB = whiteBoardForSpeedValues;
            for (int x = firstX; x < firstX + 4; x++) {
                vecXs.clear();
                vecYs.clear();
                for (int y = firstY; y <= figureCoefficientMain.getCoorMaxY(); y++) {
                    if (computeConditionForColor(_wB, &figureCoefficientMain, x, y) &&
                        !isBlack) {
                        vecXs.push_back(x);
                        vecYs.push_back(y);
                        isBlack = true;
                    }
                    if (!computeConditionForColor(_wB, &figureCoefficientMain, x, y)) {
                        isBlack = false;
                    }
                }
                vecX.push_back(vecXs);
                vecY.push_back(vecYs);
            }
    }
}

void findStartingOfNumbers(FigureCoefficient *figureCoefficientMain, cv::Mat whiteBoard, cv::Mat image,
                           cv::Mat whiteBoardForSpeedValues, std::vector<int> &vecXsr, std::vector<int> &vecYsr) {

    findPixelsOfSpeedValues(whiteBoard, image, whiteBoardForSpeedValues, *figureCoefficientMain);

    int firstX;
    int firstY;
    findFirstBlackPixel(whiteBoardForSpeedValues, &firstX, &firstY, *figureCoefficientMain);

    std::vector<std::vector<int>> vecX;
    std::vector<std::vector<int>> vecY;
    std::vector<int> vecXs;
    std::vector<int> vecYs;
    vecXs.push_back(firstX);
    vecYs.push_back(firstY);


    //findStartingPointsForAllNumbers(whiteBoardForSpeedValues, firstX, firstY, figureCoefficientMain, vecX, vecY);

    bool isBlack = true;

    CV_Assert(whiteBoardForSpeedValues.depth() != sizeof(uchar));
    //cv::Mat res(I.rows, I.cols, CV_8UC3);
    switch (whiteBoardForSpeedValues.channels()) {
        case 3:
            cv::Mat_<cv::Vec3b> _wB = whiteBoardForSpeedValues;
            for (int x = firstX; x < firstX + 4; x++) {

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
    std::cout << "vec 2: " << vecY[2].size() << std::endl;
    std::cout << "vec 2: " << vecY[3].size() << std::endl;
    for (int k = 1; k < 4; k++) {
        if (vecX[k].size() > vecXsr.size()) {
            vecXsr.clear();
            vecYsr.clear();
            vecXsr.insert(vecXsr.begin(), vecX[k].begin(), vecX[k].end());
            vecYsr.insert(vecYsr.begin(), vecY[k].begin(), vecY[k].end());
        }
    }
}


void findNeighborhood(cv::Mat whiteBoard, cv::Mat image, cv::Mat whiteBoardForSpeedValues, int firstX, int firstY,
                      FigureCoefficient figureCoefficientMain, std::vector<cv::Point> &checkedPixels) {

    std::vector<cv::Point> pixelsToCheck;
    cv::Point pixelFirst(firstX, firstY);
    pixelsToCheck.push_back(pixelFirst);

    CV_Assert(whiteBoardForSpeedValues.depth() != sizeof(uchar));
    //cv::Mat res(whiteBoardForSpeedValues.rows, whiteBoardForSpeedValues.cols, CV_8UC3);
    switch (whiteBoardForSpeedValues.channels()) {
        case 3:
            cv::Mat_<cv::Vec3b> _R = whiteBoardForSpeedValues;
            while (!pixelsToCheck.empty()) {
                cv::Point pixel = pixelsToCheck.back();
                if (!std::count(checkedPixels.begin(), checkedPixels.end(), pixel)) {
                    checkedPixels.push_back(pixel);
                }
                pixelsToCheck.pop_back();
                for (int x = pixel.x - 1; x <= pixel.x + 1; x++) {
                    for (int y = pixel.y - 1; y <= pixel.y + 1; y++) {
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


void computeField(FigureCoefficient *figureCoefficient, cv::Mat &I) {
    int field = 0;
    //int b = figureCoefficient->getColors()[0];
    //int g = figureCoefficient->getColors()[1];
    //int r = figureCoefficient->getColors()[2];
    CV_Assert(I.depth() != sizeof(uchar));
    cv::Mat res(I.rows, I.cols, CV_8UC3);
    switch (I.channels()) {
        case 3:
            cv::Mat_<cv::Vec3b> _I = I;
            for (int x = figureCoefficient->getCoorMinX(); x <= figureCoefficient->getCoorMaxX(); ++x) {
                for (int y = figureCoefficient->getCoorMinY(); y <= figureCoefficient->getCoorMaxY(); ++y) {
                    if (computeConditionForColor(_I, figureCoefficient, x, y)) {
                        field = field + 1;
                        //_I(x, y)[0] = 255;
                        //_I(x, y)[1] = 255;
                        //_I(x, y)[2] = 255;
                    }
                }
            }
    }
    figureCoefficient->setField(field);
    //std::cout <<  ": S = " << figureCoefficient->getField();
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
                            //_I(i, j)[1] = 255;
                        }
                    }
                }
            }
    }
    figureCoefficient->setCircumference(circumference);
    //std::cout << ", L = " << figureCoefficient->getCircumference();
}

void computeCoefficientOfMalinowska(FigureCoefficient *figureCoefficient) {
    float coefficientOfMalinowska =
            (figureCoefficient->getCircumference() / (2 * std::sqrtf(M_PI * figureCoefficient->getField()))) - 1;
    figureCoefficient->setW3(coefficientOfMalinowska);
    //std::cout << ", W3 = " << figureCoefficient->getW3();
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

        //cv::Mat image2 = image(cv::Rect(figureCoefficient.getCoorMinY(), figureCoefficient.getCoorMinX(),
        //                                figureCoefficient.getCoorMaxY() - figureCoefficient.getCoorMinY(),
        //                                figureCoefficient.getCoorMaxX()
        //                                - figureCoefficient.getCoorMinX()));
        //std::cout << i + 1 << ". Plik " + fileNames[i];

        computeField(&figureCoefficient, image);

        computeCircumference(&figureCoefficient, image);
        computeCoefficientOfMalinowska(&figureCoefficient);

        computeMomentsToCheckNumbers(&figureCoefficient, image);
        //std::cout << ", M1 = " << figureCoefficient.getM1() << ", M7 = " << figureCoefficient.getM7() << std::endl;
        std::cout << ", " << figureCoefficient.getM1() << ", " << figureCoefficient.getM2() << ", "
                  << figureCoefficient.getM3()
                  << ", " << figureCoefficient.getM4() << ", " << figureCoefficient.getM7() << ", "
                  << figureCoefficient.getM8() << ", " << figureCoefficient.getM10()
                  << std::endl;
        cv::imshow("Shape", image);
        cv::waitKey(-1);
    }
}

bool cmp(const cv::Point &a, const cv::Point &b) {
    if (a.x < b.x) {
        return true;
    }
    if (a.x > b.x) {
        return false;
    }
    if (a.y < b.y) {
        return true;
    }
    return false;
}

void findSpeedLimitSign() {
    cv::Mat speedLimitSign1 = cv::imread("SourceImages/road112_AgnieszkaJurkiewicz.png");
    cv::Mat speedLimitSign2 = cv::imread("SourceImages/road113_AgnieszkaJurkiewicz.png");
    cv::Mat speedLimitSign3 = cv::imread("SourceImages/road119_AgnieszkaJurkiewicz.png");
    cv::Mat imagesList[3] = {speedLimitSign1, speedLimitSign2, speedLimitSign3};

    int i = 1;
    for (cv::Mat image : imagesList) {
        cv::Mat whiteBoardWithCircle = makeWhiteBoard(image);
        cv::Mat clearWhiteBoard = makeWhiteBoard(image);

        FigureCoefficient figureCoefficientRedCircle;
        int colorsRed[6] = {0, 0, 70, 93, 82, 255}; //czerwony
        //int colorsRed [6] = {110, 100, 150, 170, 160, 200}; // różowy
        //int colorsRed [6] = {0, 0, 40, 20, 20, 250};
        figureCoefficientRedCircle.setColors(colorsRed);

        FigureCoefficient figureCoefficientWhiteInnerCircle;
        //int colorsWhite [6] = {130, 120, 155, 255, 255, 255};
        int colorsWhite[6] = {155, 155, 155, 255, 255, 255};
        figureCoefficientWhiteInnerCircle.setColors(colorsWhite);

        computeBox(image, &figureCoefficientRedCircle);

        //findNeighborhood(&figureCoefficientRedCircle, &figureCoefficientWhiteInnerCircle, speedLimitSign1, whiteBoardWithCircle, 4, 0, 0, 0);
        edgeDetect(&figureCoefficientRedCircle, &figureCoefficientWhiteInnerCircle, image, whiteBoardWithCircle, 4);
        //cutNoise(whiteBoardWithCircle, 3);

        cv::Mat whiteBoardForSpeedValues = makeWhiteBoard(image);
        //findPixelsOfSpeedValues(whiteBoardWithCircle, image, whiteBoardForSpeedValues);
        FigureCoefficient figureCoefficientMain;
        int colorsBlack[6] = {0, 0, 0, 100, 100, 100};
        figureCoefficientMain.setColors(colorsBlack);

        std::vector<int> vecXs;
        std::vector<int> vecYs;
        figureCoefficientMain.setCoorMaxY(figureCoefficientRedCircle.getCoorMaxY());
        figureCoefficientMain.setCoorMinY(figureCoefficientRedCircle.getCoorMinY());
        figureCoefficientMain.setCoorMaxX(figureCoefficientRedCircle.getCoorMaxX());
        figureCoefficientMain.setCoorMinX(figureCoefficientRedCircle.getCoorMinX());

        findStartingOfNumbers(&figureCoefficientMain, whiteBoardWithCircle, image, whiteBoardForSpeedValues, vecXs,
                              vecYs);
        std::vector<std::vector<cv::Point>> points;
        std::vector<std::vector<long>> pointsHash;
        for (int j = 0; j < vecXs.size(); ++j) {
            int firstX = vecXs[j];
            int firstY = vecYs[j];
            std::vector<cv::Point> checkedPixels;
            std::vector<long> hash;
            findNeighborhood(whiteBoardWithCircle, image, whiteBoardForSpeedValues, firstX, firstY,
                             figureCoefficientMain, checkedPixels);
            for (auto p:checkedPixels) {
                hash.push_back(p.x * 512 + p.y);
            }

            std::sort(hash.begin(), hash.end());
            if (checkedPixels.size() > 8) {
                points.push_back(checkedPixels);
                pointsHash.push_back(hash);
            }
        }

        int k = 0;
        // intersecty
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

        int ii = 0;
        std::cout << "Groups num: " << points.size() << "\n";
        if (points.size() == 3) {
            std::cout << "1";
            computeMoments(&figureCoefficientMain, points[1]);
            if (figureCoefficientMain.getM7() < 0.019) {
                std::cout << "1";
            }
            if (figureCoefficientMain.getM2() >= 0.14 && figureCoefficientMain.getM2() < 0.28) {
                std::cout << "2";
            }
            if (figureCoefficientMain.getM10() < 0 && figureCoefficientMain.getM8() > -0.001) {
                std::cout << "4";
            }
            if (figureCoefficientMain.getM8() > 0.001 && figureCoefficientMain.getM1() > 0.4) {
                std::cout << "8";
            }
            if (figureCoefficientMain.getM3() < 0.003 && figureCoefficientMain.getM2() < 0.14) {
                //if (figureCoefficientMain.getM3() < 0.0067) {
                std::cout << "0";
            }
            std::cout << "0";
        } else if (points.size() == 2) {
            computeMoments(&figureCoefficientMain, points[0]);
            if (figureCoefficientMain.getM7() < 0.019) {
                std::cout << "1";
            }
            if (figureCoefficientMain.getM2() >= 0.14 && figureCoefficientMain.getM2() < 0.28) {
                std::cout << "2";
            }
            //if (figureCoefficientMain.getM10() < 0 && figureCoefficientMain.getM8() > -0.001) {
            if (figureCoefficientMain.getM10() < 0) {
                std::cout << "4";
            }
            if (figureCoefficientMain.getM8() < 0.001 && figureCoefficientMain.getM1() > 0.4) {
                std::cout << "8";
            }
            std::cout << "0";
        } else if (points.size() == 1) {
            computeMoments(&figureCoefficientMain, points[0]);
            if (figureCoefficientMain.getM7() < 0.019) {
                std::cout << "1";
            }
            if (figureCoefficientMain.getM2() >= 0.14 && figureCoefficientMain.getM2() < 0.28) {
                std::cout << "2";
            }
            if (figureCoefficientMain.getM10() < 0 && figureCoefficientMain.getM8() > -0.001) {
                std::cout << "4";
            }
            //if (figureCoefficientMain.getM8() > 0.001 && figureCoefficientMain.getM1() > 0.57) {
            if (figureCoefficientMain.getM8() < 0 && figureCoefficientMain.getM1() > 0.4) {
                std::cout << "8";
            }
        } else {
            std::cout << "Tu nie ma znaku ograniczenia prędkości. ";
        }


        /*for (std::vector<cv::Point> checkedPixels: points) {


            int minx=image.rows;
            int miny=image.cols;
            for (auto p:checkedPixels) {
                if (p.x<minx) {
                    minx=p.x;
                }
                if (p.y<miny) {
                    miny=p.y;
                }

            }
            for (auto &p:checkedPixels) {
                p.x=p.x-minx;
                p.y=p.y-miny;
            }
            computeMoments(&figureCoefficientMain, checkedPixels);

            if (checkedPixels.size() > 10) {
                if (figureCoefficientMain.getM7() < 0.019) {
                    std::cout << "1";
                }
                if (figureCoefficientMain.getM2() >= 0.14 && figureCoefficientMain.getM2() < 0.28) {
                    std::cout << "2";
                }
                if (figureCoefficientMain.getM10() < 0 && figureCoefficientMain.getM8() > -0.001) {
                    std::cout << "4";
                }
                //if (figureCoefficientMain.getM8() > 0.001 && figureCoefficientMain.getM1() > 0.57) {
                if (figureCoefficientMain.getM8() < 0 && figureCoefficientMain.getM1() > 0.4) {
                    std::cout << "8";
                }
                if (figureCoefficientMain.getM3() < 0.003 && figureCoefficientMain.getM2() < 0.14) {
                    //if (figureCoefficientMain.getM3() < 0.0067) {
                    std::cout << "0";
                }
            }


            CV_Assert(whiteBoardWithCircle.depth() != sizeof(uchar));
            switch (whiteBoardWithCircle.channels()) {
                case 3:
                    cv::Mat_<cv::Vec3b> _I = whiteBoardWithCircle;
                    for (auto p:checkedPixels) {
                        _I(p.x, p.y)[0] = 100 + ii * 20;
                        _I(p.x, p.y)[1] = 50 + ii * 30;
                        _I(p.x, p.y)[2] = 100 + ii * 10;
                    }
            }
            ii = ii + 50;
        }*/
        std::cout << "\n";


        //computeField(&figureCoefficientRedCircle, speedLimitSign1);

        //cv::Mat image2 = speedLimitSign1(cv::Rect(figureCoefficientRedCircle.getPixels()[0], figureCoefficientRedCircle.getPixels()[1],
        //                                          figureCoefficientRedCircle.getPixels()[2] - figureCoefficientRedCircle.getPixels()[0], figureCoefficientRedCircle.getPixels()[3]
        //                                                                                                                                 - figureCoefficientRedCircle.getPixels()[1]));
        //cv::imshow("Image" + std::to_string(i), image);
        cv::imshow("Shape" + std::to_string(i), whiteBoardWithCircle);
        cv::imshow("Values" + std::to_string(i), whiteBoardForSpeedValues);

        i++;
    }
}


int main(int, char *[]) {
    std::cout << "Start ..." << std::endl;

    findSpeedLimitSign();


    cv::waitKey(-1);



    //countCoefficientFoNumbers();




    /*cv::Mat speedLimitSign1 = cv::imread("SourceImages/road112_AgnieszkaJurkiewicz.png");
    cv::Mat speedLimitSign2 = cv::imread("SourceImages/road113_AgnieszkaJurkiewicz.png");
    cv::Mat speedLimitSign3 = cv::imread("SourceImages/road119_AgnieszkaJurkiewicz.png");
    cv::Mat imagesList [3] = {speedLimitSign1, speedLimitSign2, speedLimitSign3};

    int i = 1;
    for (cv::Mat image : imagesList) {
        cv::Mat whiteBoard = makeWhiteBoard(image);
        cv::Mat clearWhiteBoard = makeWhiteBoard(image);

        FigureCoefficient figureCoefficientRedCircle;
        int colorsRed[6] = {0, 0, 70, 93, 82, 255}; //czerwony
        //int colorsRed [6] = {110, 100, 150, 170, 160, 200}; // różowy
        //int colorsRed [6] = {0, 0, 40, 20, 20, 250};
        figureCoefficientRedCircle.setColors(colorsRed);

        FigureCoefficient figureCoefficientWhiteInnerCircle;
        //int colorsWhite [6] = {130, 120, 155, 255, 255, 255};
        int colorsWhite[6] = {155, 155, 155, 255, 255, 255};
        figureCoefficientWhiteInnerCircle.setColors(colorsWhite);

        computeBox(image, &figureCoefficientRedCircle);

        //findNeighborhood(&figureCoefficientRedCircle, &figureCoefficientWhiteInnerCircle, speedLimitSign1, whiteBoard, 4, 0, 0, 0);
        edgeDetect(&figureCoefficientRedCircle, &figureCoefficientWhiteInnerCircle, image, whiteBoard, 4);
        //cutNoise(whiteBoard, 3);

        cv::Mat whiteBoardForSpeedValues = makeWhiteBoard(image);
        //findPixelsOfSpeedValues(whiteBoard, image, whiteBoardForSpeedValues);
        findStartingOfNumbers(whiteBoard, image, whiteBoardForSpeedValues);

        //computeField(&figureCoefficientRedCircle, speedLimitSign1);

        //cv::Mat image2 = speedLimitSign1(cv::Rect(figureCoefficientRedCircle.getPixels()[0], figureCoefficientRedCircle.getPixels()[1],
        //                                          figureCoefficientRedCircle.getPixels()[2] - figureCoefficientRedCircle.getPixels()[0], figureCoefficientRedCircle.getPixels()[3]
        //                                                                                                                                 - figureCoefficientRedCircle.getPixels()[1]));
        //cv::imshow("Image" + std::to_string(i), image);
        //cv::imshow("Shape" + std::to_string(i), whiteBoard);
        cv::imshow("Values" + std::to_string(i), whiteBoardForSpeedValues);
        i++;
    }*/




    return 0;
}
