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
                if (y>=I.cols) {
                    y=I.cols-1;
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
                    x = I.rows-1;
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
                if (x==I.rows) {
                    x=I.rows-1;
                }
            }

            x = I.rows - 1;
            y = I.cols - 1;
            int coorDownX = x;
            int coorDownY = y;
            while (!computeConditionForColor(_I, figureCoefficient, x, y) && x >= 0) {
                --x;
                y = I.cols - 1;
                while (!computeConditionForColor(_I, figureCoefficient, x, y) && y >= 0 ) {
                    coorDownX = x;
                    coorDownY = y;
                    y--;
                }
                if (y<0) {
                    y=0;
                }
            }

            figureCoefficient->setCoorMinX(coorUpX);
            figureCoefficient->setCoorMinY(coorLeftY);
            figureCoefficient->setCoorMaxX(coorDownX);
            figureCoefficient->setCoorMaxY(coorRightY);
            std::cout << "coorDownX = " << figureCoefficient->getCoorMinX() << ", coorUpX = " << figureCoefficient->getCoorMaxX()
                    << ", coorLeftY = " << figureCoefficient->getCoorMinY() << ", coorRightY = " << figureCoefficient->getCoorMaxY()<< std::endl;
            std::cout << "Width = " << figureCoefficient->getCoorMaxY() - figureCoefficient->getCoorMinY() << "Hight = " << figureCoefficient->getCoorMaxX()
                                                                                 - figureCoefficient->getCoorMinX() << std::endl;
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
            for (int j = figureCoefficient->getCoorMinX(); j <= figureCoefficient->getCoorMaxX(); ++j) {
                for (int i = figureCoefficient->getCoorMinY(); i <= figureCoefficient->getCoorMaxY(); ++i) {
                    if (computeConditionForColor(_I, figureCoefficient, i, j)) {
                        field = field + 1;
                        _I(i, j)[0] = 255;
                        _I(i, j)[1] = 255;
                        _I(i, j)[2] = 255;
                    }
                }
            }
    }
    figureCoefficient->setField(field);
    std::cout <<  ": S = " << figureCoefficient->getField();
}

void findNeighborhood(FigureCoefficient *figureCoefficientMain, cv::Mat whiteBoardForSpeedValues, int firstX, int firstY) {
    std::vector<cv::Point> pixelsToCheck;
    cv::Point pixelFirst(firstX, firstY);
    pixelsToCheck.push_back(pixelFirst);
    std::vector<cv::Point> checkedPixels;

    CV_Assert(whiteBoardForSpeedValues.depth() != sizeof(uchar));
    cv::Mat res(whiteBoardForSpeedValues.rows, whiteBoardForSpeedValues.cols, CV_8UC3);
    switch (whiteBoardForSpeedValues.channels()) {
        case 3:
            cv::Mat_<cv::Vec3b> _R = whiteBoardForSpeedValues;
            while (!pixelsToCheck.empty()) {
                cv::Point pixel = pixelsToCheck.back();
                pixelsToCheck.pop_back();
                checkedPixels.push_back(pixel);
                for (int x = pixel.x - 1; x <= pixel.x + 1; x++) {
                    for (int y = pixel.y - 1; y <= pixel.y + 1; y++) {
                        if (computeConditionForColor(_R, figureCoefficientMain, x, y)) {
                            cv::Point pixelNext(x,y);
                            if (!std::count(checkedPixels.begin(), checkedPixels.end(), pixelNext)) {
                                pixelsToCheck.push_back(pixelNext);
                            }
                        }
                    }
                }
            }
    }
}

void findNeighborhoodNotDiagonally (FigureCoefficient *figureCoefficientWhiteBoard, cv::Mat whiteBoard, int sizeOfNeighborhood, int r, int g, int b) {
    CV_Assert(whiteBoard.depth() != sizeof(uchar));
    cv::Mat res(whiteBoard.rows, whiteBoard.cols, CV_8UC3);
    switch (whiteBoard.channels()) {
        case 3:
            cv::Mat_<cv::Vec3b> _R = whiteBoard;
            for (int x = figureCoefficientWhiteBoard->getCoorMinX() + sizeOfNeighborhood; x <= figureCoefficientWhiteBoard->getCoorMaxX() - sizeOfNeighborhood; x++) {
                for (int y = figureCoefficientWhiteBoard->getCoorMinY() + sizeOfNeighborhood; y <= figureCoefficientWhiteBoard->getCoorMaxY() - sizeOfNeighborhood; y++) {
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

void erosion(cv::Mat whiteBoard, FigureCoefficient *figureCoefficientWhiteBoard) {
    int sizeOfNeighborhood = 1;
    int r = 100;
    int g = 1;
    int b = 100;

    findNeighborhoodNotDiagonally(figureCoefficientWhiteBoard, whiteBoard, sizeOfNeighborhood, r, g, b);
    FigureCoefficient figureCoefficientNotBlackPoints;
    int colorsNotBlack [6] = {1, 1, 1, 200, 200, 200};
    figureCoefficientNotBlackPoints.setColors(colorsNotBlack);

    CV_Assert(whiteBoard.depth() != sizeof(uchar));
    cv::Mat res(whiteBoard.rows, whiteBoard.cols, CV_8UC3);
    switch (whiteBoard.channels()) {
        case 3:
            cv::Mat_<cv::Vec3b> _R = whiteBoard;
            for (int x = 0; x <= whiteBoard.cols; ++x) {
                for (int y = 0; y <= whiteBoard.rows; ++y) {
                    if (computeConditionForColor(_R, &figureCoefficientNotBlackPoints, x, y)) {
                            _R(x, y)[0] = 255;
                            _R(x, y)[1] = 255;
                            _R(x, y)[2] = 255;
                    }
                }
            }
    }

}

void edgeDetect(FigureCoefficient *figureCoefficientMain, FigureCoefficient *figureCoefficientNeighbors, cv::Mat &I, cv::Mat whiteBoard, int sizeOfNeighborhood) {
    FigureCoefficient figureCoefficientBlackNumbers;
    int colorsBlack[6] = {0, 0, 0, 20, 20, 20};
    figureCoefficientBlackNumbers.setColors(colorsBlack);
    CV_Assert(I.depth() != sizeof(uchar));
    cv::Mat res(I.rows, I.cols, CV_8UC3);
    switch (I.channels()) {
        case 3:
            cv::Mat_<cv::Vec3b> _I = I;
            cv::Mat_<cv::Vec3b> _R = whiteBoard;
            for (int x = figureCoefficientMain->getCoorMinX() + sizeOfNeighborhood; x <= figureCoefficientMain->getCoorMaxX() - sizeOfNeighborhood; ++x) {
                for (int y = figureCoefficientMain->getCoorMinY() + sizeOfNeighborhood; y <= figureCoefficientMain->getCoorMaxY() - sizeOfNeighborhood; ++y) {
                    int cwhite=0;
                    int cred=0;
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
                        if (cred>25 and (cwhite>20 or cwhite+cblack>20)) {
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

    figureCoefficientMain->setCoorMinY(figureCoefficientFrameForSpeedValue.getCoorMinY() + (int) width/10);
    figureCoefficientMain->setCoorMaxY(figureCoefficientFrameForSpeedValue.getCoorMaxY() - (int) width/10);
    figureCoefficientMain->setCoorMinX(figureCoefficientFrameForSpeedValue.getCoorMinX() + (int) hight/6);
    figureCoefficientMain->setCoorMaxX(figureCoefficientFrameForSpeedValue.getCoorMaxX() - (int) hight/6);
}

void findPixelsOfSpeedValues(cv::Mat whiteBoard, cv::Mat image, cv::Mat whiteBoardForSpeedValues, FigureCoefficient figureCoefficientMain) {

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

void findFirstBlackPixel(cv::Mat whiteBoardForSpeedValues, int firstX, int firstY, FigureCoefficient figureCoefficientMain) {
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
                        firstX = x;
                        firstY = y;
                        isBlack = true;
                    }
                    y--;
                }
                x--;
            }
    }
}

void findStartingPointsForAllNumbers(cv::Mat whiteBoardForSpeedValues, int firstX, int firstY, FigureCoefficient figureCoefficientMain,
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

void findStartingOfNumbers(cv::Mat whiteBoard, cv::Mat image, cv::Mat whiteBoardForSpeedValues) {

    FigureCoefficient figureCoefficientMain;
    int colorsBlack [6] = {0, 0, 0, 80, 80, 80};
    figureCoefficientMain.setColors(colorsBlack);
    findPixelsOfSpeedValues(whiteBoard, image, whiteBoardForSpeedValues, figureCoefficientMain);

    int firstX;
    int firstY;
    findFirstBlackPixel(whiteBoardForSpeedValues, firstX, firstY, figureCoefficientMain);

    std::vector<std::vector<int>> vecX;
    std::vector<std::vector<int>> vecY;

    //findStartingPointsForAllNumbers(whiteBoardForSpeedValues, firstX, firstY, figureCoefficientMain, vecX, vecY);

    std::vector<int> vecXs;
    std::vector<int> vecYs;
    bool isBlack = true;

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
    vecXs = vecX[0];
    vecYs = vecY[0];
    for (int k =1; k < 4; k++) {
        if (vecX[k].size() > vecXs.size()) {
            vecXs.clear();
            vecXs = vecX[k];
            vecY.clear();
            vecYs = vecY[k];
        }
    }

}


void cutNoise(cv::Mat whiteBoard, int times) {

    FigureCoefficient figureCoefficientWhiteBoard;
    int colorsBlack [6] = {0, 0, 0, 200, 200, 200};
    figureCoefficientWhiteBoard.setColors(colorsBlack);
    computeBox(whiteBoard, &figureCoefficientWhiteBoard);

    for (int i = 0; i <= times; i++) {
        erosion(whiteBoard, &figureCoefficientWhiteBoard);
    }

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
            for (int j = figureCoefficient->getCoorMinX(); j <= figureCoefficient->getCoorMaxX(); ++j) {
                for (int i = figureCoefficient->getCoorMinY(); i <= figureCoefficient->getCoorMaxY(); ++i) {
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
    float middleX = (figureCoefficient->getCoorMaxX() - figureCoefficient->getCoorMinX())/ 2 + figureCoefficient->getCoorMinX();
    float middleY = (figureCoefficient->getCoorMaxY() - figureCoefficient->getCoorMinY())/ 2 + figureCoefficient->getCoorMinY();
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
            for (int j = figureCoefficient->getCoorMinX(); j <= figureCoefficient->getCoorMaxX(); ++j) {
                for (int i = figureCoefficient->getCoorMinY(); i <= figureCoefficient->getCoorMaxY(); ++i) {
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


    std::string fileNames[4] = {"0a", "0b", "0c", "0d"};
    for (int i = 0; i < 4; i++) {
        FigureCoefficient figureCoefficient;
        int colorsBlack[6] = {0, 0, 0, 100, 100, 100};
        figureCoefficient.setColors(colorsBlack);
        figureCoefficient.setFileName(fileNames[i]);
        cv::Mat image = cv::imread("SourceImages/" + fileNames[i] + ".png");

        computeBox(image, &figureCoefficient);

        //cv::Mat image2 = image(cv::Rect(figureCoefficient.getCoorMinY(), figureCoefficient.getCoorMinX(),
        //                                figureCoefficient.getCoorMaxY() - figureCoefficient.getCoorMinY(),
        //                                figureCoefficient.getCoorMaxX()
        //                                - figureCoefficient.getCoorMinX()));
        std::cout << i + 1 << ". Plik " + fileNames[i];
        cv::imshow("Shape", image);
        computeField(&figureCoefficient, image);

        computeCircumference(&figureCoefficient, image);
        computeCoefficientOfMalinowska(&figureCoefficient);

        computeMoments(&figureCoefficient, image);
        std::cout << ", M1 = " << figureCoefficient.getM1() << ", M7 = " << figureCoefficient.getM7() << std::endl;


        cv::waitKey(-1);
    }




    return 0;
}
