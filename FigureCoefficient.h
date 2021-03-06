//
// Created by Agnieszka Jurkiewicz on 11/05/2021.
//

#ifndef POBRLAB1_FIGURECOEFFICIENT_H
#define POBRLAB1_FIGURECOEFFICIENT_H

#include <string>


class FigureCoefficient {
public:
    FigureCoefficient();

    virtual ~FigureCoefficient();

    const std::string &getFileName() const;

    void setFileName(const std::string &fileName);

    const int *getColors() const;

    void addColor(int color, int i);

    void setColors(int colors[]);

    int getField() const;

    void setField(int field);

    int getCircumference() const;

    void setCircumference(int circumference);

    float getW3() const;

    void setW3(float w3);

    float getM1() const;

    void setM1(float m1);

    float getM7() const;

    void setM7(float m7);

    float getAngle() const;

    void setAngle(float angle);

    int getCoorMinX() const;

    void setCoorMinX(int coorMinX);

    int getCoorMaxX() const;

    void setCoorMaxX(int coorMaxX);

    int getCoorMinY() const;

    void setCoorMinY(int coorMinY);

    int getCoorMaxY() const;

    void setCoorMaxY(int coorMaxY);

    float getM2() const;

    void setM2(float m2);

    float getM3() const;

    void setM3(float m3);

    float getM4() const;

    void setM4(float m4);

    float getM8() const;

    void setM8(float m8);

    float getM10() const;

    void setM10(float m10);

private:
    std::string fileName;
    int coorMinX;
    int coorMaxX;
    int coorMinY;
    int coorMaxY;
    int colors[6];
    int field;
    int circumference;
    float W3;
    float M1;
    float M2;
    float M3;
    float M4;
    float M7;
    float M8;
    float M10;
    float angle;

};


#endif //POBRLAB1_FIGURECOEFFICIENT_H
