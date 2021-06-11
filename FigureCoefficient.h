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

    const int *getPixels() const;

    void addPixel(int pixel, int i);

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

private:
    std::string fileName;
    int pixels[4];
    int colors[6];
    int field;
    int circumference;
    float W3;
    float M1;
    float M7;
    float angle;

};


#endif //POBRLAB1_FIGURECOEFFICIENT_H
