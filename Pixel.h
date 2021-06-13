//
// Created by Agnieszka Jurkiewicz on 16/04/2021.
//

#ifndef POBRLAB1_PIXEL_H
#define POBRLAB1_PIXEL_H
#include <iostream>


class Pixel {

public:
    Pixel();

    Pixel(int coordinateX, int coordinateY, float lightness);

    Pixel(int coordinateX, int coordinateY);

    //Pixel (const Pixel &pixel);

    virtual ~Pixel();

private:

    int coordinateX;
    int coordinateY;
    float lightness;

public:
    int getCoordinateX() const;

    void setCoordinateX(int coordinateX);

    int getCoordinateY() const;

    void setCoordinateY(int coordinateY);

    float getLightness() const;

    void setLightness(float lightness);

    //friend bool operator< (const Pixel & a, const Pixel & b);

};


#endif //POBRLAB1_PIXEL_H
