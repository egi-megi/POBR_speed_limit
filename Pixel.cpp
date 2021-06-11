//
// Created by Agnieszka Jurkiewicz on 16/04/2021.
//

#include "Pixel.h"

Pixel::Pixel() {
    coordinateX = 0;
    coordinateY = 0;
    lightness = 0;
}

Pixel::Pixel(int coordinateX, int coordinateY, float lightness) : coordinateX(coordinateX), coordinateY(coordinateY),
                                                                  lightness(lightness) {}

Pixel::~Pixel() {

}

int Pixel::getCoordinateX() const {
    return coordinateX;
}

void Pixel::setCoordinateX(int coordinateX) {
    Pixel::coordinateX = coordinateX;
}

int Pixel::getCoordinateY() const {
    return coordinateY;
}

void Pixel::setCoordinateY(int coordinateY) {
    Pixel::coordinateY = coordinateY;
}

float Pixel::getLightness() const {
    return lightness;
}

void Pixel::setLightness(float lightness) {
    Pixel::lightness = lightness;
}

/*
bool operator<(const Pixel & a, const Pixel & b)
{
    return a.getLightness() <= b.getLightness();
}*/






