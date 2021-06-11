//
// Created by Agnieszka Jurkiewicz on 11/05/2021.
//

#include "FigureCoefficient.h"

FigureCoefficient::FigureCoefficient() {
    fileName = "elipsa";
    field = 0;
    circumference = 0;
    W3 = 0;
    M1 = 0;
    M7 = 0;
}
FigureCoefficient::~FigureCoefficient(){}

int FigureCoefficient::getField() const {
    return field;
}

void FigureCoefficient::setField(int field) {
    FigureCoefficient::field = field;
}

int FigureCoefficient::getCircumference() const {
    return circumference;
}

void FigureCoefficient::setCircumference(int circumference) {
    FigureCoefficient::circumference = circumference;
}

float FigureCoefficient::getW3() const {
    return W3;
}

void FigureCoefficient::setW3(float w3) {
    W3 = w3;
}

float FigureCoefficient::getM1() const {
    return M1;
}

void FigureCoefficient::setM1(float m1) {
    M1 = m1;
}

float FigureCoefficient::getM7() const {
    return M7;
}

void FigureCoefficient::setM7(float m7) {
    M7 = m7;
}

const std::string &FigureCoefficient::getFileName() const {
    return fileName;
}

void FigureCoefficient::setFileName(const std::string &fileName) {
    FigureCoefficient::fileName = fileName;
}

const int *FigureCoefficient::getPixels() const {
    return pixels;
}

void FigureCoefficient::addPixel(int pixel, int i) {
    pixels[i] = pixel;
}

const int *FigureCoefficient::getColors() const {
    return colors;
}

void FigureCoefficient::addColor(int color, int i) {
    colors[i] = color;
}

float FigureCoefficient::getAngle() const {
    return angle;
}

void FigureCoefficient::setAngle(float angle) {
    FigureCoefficient::angle = angle;
}
