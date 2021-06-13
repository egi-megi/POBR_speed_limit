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

float FigureCoefficient::getM2() const {
    return M2;
}

void FigureCoefficient::setM2(float m2) {
    M2 = m2;
}

float FigureCoefficient::getM3() const {
    return M3;
}

void FigureCoefficient::setM3(float m3) {
    M3 = m3;
}

float FigureCoefficient::getM4() const {
    return M4;
}

void FigureCoefficient::setM4(float m4) {
    M4 = m4;
}

float FigureCoefficient::getM7() const {
    return M7;
}

void FigureCoefficient::setM7(float m7) {
    M7 = m7;
}

float FigureCoefficient::getM8() const {
    return M8;
}

void FigureCoefficient::setM8(float m8) {
    M8 = m8;
}

float FigureCoefficient::getM10() const {
    return M10;
}

void FigureCoefficient::setM10(float m10) {
    M10 = m10;
}

const std::string &FigureCoefficient::getFileName() const {
    return fileName;
}

void FigureCoefficient::setFileName(const std::string &fileName) {
    FigureCoefficient::fileName = fileName;
}

const int *FigureCoefficient::getColors() const {
    return colors;
}

void FigureCoefficient::addColor(int color, int i) {
    colors[i] = color;
}

void FigureCoefficient::setColors(int colorsA[]) {
    for (int i = 0; i < 6; i++)
        colors[i] = colorsA[i];
}

float FigureCoefficient::getAngle() const {
    return angle;
}

void FigureCoefficient::setAngle(float angle) {
    FigureCoefficient::angle = angle;
}

int FigureCoefficient::getCoorMinX() const {
    return coorMinX;
}

void FigureCoefficient::setCoorMinX(int coorMinX) {
    FigureCoefficient::coorMinX = coorMinX;
}

int FigureCoefficient::getCoorMaxX() const {
    return coorMaxX;
}

void FigureCoefficient::setCoorMaxX(int coorMaxX) {
    FigureCoefficient::coorMaxX = coorMaxX;
}

int FigureCoefficient::getCoorMinY() const {
    return coorMinY;
}

void FigureCoefficient::setCoorMinY(int coorMinY) {
    FigureCoefficient::coorMinY = coorMinY;
}

int FigureCoefficient::getCoorMaxY() const {
    return coorMaxY;
}

void FigureCoefficient::setCoorMaxY(int coorMaxY) {
    FigureCoefficient::coorMaxY = coorMaxY;
}
