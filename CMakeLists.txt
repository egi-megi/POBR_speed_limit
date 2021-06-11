cmake_minimum_required(VERSION 3.15)
project(POBRlab1)

set(CMAKE_CXX_STANDARD 11)

#find_package(OpenCV REQUIRED core imgproc highgui)
find_package(OpenCV REQUIRED)
message(WARNING "After setting - ${OpenCV_INCLUDE_DIRS}")
include_directories(${OpenCV_INCLUDE_DIRS})

set(OpenCV_LIBS opencv_core opencv_imgproc opencv_highgui)


add_executable(POBRlab1 Framework3.cpp Pixel.cpp Pixel.h FigureCoefficient.cpp FigureCoefficient.h)
target_link_libraries(POBRlab1 ${OpenCV_LIBS})