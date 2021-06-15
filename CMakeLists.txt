cmake_minimum_required(VERSION 3.15)
project(TrafficSignRecognizer)

set(CMAKE_CXX_STANDARD 11)

#find_package(OpenCV REQUIRED core imgproc highgui)
find_package(OpenCV REQUIRED)
message(WARNING "After setting - ${OpenCV_INCLUDE_DIRS}")
include_directories(${OpenCV_INCLUDE_DIRS})

set(OpenCV_LIBS opencv_core opencv_imgproc opencv_highgui)


add_executable(TrafficSignRecognizer TrafficSignRecognizer.cpp FigureCoefficient.cpp FigureCoefficient.h)
target_link_libraries(TrafficSignRecognizer ${OpenCV_LIBS})