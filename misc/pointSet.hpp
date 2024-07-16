#pragma once

#include <string>
#include <random>
#include <fstream>
#include "shape.hpp"


bool loadPointSet(std::string filename, PointSet* _pointSet);
void addNoiseToPointSet(PointSet* _pointSet, double amount, char mode);
void generatePointSet(Shape* _data, PointSet* _pointSet, bool deformed = false);
