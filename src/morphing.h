#pragma once

// morphing.h
// Assignment 6

#include <array>
#include "floatimage.h"
#include "exceptions.h"
#include <iostream>
#include <math.h>

using Vec2f = std::array<float, 2>;

class Segment
{
public:
    Vec2f P;                    //!< P - 2-element vector to point (x1, y1)
    Vec2f Q;                    //!< Q - 2-element vector to pont (x2, y2)
    Vec2f PQ;                   //!< PQ - 2-element vector from P to Q
    Vec2f perpPQ;               //!< perpPQ - 2-element vector perpendicular to PQ
    Vec2f PQDivByPQlength2;     //!< PQDivByPQlength2 - 2-element vector PQ normalized by PQ2
    Vec2f perpPQDivByPQlength;  //!< perpPQDivByPQlength - 2-element vector perpPQ normalized by PQlength
    float PQ2;                  //!< PQ2 - squared distance between P and Q
    float PQlength;             //!< PQlength - distance between P and Q
    
    float getU(float x, float y);
    float getV(float x, float y);
    float dist(float x, float y);
    Vec2f UVtoX(float u, float v);
    float weight(float x, float y, float a, float b, float p);
    
    //Constructor
    Segment(float x1, float y1, float x2, float y2);
    
    // Destructor. Because there is no explicit memory management here, this doesn't do anything
    ~Segment();
    
    // The following are static helper functions
    static float dot(const Vec2f &vec1, const Vec2f &vec2);
    static Vec2f add(const Vec2f &vec1, const Vec2f &vec2);
	static Vec2f subtract(const Vec2f &vec1, const Vec2f &vec2);
    static Vec2f scalarMult(const Vec2f &vec, float factor);

};

FloatImage warpBy1(const FloatImage &im, Segment &segBefore, Segment &segAfter); 
FloatImage warp(const FloatImage &im, std::vector<Segment> &segsBefore, std::vector<Segment> &segsAfter, float a=10.0, float b=1.0, float p=1.0);
std::vector<FloatImage> morph(const FloatImage &im1, const FloatImage &im2, std::vector<Segment> &segsBefore, std::vector<Segment> &segsAfter, int N=1, float a=10.0, float b=1.0, float p=1.0);
