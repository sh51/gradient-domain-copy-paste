#include "a6.h"
#include <iostream>
#include <math.h>

using namespace std;

// create a new image that is k times bigger than the input by using nearest neighbor interpolation
FloatImage scaleNN(const FloatImage &im, float factor)
{
	// create new FloatImage that is factor times larger than input image
	FloatImage output(floor(im.sizeX()*factor),floor(im.sizeY()*factor),im.depth());

	// loop over new FloatImage pixels and set appropriate values (smartAccessor is probably overkill here)
	for(int i=0;i<output.sizeX();i++){
		for(int j=0;j<output.sizeY();j++){
			for(int k=0;k<im.depth();k++){	
				output(i,j,k)=im(round(i/factor),round(j/factor),k);
			}
		}
	}

	// return new image
	return output; // CHANGE ME
}

// using bilinear interpolation to assign the value of a location from its neighboring pixel values
float interpolateLin(const FloatImage &im, float x, float y, int z, bool clamp)
{
	// Hint: use smartAccessor() to handle coordinates outside the image
	//return final float value
	float f_x, f_y, row_1, row_2;
	f_x=x-floor(x);//dist to top
	f_y=y-floor(y);//dist to left
	//first row
	row_1=(1.f-f_y)*im.smartAccessor(floor(x),floor(y),z,clamp)+f_y*im.smartAccessor(floor(x),floor(y)+1,z,clamp);
	//second row
	row_2=(1.f-f_y)*im.smartAccessor(floor(x)+1,floor(y),z,clamp)+f_y*im.smartAccessor(floor(x)+1,floor(y)+1,z,clamp);
	//colume lerp
    return (1.f-f_x)*row_1+f_x*row_2; 
}

// create a new image that is k times bigger than the input by using bilinear interpolation
FloatImage scaleLin(const FloatImage &im, float factor)
{

	FloatImage output(floor(im.sizeX()*factor),floor(im.sizeY()*factor),im.depth());

    // loop over new FloatImage pixels and set appropriate values (use interpolateLin())
	for(int i=0;i<output.sizeX();i++){
		for(int j=0;j<output.sizeY();j++){
			for(int k=0;k<im.depth();k++){			
				output(i,j,k)=interpolateLin(im,i/factor,j/factor,k);

			}
			
		}
	}

	// return new image
    return output; // CHANGE ME
}

// rotate an image around its center by theta
FloatImage rotate(const FloatImage &im, float theta)
{
	FloatImage out(im.width(),im.height(),im.sizeZ());
	float x;
	float y;
	float x_c,y_c;
	
	//rotate im around its center by theta
	for(int i=0;i< out.sizeX();i++){
		for(int j=0;j<out.sizeY();j++){
			for(int k=0;k<im.depth();k++){
				//translate the rotation:x y-c
				x_c=i-im.width()/2;
				y_c=j-im.height()/2;
				//apply rotation
				x=x_c*cos(theta)-y_c*sin(theta);
				y=x_c*sin(theta)+y_c*cos(theta);
				//translate back: xy +c
				x_c=x+im.width()/2;
				y_c=y+im.height()/2;
				//transfomed position, bilinear
				out(i,j,k)= interpolateLin(im,x_c,y_c,k,false);
			}
		}
	}

	// return rotated image
    return out; 
}