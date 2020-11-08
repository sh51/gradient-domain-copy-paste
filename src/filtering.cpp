// filtering.cpp
// Assignment 4

#include "filtering.h"
#include "a2.h"
#include <math.h>

using namespace std;

/**************************************************************
 //            IMAGE CONVOLUTION AND FILTERING               //
 *************************************************************/


// convolve an image with a box filter of size k by k
FloatImage boxBlur(const FloatImage &im, const int &k, bool clamp)
{
	// create a new empty image
	FloatImage imFilter(im.width(), im.height(), im.channels());


	// convolve the image with the box filter


	return imFilter;
}

// reimeplement the box filter using the filter class.
// check that your results math those in the previous function "boxBlur"
FloatImage boxBlur_filterClass(const FloatImage &im, const int &k, bool clamp)
{
	// use Convolve() to apply convolution

	return FloatImage();// CHANGEME
}

// uses a Sobel kernel to compute the horizontal and vertical
// components of the gradient of an image and returns the gradient magnitude.
FloatImage gradientMagnitude(const FloatImage &im, bool clamp)
{
	// sobel filtering in x direction


	// sobel filtering in y direction


	// compute squared magnitude



	return FloatImage();//  CHANGEME
}

// create a vector containing the normalized values in a 1D Gaussian filter
vector<float> gauss1DFilterValues(float sigma, float truncate)
{
	// calculate the size of the filter

	// compute the un-normalized value of the gaussian

	// normalize


	return vector <float>();//  CHANGEME
}

// blur across the rows of an image
FloatImage gaussianBlur_horizontal(const FloatImage &im, float sigma, float truncate, bool clamp)
{
	// filter in the x direction

	return FloatImage();// CHANGEME
}

// create a vector containing the normalized values in a 2D Gaussian filter
vector<float> gauss2DFilterValues(float sigma, float truncate)
{
	// compute the filter size

	// compute the unnormalized value of the gaussian and put it in a row-major vector

	// normalize

	return vector<float>();// CHANGEME
}

// Blur an image with a full  full 2D rotationally symmetric Gaussian kernel
FloatImage gaussianBlur_2D(const FloatImage &im, float sigma, float truncate, bool clamp)
{
	// blur using a 2D gaussian filter (use gauss2DFilterValues())

	return FloatImage();//  CHANGEME
}

// Use principles of seperabiltity to blur an image using 2 1D Gaussian Filters
FloatImage gaussianBlur_seperable(const FloatImage &im, float sigma, float truncate, bool clamp)
{
	// blur using 2, 1D filters in the x and y directions

	return FloatImage();//  CHANGEME
}

// sharpen an image
FloatImage unsharpMask(const FloatImage &im, float sigma, float truncate, float strength, bool clamp)
{
	// get the low pass image and subtract it from the original image to get the high pass image

	return FloatImage();//  CHANGEME
}

// Denoise an image using bilateral filtering
FloatImage bilateral(const FloatImage &im, float sigmaRange, float sigmaDomain, float truncateDomain, bool clamp)
{
	// calculate the filter size

	// for every pixel (x,y) in the image set value to weighted sum of values in the filter region

	return FloatImage();//  CHANGEME
}

// Bilaterial Filter an image seperately for the Y and UV components of an image
FloatImage bilaYUV(const FloatImage &im, float sigmaRange, float sigmaY, float sigmaUV, float truncateDomain,
                   bool clamp)
{
	// convert from RGB to YUV

	// denoise Y and UV channels using different domain sigmas

	// convert from YUV back to RGB


	return FloatImage();//  CHANGEME
}


/**************************************************************
 //                 FILTER CLASS FUNCTIONS                  //
 *************************************************************/


// write a convolution function for the filter class
FloatImage Filter::Convolve(const FloatImage &im, bool clamp) const
{
	FloatImage imFilter(im.width(), im.height(), im.channels());

	// implement convultion
	// Hint: use use Filter::operator()(x, y) to access (x,y) kernel location

	return imFilter;
}


/**************************************************************
 //               DON'T EDIT BELOW THIS LINE                //
 *************************************************************/

// Create an image of 0's with a value of 1 in the middle. This function
// can be used to test that you have properly set the kernel values in your
// Filter object. Make sure to set k to be larger than the size of your kernel
FloatImage impulseImg(const int &k)
{
	// initlize a kxkx1 image of all 0's
	FloatImage impulse(k, k, 1);

	// set the center pixel to have intensity 1
	int center = floor(k / 2);
	impulse(center, center, 0) = 1;

	return impulse;
}

Filter::Filter(const vector<float> &fData, const int &fWidth, const int &fHeight)
{
	// TODO: check that width*height = length of filterVals and that width and height are odd

	kernel = fData;
	width = fWidth;
	height = fHeight;

}

Filter::Filter(const int &fWidth, const int &fHeight)
{
	width = fWidth;
	height = fHeight;
	kernel = std::vector<float>(width * height, 0);
}

const float &Filter::operator()(int x, int y) const
{
	if (x < 0 || x >= width)
		throw OutOfBoundsException();
	if (y < 0 || y >= height)
		throw OutOfBoundsException();

	return kernel[x + y * width];
}

float &Filter::operator()(int x, int y)
{
	if (x < 0 || x >= width)
		throw OutOfBoundsException();
	if (y < 0 || y >= height)
		throw OutOfBoundsException();

	return kernel[x + y * width];
}
Filter::~Filter() {}
