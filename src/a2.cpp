#include "a2.h"
#include <math.h>
#include "utils.h"
#include <assert.h>
#include <iostream>

using namespace std;

/* Replace this file with the one that includes your implementations in assignment 2 */

/* 2 functions you've already implemented in assignment 1 */

// Change the brightness of the image
// const FloatImage & means a reference to im will get passed to you, 
// but the compiler won't let you modify it
FloatImage brightness(const FloatImage &im, float factor)
{
	// Create output FloatImage
	// Modify image brightness
	// return output;
	return FloatImage(1,1,1); // Change this
}

FloatImage contrast(const FloatImage &im, float factor, float midpoint)
{
	// Create output FloatImage
	// Modify image contrast
	// return output;
	return FloatImage(1,1,1); // Change this
}


/*New tasks in assignment 2 */

FloatImage changeGamma(const FloatImage &im, float old_gamma, float new_gamma)
{
	// create an appropriately sized output FloatImage
	FloatImage output(im);

	// Figure out what power to take the values of im, to get the values of output

	return output;
}

// Change the exposure of the image. This is different than brightness because
// it means you multiply an image's intensity in the linear domain.
FloatImage exposure(const FloatImage &im, float factor)
{
	// create an appropriately sized output FloatImage
	FloatImage output(im);
	return output;
}

FloatImage color2gray(const FloatImage &im, const vector<float> &weights)
{
	if ((int) weights.size() != im.channels())
		throw MismatchedDimensionsException();

	FloatImage output(im);  // CHANGEME

	// Convert to grayscale

	return output;
}

// For this function, we want two outputs, a single channel luminance image
// and a three channel chrominance image. Return them in a vector with luminance first
vector<FloatImage> lumiChromi(const FloatImage &im)
{
	vector<FloatImage> output;

	// Create luminance and chrominance images

	// push the luminance and chrominance images onto the vector

	return output;
}

// go from a luminance/chrominance images back to an rgb image
FloatImage lumiChromi2rgb(const FloatImage &lumi, const FloatImage &chromi)
{
	FloatImage output;

	// combine the luminance and chrominance images into an rgb image

	return output;
}

// Modify brightness then contrast
FloatImage brightnessContrastLumi(const FloatImage &im, float brightF, float contrastF, float midpoint)
{
	// Create output image of appropriate size
	FloatImage output(im);

	// Modify brightness, then contrast of luminance image

	return output;
}

FloatImage rgb2yuv(const FloatImage &im)
{
	// Create output image of appropriate size
	FloatImage output(im);

	// Change colorspace

	return output;
}

FloatImage yuv2rgb(const FloatImage &im)
{
	// Create output image of appropriate size
	FloatImage output(im);

	// Change colorspace

	return output;
}

FloatImage saturate(const FloatImage &im, float factor)
{
	// Create output image of appropriate size
	FloatImage output(im);

	// Saturate image

	return output;
}

// Return two images in a C++ vector
vector<FloatImage> spanish(const FloatImage &im)
{
	// Remember to create the output images and the output vector
	vector<FloatImage> output;

	// Do all the required processing

	// Push the images onto the vector

	// Return the vector
	return output;
}

// White balances an image using the gray world assumption
FloatImage grayworld(const FloatImage &im)
{
	// Your code goes here
	FloatImage output(im);  // CHANGEME

	return output;
}


// Histograms

// Stretches the pixel values of channel c of im so that the minimum value maps to 0,
// and the maximum value maps to 1
void autoLevels(FloatImage &im, int c)
{
	// stretch pixel values of channel c to fill 0..1 range
	// you may find FloatImage::min and FloatImage::max useful
}

// initialize a histogram with numBins bins from the c-th channel of im
Histogram::Histogram(const FloatImage &im, int c, int numBins) :
	m_pdf(numBins, 0.0f), m_cdf(numBins, 0.0f)
{
	// populate m_pdf with histogram values

	// Grad/extra credit: populate m_cdf with running sum of m_pdf
}

// return the histogram bin that value falls within
int Histogram::inverseCDF(double value) const
{
	return 0; // CHANGEME
}

// Produce a numBins() x 100 x 3 image containing a visualization of the
// red, green and blue histogram pdfs
FloatImage visualizeRGBHistogram(const Histogram &histR,
                                 const Histogram &histG,
                                 const Histogram &histB)
{
	assert(histR.numBins() == histG.numBins() && histR.numBins() == histB.numBins());

	// create an image of appropriate size
	FloatImage im(histR.numBins(), 100, 3);

	// populate im with RGB histograms
	return im;
}

// Return a FloatImage which is the result of applying histogram equalization to im
FloatImage equalizeRGBHistograms(const FloatImage &im)
{
	int numLevels = 256;

	// create an image of appropriate size
	FloatImage output = im;

	// apply histogram equalization to each channel

	return output;
}

// Return a FloatImage which is the result of transfering the histogram of F2 onto the image data in F1
FloatImage matchRGBHistograms(const FloatImage &F1, const FloatImage &F2)
{
	int numBins = 256;

	FloatImage output = F1;

	// perform histogram matching

	return output;
}

