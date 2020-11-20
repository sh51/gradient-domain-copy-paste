#include "filtering.h"
#include "utils.h"

using namespace std;

// test the smart accessor function
void testSmartAccessor()
{
    // load an image and create 2 images that will test the smart accessor
    const FloatImage input(DATA_DIR "/input/a4/bear.png");
    input.write(DATA_DIR "/output/bear.png");

    FloatImage clampTrue(input.width(), input.height(), input.channels());
    FloatImage clampFalse(input.width(), input.height(), input.channels());
    
    for (int z=0; z<input.channels(); z++)
    {
        for (int x=0; x<input.width(); x++)
        {
            for (int y=0; y<input.height(); y++)
            {
                // replace non-valid pixel values with the value of the nearest pixel
                clampTrue(x,y,z) = input.smartAccessor(x-10, y-10, z, true);
                // replace non-valid pixel values with black (value=0)
                clampFalse(x,y,z) = input.smartAccessor(x-10, y-10, z);
            }
        }
    }
    
    clampTrue.write(DATA_DIR "/output/smartAccessor_clampTrue.png");
    clampFalse.write(DATA_DIR "/output/smartAccessor_clampFalse.png");
}

// test your box blurring fuctions
void testBoxBlur() {
    
    FloatImage im(DATA_DIR "/input/a4/Cambridge2.png");
    int k=9; //set the width of your box to be 9 pixels (must be odd)
    
    // blur your image by doing a moving average
    FloatImage boxBlur1 = boxBlur(im, k);
    boxBlur1.write(DATA_DIR "/output/boxBlurred1.png");
    
    // blur your image in the same way using the Filter class
    FloatImage boxBlur2 = boxBlur_filterClass(im, k);
    boxBlur2.write(DATA_DIR "/output/boxBlurred2.png");
    
    // verify that your two funtion implementations obtain the
    // same result for the same value of k
    FloatImage diffImg = (boxBlur1 - boxBlur2)/2 + 0.5;
    diffImg.write(DATA_DIR "/output/boxBlurDiff.png");
    
}

// test that your convolution is properly flipping your kernel
void testShiftedImpulse()
{
    FloatImage im(DATA_DIR "/input/a4/pru.png");
    
    // create an impulse kernel shifted 1 pixel to the right and down
    // [ 0 0 0 ]
    // [ 0 0 0 ]
    // [ 0 0 1 ]
    float k = 5;
    Filter impulse(k, k);
    impulse(k-1,k-2) = 1.0;
    
    // filter your image with this impulse kernel. Remmber that in convolution
    // you must flip your kernel. Thus, this should result in a copy of your image
    // shifted 1 pixel to the right and down.
    FloatImage imFilter = impulse.Convolve(im);
    imFilter.write(DATA_DIR "/output/impulseFiltered.png");
    
}

// test functions that extract gradient information from your image
void testGradient()
{
    // create Sobel Filter that extracts horizontal gradients
    // [ -1 0 1 ]
    // [ -2 0 2 ]
    // [ -1 0 1 ]
    Filter sobelX(3, 3);
    sobelX(0,0) = -1.0; sobelX(1,0) = 0.0; sobelX(2,0) = 1.0;
    sobelX(0,1) = -2.0; sobelX(1,1) = 0.0; sobelX(2,1) = 2.0;
    sobelX(0,2) = -1.0; sobelX(1,2) = 0.0; sobelX(2,2) = 1.0;
    
    // verify that your filter is correct by using it to filter an impulse image
    FloatImage impulse = impulseImg(11); //create an image containing an impulse
    // convolve the impulse image with the Sobel kernel. We divide the output by 4 and
    // add 0.5 to make the range of the image between 0 and 1
    FloatImage verifyKernel = sobelX.Convolve(impulse)/4 + 0.5;
    verifyKernel.write(DATA_DIR "/output/verifySobelKernel.png");
    
    // filter an image using the sobel kernel
    FloatImage im(DATA_DIR "/input/a4/lounge_view.png");
    FloatImage sobelFiltered = sobelX.Convolve(im,true);

    FloatImage im2(DATA_DIR "/input/a4/dart.png");
    FloatImage sobelFiltered2 = sobelX.Convolve(im2,true);
    
    // make the range of the output image from 0 to 1
    // since the Sobel filter changes the range of a (0,1) image to (-2,2)
    FloatImage sobelOut = sobelFiltered/4 + 0.5;
    sobelOut.write(DATA_DIR "/output/sobelFiltered.png");

    FloatImage sobelOut2 = sobelFiltered2/4 + 0.5;
    sobelOut2.write(DATA_DIR "/output/sobelFiltered2.png");


    // compute the magnitude of the image using Sobel filters
    // that extract horizontal and vertical gradients
    FloatImage imMag = gradientMagnitude(im);
    imMag.write(DATA_DIR "/output/imageMagnitude.png");

    FloatImage imMag2 = gradientMagnitude(im2);
    imMag2.write(DATA_DIR "/output/imageMagnitude2.png");
}

// test filtering by a Gaussian kernel
void testGaussianFilters()
{
    float sigma = 3.0; //set the standard deviation of the Gaussians
    Timer timer;
    
    FloatImage im(DATA_DIR "/input/a4/Cambridge2.png");
    
    // blur an image in just the X direction
    FloatImage blurHorizontal = gaussianBlur_horizontal(im, sigma);
    blurHorizontal.write(DATA_DIR "/output/gaussBlurHorizontal.png");
    
    // blur an image in 2D using a full 2D kernel
    timer.reset();
    FloatImage blur2D = gaussianBlur_2D(im, sigma);
    // print the time it takes to run this function
    printf("Filtering with 2D Gaussian kernel took %3.5f seconds\n", timer.elapsed()/1000.f);
    blur2D.write(DATA_DIR "/output/gaussBlur2D.png");

    // blur an image in 2D using 2 1D Gaussian kernels
    timer.reset();
    FloatImage blur2DSeperable = gaussianBlur_seperable(im, sigma);
    // print the time it takes to run this function
    printf("Filtering using seperable Gaussian kernels took %3.5f seconds\n", timer.elapsed()/1000.f);
    blur2DSeperable.write(DATA_DIR "/output/gaussBlur2D_seperable.png");
    
    // verify that both methods result in the same image
    FloatImage diffImg = (blur2D - blur2DSeperable)/2 + 0.5;
    diffImg.write(DATA_DIR "/output/gaussBlurDiff.png");
}

// test sharpening function
void testSharpen()
{
    FloatImage im(DATA_DIR "/input/a4/Cambridge1.png");
    
    FloatImage sharp = unsharpMask(im,2.0,3.0,1.5);
    sharp.write(DATA_DIR "/output/sharp.png");
}

// test bilaterial filtering
void testBilaterial()
{
    FloatImage im(DATA_DIR "/input/a4/lens.png");
    
    // Perform bilaterial filtering on an RGB image
    FloatImage rgbBilatIm = bilateral(im);
    rgbBilatIm.write(DATA_DIR "/output/bilaterial_RGB.png");
    
    // NOTE: Uncomment the code below if you have implemented bilaYUV
    // Perform bilaterial filtering with different domain sigmas for Y and UV
//    FloatImage yuvBilatIm = bilaYUV(im);
//    yuvBilatIm.write(DATA_DIR "/output/bilaterial_YUV.png");
    
}


// This is a way for you to test your functions.
// We will not grade the contents of the main function
int main()
{
    // uncomment these test functions as you complete the assignment
    // try {testSmartAccessor();}      catch(...) { cout << "EXCEPTION: Smart Accessor failed" << endl; }
    //try {testBoxBlur(); }           catch(...) { cout << "EXCEPTION: Box Blur failed" << endl; }
    //try {testShiftedImpulse();}     catch(...) { cout << "EXCEPTION: Box Blur failed" << endl; }
    //try {testGradient();}           catch(...) { cout << "EXCEPTION: Box Blur failed" << endl; }
    //try {testGaussianFilters();}    catch(...) { cout << "EXCEPTION: Box Blur failed" << endl; }
    //try {testSharpen();}            catch(...) { cout << "EXCEPTION: Box Blur failed" << endl; }
    //try {testBilaterial();}         catch(...) { cout << "EXCEPTION: Box Blur failed" << endl; }
}
