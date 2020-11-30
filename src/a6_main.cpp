#include "morphing.h"
#include "a6.h"
#include "a2.h"
#include <map>
#include "utils.h"
#include "Eigen/Sparse"
typedef Eigen::SparseMatrix<double> SpMatrix;
typedef Eigen::Triplet<double> Triplet;
typedef Eigen::VectorXd VectorD;


using namespace std;

// min non-zero pixel value of image
float image_minnonzero(const FloatImage &im)
{
	float minf = 2;
	long imsize = im.size();
	for (int i = 0; i < imsize; i++)
		if (im(i) > 0)
			minf = min(minf, im(i));

	return minf;
}

// image --> log10FloatImage
FloatImage log10FloatImage(const FloatImage &im)
{
	// Taking a linear image im, transform to log10 scale.
	// To avoid infinity issues, make any 0-valued pixel be equal the the minimum
	// non-zero value. See image_minnonzero(im).
	float minf = image_minnonzero(im);
	FloatImage loglumi(im.width(), im.height(), im.channels());
	for (int i = 0; i < im.size(); i++)
		loglumi(i) = log10(max(minf, im(i)));
	return loglumi;
}

// FloatImage --> 10^FloatImage
FloatImage exp10FloatImage(const FloatImage &im)
{
	// take an image in log10 domain and transform it back to linear domain.
	// see pow(a, b)
	FloatImage out=im;
	for(int i=0;i<im.size();i++){
		out(i)=pow(10,im(i));

	}

	return out; // change this
}
// test the smart accessor function
void testSmartAccessor()
{
	// load an image and create 2 images that will test the smart accessor
	const FloatImage input(DATA_DIR "/input/a6/bear.png");
	input.write(DATA_DIR "/output/bear.png");

	FloatImage clampTrue(input.width(), input.height(), input.channels());
	FloatImage clampFalse(input.width(), input.height(), input.channels());

	for (int z = 0; z < input.channels(); z++)
	{
		for (int x = 0; x < input.width(); x++)
		{
			for (int y = 0; y < input.height(); y++)
			{
				// replace non-valid pixel values with the value of the nearest pixel
				clampTrue(x, y, z) = input.smartAccessor(x - 10, y - 10, z, true);
				// replace non-valid pixel values with black (value=0)
				clampFalse(x, y, z) = input.smartAccessor(x - 10, y - 10, z);
			}
		}
	}

	clampTrue.write(DATA_DIR "/output/smartAccessor_clampTrue.png");
	clampFalse.write(DATA_DIR "/output/smartAccessor_clampFalse.png");
}

// a function to test scaling
void testScaling()
{
	// load in the image and print out original size
	float fact = 2.253;
	const FloatImage bostonim(DATA_DIR "/input/a6/BostonRainbow-crop-400.png");
	printf("Boston image is %dx%dx%d\n", bostonim.width(), bostonim.height(), bostonim.channels());

	// scale using NN interpolation and print the size of the new image
	FloatImage scaledNN = scaleNN(bostonim, fact);
	scaledNN.write(DATA_DIR "/output/BostonRainbow-scaled-NN.png");
	printf("Scaled-NN image is %dx%dx%d\n", scaledNN.width(), scaledNN.height(), scaledNN.channels());

	// scale using bilinear interpolation and print the size of the new image
	FloatImage scaledLin = scaleLin(bostonim, fact);
	scaledLin.write(DATA_DIR "/output/BostonRainbow-scaled-Lin.png");
	printf("Scaled-Lin image is %dx%dx%d\n", scaledLin.width(), scaledLin.height(), scaledLin.channels());
}

// a function to test rotation for those who have done it
void testRotation()
{
	float theta = 3.141 / 4;

	const FloatImage bostonim(DATA_DIR "/input/a6/BostonRainbow-crop-400.png");

	FloatImage rot = rotate(bostonim, theta);
	rot.write(DATA_DIR "/output/BostonRainbow-rotated.png");
	printf("Scaled-Lin image is %dx%dx%d\n", rot.width(), rot.height(), rot.channels());
}

// test warp by 1
void testWarpBy1()
{
	FloatImage bearim(DATA_DIR "/input/a6/bear.png");

	// define before and after segments
	Segment segBefore(0, 0, 10, 0);
	Segment segAfter(10, 10, 30, 15);

	FloatImage warp1 = warpBy1(bearim, segBefore, segAfter);

	warp1.write(DATA_DIR "/output/warpBy1.png");
}

// a function to test your morphing function with the fredo and werewolf images
void testMorph()
{
	// load the images
	FloatImage fredo(DATA_DIR "/input/a6/fredo.png");
	FloatImage werewolf(DATA_DIR "/input/a6/werewolf.png");

	// print out the size of the two images
	printf("Fredo image is %dx%dx%d\n", fredo.width(), fredo.height(), fredo.channels());
	printf("Werewolf image is %dx%dx%d\n", werewolf.width(), werewolf.height(), werewolf.channels());

	// define the segments corresponding to fredo (segsBefore) and the werewolf (segsAfter)
	vector<Segment> segsBefore, segsAfter;

	segsBefore.push_back(Segment(87, 128, 109, 127));
	segsBefore.push_back(Segment(143, 127, 162, 131));
	segsBefore.push_back(Segment(96, 197, 129, 190));
	segsBefore.push_back(Segment(118, 221, 132, 200));
	segsBefore.push_back(Segment(140, 238, 165, 170));
	segsBefore.push_back(Segment(71, 242, 44, 196));
	segsBefore.push_back(Segment(9, 46, 34, 14));
	segsBefore.push_back(Segment(175, 66, 154, 27));

	segsAfter.push_back(Segment(83, 114, 107, 109));
	segsAfter.push_back(Segment(139, 103, 157, 101));
	segsAfter.push_back(Segment(100, 170, 132, 151));
	segsAfter.push_back(Segment(125, 198, 145, 171));
	segsAfter.push_back(Segment(142, 196, 158, 139));
	segsAfter.push_back(Segment(96, 211, 75, 180));
	segsAfter.push_back(Segment(11, 41, 33, 7));
	segsAfter.push_back(Segment(175, 56, 155, 13));

	int numFrames = 10;
	//numFrames=1;
	// create an image morphing between fredo and werewolf at time t=0.5
	vector<FloatImage> imMorph = morph(fredo, werewolf, segsBefore, segsAfter, numFrames);

	// write out images
	char buffer[255];
	for (int n = 0; n < numFrames + 2; n++)
	{
		FloatImage im = imMorph[n];
		sprintf(buffer, DATA_DIR "/output/fredo_werewolf_morph_%d.png", n);

		im.write(buffer);
	}
}

void classMorph(int studentNumber) {
	char buffer1[255];
	char buffer2[255];

	int studentNumberTwo = studentNumber + 1;

	if (studentNumber == 34)  studentNumberTwo = 0;

	sprintf(buffer1, DATA_DIR "/input/a6/class-morph/class-morph-%02d.jpg", studentNumber);
	sprintf(buffer2, DATA_DIR "/input/a6/class-morph/class-morph-%02d.jpg", studentNumberTwo);

	FloatImage start(buffer1);
	FloatImage end(buffer2);

	vector<Segment> segsBefore, segsAfter;

	// TODO :: REPLACE WITH YOUR NEW SEGMENTS HERE
	segsBefore.push_back(Segment(235, 225, 183, 201)); 
segsBefore.push_back(Segment(298, 230, 351, 205)); 
segsBefore.push_back(Segment(161, 218, 179, 203)); 
segsBefore.push_back(Segment(376, 219, 355, 208)); 
segsBefore.push_back(Segment(269, 252, 266, 346)); 
segsBefore.push_back(Segment(273, 343, 291, 338)); 
segsBefore.push_back(Segment(258, 340, 243, 339)); 
segsBefore.push_back(Segment(230, 327, 208, 354)); 
segsBefore.push_back(Segment(302, 327, 325, 352)); 
segsBefore.push_back(Segment(252, 407, 280, 406)); 
segsBefore.push_back(Segment(219, 380, 245, 404)); 
segsBefore.push_back(Segment(310, 379, 285, 403)); 
segsBefore.push_back(Segment(267, 370, 254, 368)); 
segsBefore.push_back(Segment(269, 372, 279, 369)); 
segsBefore.push_back(Segment(222, 378, 249, 369)); 
segsBefore.push_back(Segment(283, 368, 307, 378)); 
segsBefore.push_back(Segment(235, 450, 292, 451)); 
segsBefore.push_back(Segment(221, 101, 312, 99)); 
segsBefore.push_back(Segment(315, 101, 390, 219)); 
segsBefore.push_back(Segment(217, 104, 163, 202)); 
segsBefore.push_back(Segment(389, 226, 358, 391)); 
segsBefore.push_back(Segment(354, 394, 299, 447)); 
segsBefore.push_back(Segment(170, 391, 231, 446)); 
segsBefore.push_back(Segment(158, 247, 168, 383)); 

segsBefore.push_back(Segment(223, 378, 265, 383)); 
segsBefore.push_back(Segment(270, 384, 302, 379)); 
segsBefore.push_back(Segment(225, 381, 266, 389)); 
segsBefore.push_back(Segment(269, 390, 302, 383)); 

segsBefore.push_back(Segment(304, 260, 346, 251)); 
segsBefore.push_back(Segment(231, 258, 185, 248)); 
//verical eye
segsBefore.push_back(Segment(210, 241, 211, 257)); 
segsBefore.push_back(Segment(324, 244, 324, 259)); 



segsAfter.push_back(Segment(212, 247, 161, 241)); 
segsAfter.push_back(Segment(257, 247, 311, 237)); 
segsAfter.push_back(Segment(148, 252, 157, 242)); 
segsAfter.push_back(Segment(330, 252, 313, 239)); 
// segsAfter.push_back(Segment(207, 277, 171, 275)); //eye
// segsAfter.push_back(Segment(266, 274, 304, 271)); //eye
segsAfter.push_back(Segment(233, 272, 234, 337)); 
segsAfter.push_back(Segment(244, 333, 260, 328)); 
segsAfter.push_back(Segment(228, 335, 211, 329)); 
segsAfter.push_back(Segment(203, 319, 172, 350)); 
segsAfter.push_back(Segment(267, 312, 301, 341)); 
segsAfter.push_back(Segment(226, 395, 255, 396)); 
segsAfter.push_back(Segment(185, 361, 220, 393)); 
segsAfter.push_back(Segment(295, 359, 260, 394)); 
segsAfter.push_back(Segment(237, 349, 221, 346)); 
segsAfter.push_back(Segment(240, 350, 254, 346)); 
segsAfter.push_back(Segment(187, 359, 218, 347)); 
segsAfter.push_back(Segment(257, 345, 292, 355)); 
segsAfter.push_back(Segment(214, 435, 279, 433)); 
segsAfter.push_back(Segment(183, 151, 263, 142)); 
segsAfter.push_back(Segment(266, 145, 345, 249)); 
segsAfter.push_back(Segment(179, 155, 146, 239)); 
segsAfter.push_back(Segment(346, 255, 329, 391)); 
segsAfter.push_back(Segment(323, 393, 285, 430)); 
segsAfter.push_back(Segment(157, 383, 210, 432)); 
segsAfter.push_back(Segment(147, 271, 155, 376)); 

segsAfter.push_back(Segment(186, 359, 239, 359)); 
segsAfter.push_back(Segment(243, 358, 292, 358)); 
segsAfter.push_back(Segment(189, 363, 238, 382)); 
segsAfter.push_back(Segment(242, 382, 293, 361)); 

segsAfter.push_back(Segment(269, 275, 303, 270)); 
segsAfter.push_back(Segment(204, 275, 171, 274));
//eye vertical 
segsAfter.push_back(Segment(188, 267, 189, 276)); 
segsAfter.push_back(Segment(289, 264, 290, 275)); 




	int numFrames = 19;

	vector<FloatImage> imMorph = morph(start, end, segsBefore, segsAfter, numFrames);


	// write out images
	char buffer[255];
	for (int n = 0; n < 20; n++)
	{
		FloatImage im = imMorph[n];
		int classIndex = n + 20 * studentNumber;

		sprintf(buffer, DATA_DIR "/output/class-morph%03d.jpg", classIndex);

		im.write(buffer);
	}
}

//mixed blending
void testMixedBlend(){

	//load images	
	FloatImage sourceImage(DATA_DIR "/input/a6/source_r.png");
	FloatImage targetImage(DATA_DIR "/input/a6/target_r.png");
	FloatImage maskImage(DATA_DIR "/input/a6/mask_r.png");
	FloatImage targetImage_copy=targetImage;
	// sourceImage=log10FloatImage(sourceImage);
	// targetImage=log10FloatImage(targetImage);
	


	//check if three images have the same dimension


	//save every pixel inside the mask into a map
	//map(pixel,int number)
	map<unsigned int, unsigned int> varMap;
	{
		int number=0;
		for(int x=0;x<maskImage.width();x++){
			for(int y=0;y<maskImage.height();y++){
				//white part of the mask
				if(maskImage(x,y,0)>0.9f){
					varMap[y*maskImage.width()+x]=number;
					number++;
				}
			}
		}
	}


	//construct sparse matrix M in Mx = b
	vector<Triplet> mt;
	{			
		int i=0;		

		for(int x=0;x<maskImage.width();x++){
			for(int y=0;y<maskImage.height();y++){
				
				if(maskImage(x,y,0)>0.9f){
					int j=y*maskImage.width()+x;
					mt.push_back(Triplet(i,varMap[j],4));
				
					if(maskImage(x+1,y,0)>0.9f){
						mt.push_back(Triplet(i,varMap[j+1],-1));
					}
					if(maskImage(x-1,y,0)>0.9f){
						mt.push_back(Triplet(i,varMap[j-1],-1));
					}
					if(maskImage(x,y+1,0)>0.9f){
						mt.push_back(Triplet(i,varMap[j+maskImage.width()],-1));
					}
					if(maskImage(x,y-1,0)>0.9f){
						mt.push_back(Triplet(i,varMap[j-maskImage.width()],-1));
					}
					i++;
				}

			}
		}	

	}
	
	

	Eigen::SimplicialCholesky<SpMatrix> solver;
	{
		SpMatrix mat(varMap.size(),varMap.size());
		mat.setFromTriplets(mt.begin(),mt.end());
		//cout<<mat<<endl;
		solver.compute(mat);		
	}
	


	//construct b
	VectorD b(varMap.size());
	VectorD sol[3];
	for(int ic=0;ic<3;ic++){		
		int i=0;
		for(int x=0;x<maskImage.width();x++){
			for(int y=0;y<maskImage.height();y++){


				if(maskImage(x,y,ic)>0.9f){
					//graident from source image
					float pixel=0.f;
					// float pixel_s=4*sourceImage(x,y,ic)-sourceImage(x-1,y,ic)-sourceImage(x+1,y,ic)
					// 			-sourceImage(x,y-1,ic)-sourceImage(x,y+1,ic);		
					//pixel=pixel_t;			
					if(abs(sourceImage(x,y,ic)-sourceImage(x-1,y,ic))<abs(targetImage(x,y,ic)-targetImage(x-1,y,ic))){
						pixel+=targetImage(x,y,ic)-targetImage(x-1,y,ic);
					}else{
						pixel+=sourceImage(x,y,ic)-sourceImage(x-1,y,ic);
					}

					if(abs(sourceImage(x,y,ic)-sourceImage(x+1,y,ic))<abs(targetImage(x,y,ic)-targetImage(x+1,y,ic))){
						pixel+=targetImage(x,y,ic)-targetImage(x+1,y,ic);
					}else{
						pixel+=sourceImage(x,y,ic)-sourceImage(x+1,y,ic);
					}

					if(abs(sourceImage(x,y,ic)-sourceImage(x,y-1,ic))<abs(targetImage(x,y,ic)-targetImage(x,y-1,ic))){
						pixel+=targetImage(x,y,ic)-targetImage(x,y-1,ic);
					}else{
						pixel+=sourceImage(x,y,ic)-sourceImage(x,y-1,ic);
					}

					if(abs(sourceImage(x,y,ic)-sourceImage(x,y+1,ic))<abs(targetImage(x,y,ic)-targetImage(x,y+1,ic))){
						pixel+=targetImage(x,y,ic)-targetImage(x,y+1,ic);
					}else{
						pixel+=sourceImage(x,y,ic)-sourceImage(x,y+1,ic);
					}
					

					//add boundary from target image
					if(maskImage(x+1,y,ic)<0.1f){
						pixel+=targetImage(x+1,y,ic);
					}
					if(maskImage(x-1,y,ic)<0.1f){
						pixel+=targetImage(x-1,y,ic);
					}
					if(maskImage(x,y+1,ic)<0.1f){
						pixel+=targetImage(x,y+1,ic);
					}
					if(maskImage(x,y-1,ic)<0.1f){
						pixel+=targetImage(x,y-1,ic);
					}
					b[i]=pixel;
					i++;
				}


			}
		}	
		//solve x and save each channel in sol
		sol[ic]=solver.solve(b);
	}


	//write image to output
	FloatImage output(maskImage.width(),maskImage.height(),3);
	
	for (int ic = 0; ic < 3; ic++){
		for (int x = 0; x < maskImage.width(); x++){
			for (int y = 0; y < maskImage.height(); y++){
				if(maskImage(x,y,0)<0.1f){
					output(x,y,ic)=targetImage_copy(x,y,ic);
				}else{
					int j=y*maskImage.width()+x;
					//float value=(float)pow(10,sol[ic][varMap[j]]);
					float value=(float)sol[ic][varMap[j]];
					//output(x,y,ic)=clamp(value,0.f,1.f);
					output(x,y,ic)=value;
					//output(x,y,ic)=sourceImage(x,y,ic);
				}				
			}			
		}		
	}
	
	output.write(DATA_DIR "/output/linear_r_mixed.png");
}


//flatten local texture
void testFlatten(){

	//load images	
	FloatImage sourceImage(DATA_DIR "/input/a6/source_b.png");//baby img
	//FloatImage targetImage(DATA_DIR "/input/a6/target_b.png");//edge detector
	FloatImage maskImage(DATA_DIR "/input/a6/mask_b.png");
	FloatImage targetImage_copy=sourceImage;
	//sourceImage=log10FloatImage(sourceImage);
	


	//check if three images have the same dimension


	//save every pixel inside the mask into a map
	//map(pixel,int number)
	map<unsigned int, unsigned int> varMap;
	{
		int number=0;
		for(int x=0;x<maskImage.width();x++){
			for(int y=0;y<maskImage.height();y++){
				//white part of the mask
				if(maskImage(x,y,0)>0.9f){
					varMap[y*maskImage.width()+x]=number;
					number++;
				}
			}
		}
	}


	//construct sparse matrix M in Mx = b
	vector<Triplet> mt;
	{			
		int i=0;		

		for(int x=0;x<maskImage.width();x++){
			for(int y=0;y<maskImage.height();y++){
				
				if(maskImage(x,y,0)>0.9f){
					int j=y*maskImage.width()+x;
					mt.push_back(Triplet(i,varMap[j],4));
				
					if(maskImage(x+1,y,0)>0.9f){
						mt.push_back(Triplet(i,varMap[j+1],-1));
					}
					if(maskImage(x-1,y,0)>0.9f){
						mt.push_back(Triplet(i,varMap[j-1],-1));
					}
					if(maskImage(x,y+1,0)>0.9f){
						mt.push_back(Triplet(i,varMap[j+maskImage.width()],-1));
					}
					if(maskImage(x,y-1,0)>0.9f){
						mt.push_back(Triplet(i,varMap[j-maskImage.width()],-1));
					}
					i++;
				}

			}
		}	

	}
	
	

	Eigen::SimplicialCholesky<SpMatrix> solver;
	{
		SpMatrix mat(varMap.size(),varMap.size());
		mat.setFromTriplets(mt.begin(),mt.end());
		//cout<<mat<<endl;
		solver.compute(mat);		
	}
	


	//construct b
	VectorD b(varMap.size());
	VectorD sol[3];
	for(int ic=0;ic<3;ic++){		
		int i=0;
		for(int x=0;x<maskImage.width();x++){
			for(int y=0;y<maskImage.height();y++){


				if(maskImage(x,y,ic)>0.9f){
					//graident from source image
					float pixel=0.f;
					///////////////////////////convert to grey scale?
					//0.3 0.59 0.11
					if(0.3*abs(sourceImage(x,y,0)-sourceImage(x-1,y,0))
						+0.59*abs(sourceImage(x,y,1)-sourceImage(x-1,y,1))
						+0.11*abs(sourceImage(x,y,2)-sourceImage(x-1,y,2))>0.006f){
						pixel+=sourceImage(x,y,ic)-sourceImage(x-1,y,ic);
					}
					if(0.3*abs(sourceImage(x,y,0)-sourceImage(x+1,y,0))
						+0.59*abs(sourceImage(x,y,1)-sourceImage(x+1,y,1))
						+0.11*abs(sourceImage(x,y,2)-sourceImage(x+1,y,2))>0.006f){
						pixel+=sourceImage(x,y,ic)-sourceImage(x+1,y,ic);
					}
					if(0.3*abs(sourceImage(x,y,0)-sourceImage(x,y-1,0))
						+0.59*abs(sourceImage(x,y,1)-sourceImage(x,y-1,1))
						+0.11*abs(sourceImage(x,y,2)-sourceImage(x,y-1,2))>0.006f){
						pixel+=sourceImage(x,y,ic)-sourceImage(x,y-1,ic);
					}
					if(0.3*abs(sourceImage(x,y,0)-sourceImage(x,y+1,0))
						+0.59*abs(sourceImage(x,y,1)-sourceImage(x,y+1,1))
						+0.11*abs(sourceImage(x,y,2)-sourceImage(x,y+1,2))>0.006f){
						pixel+=sourceImage(x,y,ic)-sourceImage(x,y+1,ic);
					}

					//add boundary from target image
					if(maskImage(x+1,y,ic)<0.1f){
						pixel+=sourceImage(x+1,y,ic);
					}
					if(maskImage(x-1,y,ic)<0.1f){
						pixel+=sourceImage(x-1,y,ic);
					}
					if(maskImage(x,y+1,ic)<0.1f){
						pixel+=sourceImage(x,y+1,ic);
					}
					if(maskImage(x,y-1,ic)<0.1f){
						pixel+=sourceImage(x,y-1,ic);
					}
					b[i]=pixel;
					i++;
				}


			}
		}	
		//solve x and save each channel in sol
		sol[ic]=solver.solve(b);
	}


	//write image to output
	FloatImage output(maskImage.width(),maskImage.height(),3);
	
	for (int ic = 0; ic < 3; ic++){
		for (int x = 0; x < maskImage.width(); x++){
			for (int y = 0; y < maskImage.height(); y++){
				if(maskImage(x,y,0)<0.1f){
					output(x,y,ic)=targetImage_copy(x,y,ic);
				}else{
					int j=y*maskImage.width()+x;
					//float value=(float)pow(10,sol[ic][varMap[j]]);
					float value=(float)sol[ic][varMap[j]];
					//output(x,y,ic)=clamp(value,0.f,1.f);
					output(x,y,ic)=value;
					//output(x,y,ic)=sourceImage(x,y,ic);
				}				
			}			
		}		
	}
	
	output.write(DATA_DIR "/output/flatten_linear_b.png");
}


//change local color
void testColor(){

	//load images	
	FloatImage targetImage(DATA_DIR "/input/a6/source_cat.png");
	//FloatImage targetImage(DATA_DIR "/input/a6/target_f.png");
	FloatImage sourceImage=targetImage;
	for(int i=0;i<sourceImage.width();i++){
		for(int j=0;j<sourceImage.height();j++){
			sourceImage(i,j,0)*=0.5;
			sourceImage(i,j,1)*=1.5;
			sourceImage(i,j,2)*=1.0;
		}
	}
	
	FloatImage maskImage(DATA_DIR "/input/a6/cat_mask.png");
	FloatImage targetImage_copy=targetImage;
	// sourceImage=log10FloatImage(sourceImage);
	// targetImage=log10FloatImage(targetImage);
	


	//check if three images have the same dimension


	//save every pixel inside the mask into a map
	//map(pixel,int number)
	map<unsigned int, unsigned int> varMap;
	{
		int number=0;
		for(int x=0;x<maskImage.width();x++){
			for(int y=0;y<maskImage.height();y++){
				//white part of the mask
				if(maskImage(x,y,0)>0.9f){
					varMap[y*maskImage.width()+x]=number;
					number++;
				}
			}
		}
	}


	//construct sparse matrix M in Mx = b
	vector<Triplet> mt;
	{			
		int i=0;		

		for(int x=0;x<maskImage.width();x++){
			for(int y=0;y<maskImage.height();y++){
				
				if(maskImage(x,y,0)>0.9f){
					int j=y*maskImage.width()+x;
					mt.push_back(Triplet(i,varMap[j],4));
				
					if(maskImage(x+1,y,0)>0.9f){
						mt.push_back(Triplet(i,varMap[j+1],-1));
					}
					if(maskImage(x-1,y,0)>0.9f){
						mt.push_back(Triplet(i,varMap[j-1],-1));
					}
					if(maskImage(x,y+1,0)>0.9f){
						mt.push_back(Triplet(i,varMap[j+maskImage.width()],-1));
					}
					if(maskImage(x,y-1,0)>0.9f){
						mt.push_back(Triplet(i,varMap[j-maskImage.width()],-1));
					}
					i++;
				}

			}
		}	

	}
	
	

	Eigen::SimplicialCholesky<SpMatrix> solver;
	{
		SpMatrix mat(varMap.size(),varMap.size());
		mat.setFromTriplets(mt.begin(),mt.end());
		//cout<<mat<<endl;
		solver.compute(mat);		
	}
	


	//construct b
	VectorD b(varMap.size());
	VectorD sol[3];
	for(int ic=0;ic<3;ic++){		
		int i=0;
		for(int x=0;x<maskImage.width();x++){
			for(int y=0;y<maskImage.height();y++){


				if(maskImage(x,y,ic)>0.9f){
					//graident from source image
					float pixel=4*sourceImage(x,y,ic)-sourceImage(x-1,y,ic)-sourceImage(x+1,y,ic)
								-sourceImage(x,y-1,ic)-sourceImage(x,y+1,ic);

					//add boundary from target image
					if(maskImage(x+1,y,ic)<0.1f){
						pixel+=targetImage(x+1,y,ic);
					}
					if(maskImage(x-1,y,ic)<0.1f){
						pixel+=targetImage(x-1,y,ic);
					}
					if(maskImage(x,y+1,ic)<0.1f){
						pixel+=targetImage(x,y+1,ic);
					}
					if(maskImage(x,y-1,ic)<0.1f){
						pixel+=targetImage(x,y-1,ic);
					}
					b[i]=pixel;
					i++;
				}


			}
		}	
		//solve x and save each channel in sol
		sol[ic]=solver.solve(b);
	}


	//write image to output
	FloatImage output(maskImage.width(),maskImage.height(),3);
	
	for (int ic = 0; ic < 3; ic++){
		for (int x = 0; x < maskImage.width(); x++){
			for (int y = 0; y < maskImage.height(); y++){
				if(maskImage(x,y,0)<0.1f){
					output(x,y,ic)=targetImage_copy(x,y,ic);
				}else{
					int j=y*maskImage.width()+x;
					//float value=(float)pow(10,sol[ic][varMap[j]]);
					float value=(float)sol[ic][varMap[j]];
					//output(x,y,ic)=clamp(value,0.f,1.f);
					output(x,y,ic)=value;
					//output(x,y,ic)=sourceImage(x,y,ic);
				}				
			}			
		}		
	}
	
	output.write(DATA_DIR "/output/change_color_cat.png");
}


//change local lumination
void testLumi(){

	//load images	
	// FloatImage sourceImage(DATA_DIR "/input/a6/orange.png");
	// FloatImage targetImage(DATA_DIR "/input/a6/orange.png");
	// FloatImage maskImage(DATA_DIR "/input/a6/mask_orange.png");
	FloatImage sourceImage(DATA_DIR "/input/a6/man.png");
	FloatImage targetImage(DATA_DIR "/input/a6/man.png");
	FloatImage maskImage(DATA_DIR "/input/a6/man_mask.png");
	FloatImage targetImage_copy=targetImage;
	sourceImage=log10FloatImage(sourceImage);
	targetImage=log10FloatImage(targetImage);
	


	//check if three images have the same dimension


	//save every pixel inside the mask into a map
	//map(pixel,int number)
	map<unsigned int, unsigned int> varMap;
	{
		int number=0;
		for(int x=0;x<maskImage.width();x++){
			for(int y=0;y<maskImage.height();y++){
				//white part of the mask
				if(maskImage(x,y,0)>0.9f){
					varMap[y*maskImage.width()+x]=number;
					number++;
				}
			}
		}
	}


	//construct sparse matrix M in Mx = b
	vector<Triplet> mt;
	{			
		int i=0;		

		for(int x=0;x<maskImage.width();x++){
			for(int y=0;y<maskImage.height();y++){
				
				if(maskImage(x,y,0)>0.9f){
					int j=y*maskImage.width()+x;
					mt.push_back(Triplet(i,varMap[j],4));
				
					if(maskImage(x+1,y,0)>0.9f){
						mt.push_back(Triplet(i,varMap[j+1],-1));
					}
					if(maskImage(x-1,y,0)>0.9f){
						mt.push_back(Triplet(i,varMap[j-1],-1));
					}
					if(maskImage(x,y+1,0)>0.9f){
						mt.push_back(Triplet(i,varMap[j+maskImage.width()],-1));
					}
					if(maskImage(x,y-1,0)>0.9f){
						mt.push_back(Triplet(i,varMap[j-maskImage.width()],-1));
					}
					i++;
				}

			}
		}	

	}
	
	

	Eigen::SimplicialCholesky<SpMatrix> solver;
	{
		SpMatrix mat(varMap.size(),varMap.size());
		mat.setFromTriplets(mt.begin(),mt.end());
		//cout<<mat<<endl;
		solver.compute(mat);		
	}
	


	//construct b
	VectorD b(varMap.size());
	VectorD sol[3];
	for(int ic=0;ic<3;ic++){		
		int i=0;
		for(int x=0;x<maskImage.width();x++){
			for(int y=0;y<maskImage.height();y++){


				if(maskImage(x,y,ic)>0.9f){
					//graident from source image
					// float pixel=4*sourceImage(x,y,ic)-sourceImage(x-1,y,ic)-sourceImage(x+1,y,ic)
					// 			-sourceImage(x,y-1,ic)-sourceImage(x,y+1,ic);
					float v1=sourceImage(x,y,ic)-sourceImage(x-1,y,ic);
					float v2=sourceImage(x,y,ic)-sourceImage(x+1,y,ic);
					float v3=sourceImage(x,y,ic)-sourceImage(x,y-1,ic);
					float v4=sourceImage(x,y,ic)-sourceImage(x,y+1,ic);
					
					float pixel=clamp(pow(0.2/abs(v1),0.2),0.,0.4)*v1
								+clamp(pow(0.2/abs(v2),0.2),0.,0.4)*v2
								+clamp(pow(0.2/abs(v3),0.2),0.,0.4)*v3
								+clamp(pow(0.2/abs(v4),0.2),0.,0.4)*v4;
					//float pixel=v1+v2+v3+v4;



					//add boundary from target image
					if(maskImage(x+1,y,ic)<0.1f){
						pixel+=targetImage(x+1,y,ic);
					}
					if(maskImage(x-1,y,ic)<0.1f){
						pixel+=targetImage(x-1,y,ic);
					}
					if(maskImage(x,y+1,ic)<0.1f){
						pixel+=targetImage(x,y+1,ic);
					}
					if(maskImage(x,y-1,ic)<0.1f){
						pixel+=targetImage(x,y-1,ic);
					}
					b[i]=pixel;
					i++;
				}


			}
		}	
		//solve x and save each channel in sol
		sol[ic]=solver.solve(b);
	}
	//cout<<sol[0]<<endl;


	//write image to output
	FloatImage output(maskImage.width(),maskImage.height(),3);
	
	for (int ic = 0; ic < 3; ic++){
		for (int x = 0; x < maskImage.width(); x++){
			for (int y = 0; y < maskImage.height(); y++){
				if(maskImage(x,y,0)<0.1f){
					output(x,y,ic)=targetImage_copy(x,y,ic);
				}else{
					int j=y*maskImage.width()+x;
					float value=(float)pow(10,sol[ic][varMap[j]]);
					//float value=(float)sol[ic][varMap[j]];
					//output(x,y,ic)=clamp(value,0.f,1.f);
					output(x,y,ic)=value;
					//output(x,y,ic)=sourceImage(x,y,ic);
				}				
			}			
		}		
	}
	
	output.write(DATA_DIR "/output/gradient_lumi_man.png");
}




//poisson blending
void testBlend(){

	//load images	
	FloatImage sourceImage(DATA_DIR "/input/a6/source.png");
	FloatImage targetImage(DATA_DIR "/input/a6/target.png");
	FloatImage maskImage(DATA_DIR "/input/a6/mask.png");
	FloatImage targetImage_copy=targetImage;
	sourceImage=log10FloatImage(sourceImage);
	targetImage=log10FloatImage(targetImage);
	


	//check if three images have the same dimension


	//save every pixel inside the mask into a map
	//map(pixel,int number)
	map<unsigned int, unsigned int> varMap;
	{
		int number=0;
		for(int x=0;x<maskImage.width();x++){
			for(int y=0;y<maskImage.height();y++){
				//white part of the mask
				if(maskImage(x,y,0)>0.9f){
					varMap[y*maskImage.width()+x]=number;
					number++;
				}
			}
		}
	}


	//construct sparse matrix M in Mx = b
	vector<Triplet> mt;
	{			
		int i=0;		

		for(int x=0;x<maskImage.width();x++){
			for(int y=0;y<maskImage.height();y++){
				
				if(maskImage(x,y,0)>0.9f){
					int j=y*maskImage.width()+x;
					mt.push_back(Triplet(i,varMap[j],4));
				
					if(maskImage(x+1,y,0)>0.9f){
						mt.push_back(Triplet(i,varMap[j+1],-1));
					}
					if(maskImage(x-1,y,0)>0.9f){
						mt.push_back(Triplet(i,varMap[j-1],-1));
					}
					if(maskImage(x,y+1,0)>0.9f){
						mt.push_back(Triplet(i,varMap[j+maskImage.width()],-1));
					}
					if(maskImage(x,y-1,0)>0.9f){
						mt.push_back(Triplet(i,varMap[j-maskImage.width()],-1));
					}
					i++;
				}

			}
		}	

	}
	
	

	Eigen::SimplicialCholesky<SpMatrix> solver;
	{
		SpMatrix mat(varMap.size(),varMap.size());
		mat.setFromTriplets(mt.begin(),mt.end());
		//cout<<mat<<endl;
		solver.compute(mat);		
	}
	


	//construct b
	VectorD b(varMap.size());
	VectorD sol[3];
	for(int ic=0;ic<3;ic++){		
		int i=0;
		for(int x=0;x<maskImage.width();x++){
			for(int y=0;y<maskImage.height();y++){


				if(maskImage(x,y,ic)>0.9f){
					//graident from source image
					float pixel=4*sourceImage(x,y,ic)-sourceImage(x-1,y,ic)-sourceImage(x+1,y,ic)
								-sourceImage(x,y-1,ic)-sourceImage(x,y+1,ic);

					//add boundary from target image
					if(maskImage(x+1,y,ic)<0.1f){
						pixel+=targetImage(x+1,y,ic);
					}
					if(maskImage(x-1,y,ic)<0.1f){
						pixel+=targetImage(x-1,y,ic);
					}
					if(maskImage(x,y+1,ic)<0.1f){
						pixel+=targetImage(x,y+1,ic);
					}
					if(maskImage(x,y-1,ic)<0.1f){
						pixel+=targetImage(x,y-1,ic);
					}
					b[i]=pixel;
					i++;
				}


			}
		}	
		//solve x and save each channel in sol
		sol[ic]=solver.solve(b);
	}


	//write image to output
	FloatImage output(maskImage.width(),maskImage.height(),3);
	
	for (int ic = 0; ic < 3; ic++){
		for (int x = 0; x < maskImage.width(); x++){
			for (int y = 0; y < maskImage.height(); y++){
				if(maskImage(x,y,0)<0.1f){
					output(x,y,ic)=targetImage_copy(x,y,ic);
				}else{
					int j=y*maskImage.width()+x;
					float value=(float)pow(10,sol[ic][varMap[j]]);
					//float value=(float)sol[ic][varMap[j]];
					//output(x,y,ic)=clamp(value,0.f,1.f);
					output(x,y,ic)=value;
					//output(x,y,ic)=sourceImage(x,y,ic);
				}				
			}			
		}		
	}
	
	output.write(DATA_DIR "/output/gradient_copy_p.png");
}



// This is a way for you to test your functions.
// We will not grade the contents of the main function
int main()
{
	// uncomment these test functions as you complete the assignment to test your code

    // try { testSmartAccessor();}   catch(...) {cout << "testSmartAccessor Failed!" << endl;}
    //try { testScaling(); }        catch(...) {cout << "testScaling Failed!" << endl;}
    //try { testRotation(); }       catch(...) {cout << "testRotation Failed!" << endl;}
    //try { testWarpBy1(); }        catch(...) {cout << "testWarpBy1 Failed!" << endl;}

	// // TODO: replace '1' with your own student number
	//try { classMorph(34); }        catch(...) {cout << "classMorph Failed!" << endl;}
    //try { testMorph(); }          catch(...) {cout << "testMorph Failed!" << endl;}
	//try { testBlend(); }          catch(...) {cout << "testBlend Failed!" << endl;}
	//try { testMixedBlend(); }          catch(...) {cout << "testmixedBlend Failed!" << endl;}
	//try { testFlatten(); }          catch(...) {cout << "testmixedBlend Failed!" << endl;}
	try { testColor(); }          catch(...) {cout << "testColor Failed!" << endl;}
	//try { testLumi(); }          catch(...) {cout << "testLumi Failed!" << endl;}


	


}
