#include "floatimage.h"
#include <map>
#include "utils.h"
#include "Eigen/Sparse"
#include <iostream>
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

// util functions for gradient descent
// binary search on pos_map
int bs(const vector<pair<unsigned int, pair<vector<unsigned int>, vector<float>>>>& pos_map, int l, int r, int x) {
		if (r >= l) { 
        int mid = l + (r - l) / 2; 
        if (pos_map[mid].first == x) 
            return mid; 
        if (pos_map[mid].first > x) 
            return bs(pos_map, l, mid - 1, x); 
        return bs(pos_map, mid + 1, r, x); 
    } 
    return -1; 
} 
void set_vector(vector<VectorD>& v, const FloatImage& im, const vector<vector<pair<unsigned int, pair<vector<unsigned int>, vector<float>>>>>& pos_map) {
	for (int k = 0; k < im.sizeZ(); k++)
		for (int i = 0; i < pos_map[k].size(); i++)
			v[k](i) = im(pos_map[k][i].first % im.width(), pos_map[k][i].first / im.width(), k);
}
void fill_image(const vector<VectorD>& v, FloatImage& im, const vector<vector<pair<unsigned int, pair<vector<unsigned int>, vector<float>>>>>& pos_map) {
	for (int k = 0; k < im.sizeZ(); k++) 
		for (int i = 0; i < pos_map[k].size(); i++) im(pos_map[k][i].first % im.width(), pos_map[k][i].first / im.width(), k) = v[k](i);
}
VectorD conv(const VectorD &v, const vector<pair<unsigned int, pair<vector<unsigned int>, vector<float>>>>& pos_map) {
	VectorD output = v;
	float val;
	for (int i = 0; i < pos_map.size(); i++) {
		val = 4 * v(i);
		for (int j = 0; j < pos_map[i].second.first.size(); j++) val -= v(pos_map[i].second.first[j]);
		for (int j = 0; j < pos_map[i].second.second.size(); j++) val -= pos_map[i].second.second[j];
		output(i) = val;
	}
	return output;
}
 vector<VectorD> conv_3d(const vector<VectorD> &vs, vector<vector<pair<unsigned int, pair<vector<unsigned int>, vector<float>>>>>& pos_map) {
	vector<VectorD> output;
	for (int k = 0; k < vs.size(); k++) output.push_back(conv(vs[k], pos_map[k]));
	return output;
}
void fill_map(vector<vector<pair<unsigned int, pair<vector<unsigned int>, vector<float>>>>>& pos_map, const FloatImage& maskImage, const FloatImage& im) {
	int width = maskImage.width();
	for (int i = 0; i < pos_map[0].size(); i++) 
		for (int k = 0; k < pos_map.size(); k++) {
			pos_map[k][i].second = make_pair(vector<unsigned int>(), vector<float>());
			int mx = pos_map[k][i].first % width;
			int my = pos_map[k][i].first / width;
			if (maskImage(mx-1,my,k)<0.9f) pos_map[k][i].second.second.push_back(im(mx - 1,my,k));
			else pos_map[k][i].second.first.push_back(i - 1);
			if (maskImage(mx+1,my,k)<0.9f) pos_map[k][i].second.second.push_back(im(mx + 1,my,k));
			else pos_map[k][i].second.first.push_back(i + 1);
			if (maskImage(mx,my-1,k)<0.9f) pos_map[k][i].second.second.push_back(im(mx,my - 1,k));
			else pos_map[k][i].second.first.push_back(bs(pos_map[k], 0, pos_map[0].size() - 1, pos_map[k][i].first - width));
			if (maskImage(mx,my+1,k)<0.9f) pos_map[k][i].second.second.push_back(im(mx,my + 1,k));
			else pos_map[k][i].second.first.push_back(bs(pos_map[k], 0, pos_map[0].size() - 1, pos_map[k][i].first + width));
		}
}

//poisson blending: assuming boundary of the masked area does not overlap with boundary of the target image
void poisson_blending_gradient_descent() {
	//load images	
	FloatImage sourceImage(DATA_DIR "/input/a6/source_p.png");
	FloatImage targetImage(DATA_DIR "/input/a6/target_p.png");
	FloatImage maskImage(DATA_DIR "/input/a6/mask_p.png");
	sourceImage = log10FloatImage(sourceImage);
	targetImage = log10FloatImage(targetImage);
	// FloatImage sourceImage(DATA_DIR "/input/a6/src.png");
	// FloatImage targetImage(DATA_DIR "/input/a6/tar.png");
	// FloatImage maskImage(DATA_DIR "/input/a6/msk.png");
	FloatImage output(targetImage);

	int iter = 10000, start_from = 0;
	// real coord, internal coords, boundary values
	vector<vector<pair<unsigned int, pair<vector<unsigned int>, vector<float>>>>> pos_map;
	int marked_pixel, width = maskImage.width();
	vector<VectorD> b, t, r, tmp;

	// construct the position map for the target area
		// initialize pos_map
	for (int k = 0; k < targetImage.sizeZ(); k++) pos_map.push_back(vector<pair<unsigned int, pair<vector<unsigned int>, vector<float>>>>());
		// save all the masked coords at first slots
	for(int y=0;y<maskImage.height();y++)
		for(int x=0;x<maskImage.width();x++) {
			marked_pixel = y*maskImage.width()+x;
			if(maskImage(x,y,0)>0.9f)
				for (int k = 0; k < targetImage.sizeZ(); k++)
					pos_map[k].push_back(make_pair(marked_pixel, make_pair(vector<unsigned int>(), vector<float>())));
		}
	cout<< "pos_map constructed: "<< pos_map.size()<< " channels."<< endl;

	// initialize vectors
	for (int k = 0; k < targetImage.sizeZ(); k++) {
		b.push_back(VectorD());
		b[k].resize(pos_map[0].size());
		t.push_back(VectorD());
		t[k].resize(pos_map[0].size());
		r.push_back(VectorD());
		r[k].resize(pos_map[0].size());
		tmp.push_back(VectorD());
		tmp[k].resize(pos_map[0].size());
	}

	// calculate b
	set_vector(b, sourceImage, pos_map);
		// fill in pos_map for b
	fill_map(pos_map, maskImage, sourceImage);
	b = conv_3d(b, pos_map);
	// save a preview of b
	FloatImage sourceImage_conv(sourceImage);
	fill_image(b, sourceImage_conv, pos_map);
	exp10FloatImage(sourceImage_conv).write(DATA_DIR "/output/gradient_b.png");
	
	// initiate t_0: if start_from is specified, then load t_start_from
	if (start_from) {
		output = FloatImage(DATA_DIR "/output/log_gradient_descent_" + to_string(start_from) + ".png");
		set_vector(t, output, pos_map);
	} else {
		for (int k = 0; k < targetImage.sizeZ(); k++) 
			for (int i = 0; i < pos_map[k].size(); i++) t[k](i) = 0;	
	// save a preview of t_0
		fill_image(t, output, pos_map);
		exp10FloatImage(output).write(DATA_DIR "/output/log_gradient_descent_0.png");
	}
	// gradient descent
		// fill in pos_map for t
	fill_map(pos_map, maskImage, targetImage);
	for (size_t i = start_from + 1; i <= iter; i++) {
		// At_i -> tmp
		tmp = conv_3d(t, pos_map);
		// r_i
		for (int k = 0; k < output.sizeZ(); k++) r[k] = b[k] - tmp[k];
		// A^Tr_i -> tmp
		tmp = conv_3d(r, pos_map);
		// t_i -> t_i+1
		for (int k = 0; k < output.sizeZ(); k++) 
			t[k] = t[k] + (float)(r[k].transpose() * r[k])/(float)(r[k].transpose() * tmp[k]) * r[k];
		cout<< "iter: "<< i<< endl;
		if (!(i % 100)) {
			fill_image(t, output, pos_map);
			exp10FloatImage(output).write(DATA_DIR "/output/log_gradient_descent_" + to_string(i) + ".png");
		}
	}
}
// conjugate gradient
void poisson_blending_conjugate_gradient() {
	//load images	
	FloatImage sourceImage(DATA_DIR "/input/a6/source_p.png");
	FloatImage targetImage(DATA_DIR "/input/a6/target_p.png");
	FloatImage maskImage(DATA_DIR "/input/a6/mask_p.png");
	sourceImage = log10FloatImage(sourceImage);
	targetImage = log10FloatImage(targetImage);
	FloatImage output(targetImage);

	int iter = 6000, start_from = 0;
	// real coord, internal coords, boundary values
	vector<vector<pair<unsigned int, pair<vector<unsigned int>, vector<float>>>>> pos_map;
	int marked_pixel, width = maskImage.width();
	vector<VectorD> b, t, r, tmp, w, z;
	vector<float> a;

	// construct the position map for the target area
		// initialize pos_map
	for (int k = 0; k < targetImage.sizeZ(); k++) pos_map.push_back(vector<pair<unsigned int, pair<vector<unsigned int>, vector<float>>>>());
		// save all the masked coords at first slots
	for(int y=0;y<maskImage.height();y++)
		for(int x=0;x<maskImage.width();x++) {
			marked_pixel = y*maskImage.width()+x;
			if(maskImage(x,y,0)>0.9f)
				for (int k = 0; k < targetImage.sizeZ(); k++)
					pos_map[k].push_back(make_pair(marked_pixel, make_pair(vector<unsigned int>(), vector<float>())));
		}
	cout<< "pos_map constructed: "<< pos_map.size()<< " channels."<< endl;

	// initialize vectors
	for (int k = 0; k < targetImage.sizeZ(); k++) {
		b.push_back(VectorD());
		b[k].resize(pos_map[0].size());
		t.push_back(VectorD());
		t[k].resize(pos_map[0].size());
		r.push_back(VectorD());
		r[k].resize(pos_map[0].size());
		tmp.push_back(VectorD());
		tmp[k].resize(pos_map[0].size());
		w.push_back(VectorD());
		w[k].resize(pos_map[0].size());
		z.push_back(VectorD());
		z[k].resize(pos_map[0].size());
		a.push_back(0);
	}

	// calculate b
	set_vector(b, sourceImage, pos_map);
		// fill in pos_map for b
	fill_map(pos_map, maskImage, sourceImage);
	b = conv_3d(b, pos_map);
	
	// initial t: a random guess
	if (start_from) {
		output = FloatImage(DATA_DIR "/output/cg_gradient_descent_" + to_string(start_from) + ".png");
		set_vector(t, output, pos_map);
	} else {
		for (int k = 0; k < targetImage.sizeZ(); k++) 
			for (int i = 0; i < pos_map[k].size(); i++) t[k](i) = 0;	
	}

	// fill in pos_map for t
	fill_map(pos_map, maskImage, targetImage);
	// r0, w0
		// At_i -> tmp
	tmp = conv_3d(t, pos_map);
	for (int k = 0; k < output.sizeZ(); k++) {
		r[k] = b[k] - tmp[k];
		w[k] = -r[k];
	}
	// z
	z = conv_3d(w, pos_map);
	// a
	cout<< "getting a!"<< endl;
	for (int k = 0; k < output.sizeZ(); k++) a[k] = (float)(r[k].transpose() * w[k])/(float)(w[k].transpose() * z[k]);
	cout<< "got a!"<< endl;
	// t0
	for (int k = 0; k < output.sizeZ(); k++) t[k] = t[k] + a[k] * w[k];
	cout<< "entering cg"<< endl;
	// conjugate gradient
	for (size_t i = start_from + 1; i <= iter; i++) {
		unsigned int err_sum = 0;
		// r_i+1
		for (int k = 0; k < output.sizeZ(); k++) {
			r[k] = r[k] - a[k] * z[k];
			err_sum += r[k].norm();
		}
		// return condition
		if (err_sum < 1e-10) break;
		// w_i+1
		for (int k = 0; k < output.sizeZ(); k++)
			w[k] = -r[k] + (float)(r[k].transpose()*z[k])/(float)(w[k].transpose()*z[k])*w[k];
		// z_i+1
		z = conv_3d(w, pos_map);
		// a_i+1
		for (int k = 0; k < output.sizeZ(); k++) 
			a[k] = (float)(r[k].transpose() * w[k])/(float)(w[k].transpose() * z[k]);
		// t_i -> t_i+1
		for (int k = 0; k < output.sizeZ(); k++) 
			t[k] = t[k] + a[k]*w[k];
		cout<< "iter: "<< i<< endl;
		if (!(i % 100)) {
			fill_image(t, output, pos_map);
			exp10FloatImage(output).write(DATA_DIR "/output/cg_gradient_descent_" + to_string(i) + ".png");
		}
	}
}

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
	// try { testBlend(); }          catch(...) {cout << "testBlend Failed!" << endl;}
	// try { testMixedBlend(); }          catch(...) {cout << "testmixedBlend Failed!" << endl;}
	// try { testFlatten(); }          catch(...) {cout << "testmixedBlend Failed!" << endl;}
	// try { testColor(); }          catch(...) {cout << "testColor Failed!" << endl;}
	// try { testLumi(); }          catch(...) {cout << "testLumi Failed!" << endl;}
	try { poisson_blending_gradient_descent(); }	catch(...) {cout << "poisson_blending_gradient_descent Failed!" << endl;}

	


}
