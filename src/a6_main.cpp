#include "floatimage.h"
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
	//return std::numeric_limits<float>::min();
	return 1.17549e-37;
}

// image --> log10FloatImage
FloatImage log10FloatImage(const FloatImage &im)
{
	// Taking a linear image im, transform to log10 scale.
	// To avoid infinity issues, make any 0-valued pixel be equal the the minimum
	// non-zero value. See image_minnonzero(im).
	FloatImage out=im;
	for(int i=0;i<im.size();i++){
		if(im(i)==0.f){
			out(i)=image_minnonzero(im);
		}else{
			out(i)=log10(im(i));
		}

		
	}
	return out; // change this
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

//poisson blending
void testMixedBlend(){

	//load images	
	FloatImage sourceImage(DATA_DIR "/input/a6/source_t.png");
	FloatImage targetImage(DATA_DIR "/input/a6/target_t.png");
	FloatImage maskImage(DATA_DIR "/input/a6/mask_t.png");
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
					output(x,y,ic)=clamp(value,0.f,1.f);
					//output(x,y,ic)=value;
					//output(x,y,ic)=sourceImage(x,y,ic);
				}				
			}			
		}		
	}
	
	output.write(DATA_DIR "/output/gradient_linear_t_mixed.png");
}


//poisson blending
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


//poisson blending
void testColor(){

	//load images	
	FloatImage targetImage(DATA_DIR "/input/a6/source_f.png");
	//FloatImage targetImage(DATA_DIR "/input/a6/target_f.png");
	FloatImage sourceImage=targetImage;
	for(int i=0;i<sourceImage.width();i++){
		for(int j=0;j<sourceImage.height();j++){
			sourceImage(i,j,0)*=1.5;
			sourceImage(i,j,1)*=0.5;
			sourceImage(i,j,2)*=0.5;
		}
	}
	
	FloatImage maskImage(DATA_DIR "/input/a6/mask_f.png");
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
	
	output.write(DATA_DIR "/output/change_color_f.png");
}



//poisson blending
void testBlend(){

	//load images	
	FloatImage sourceImage(DATA_DIR "/input/a6/source_p.png");
	FloatImage targetImage(DATA_DIR "/input/a6/target_p.png");
	FloatImage maskImage(DATA_DIR "/input/a6/mask_p.png");
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



int main()
{
	//testBlend();
	//testMixedBlend();
	//testFlatten();
	testColor();
}
