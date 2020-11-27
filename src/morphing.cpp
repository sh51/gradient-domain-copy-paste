// morphing.cpp
// Assignment 6


#include "morphing.h"
#include "utils.h"
#include "a2.h"
#include "a6.h"

using namespace std;

/**************************************************************
 //            IMAGE WARPING/MORPHING FUNCTIONS              //
 *************************************************************/

// warp an entire image according to a pair of segments.
FloatImage warpBy1(const FloatImage &im, Segment &segBefore, Segment &segAfter)
{

	FloatImage imOut(im.width(), im.height(), im.channels());

	for (int y = 0; y < im.height(); y++)
	{
		for (int x = 0; x < im.width(); x++)
		{
			float floatX = (float) x;
			float floatY = (float) y;
			float u = segAfter.getU(floatX, floatY);
			float v = segAfter.getV(floatX, floatY);
			auto s = segBefore.UVtoX(u, v);

			for (int z = 0; z < im.channels(); z++)
				imOut(x, y, z) = interpolateLin(im, s[0], s[1], z, true);
		}
	}

	return imOut;

}

// warp an image according to a vector of before and after segments using segment weighting
FloatImage warp(const FloatImage &im, vector<Segment> &segsBefore, vector<Segment> &segsAfter, float a, float b,
                float p)
{

	float dist, u, v, w, sumW;
	float floatX, floatY;
	Vec2f newPoint;
	FloatImage imOut(im.width(), im.height(), im.channels());

	//todo: check segsBefore size is the same as segsAfter size
	for (auto y : range(im.height()))
	{
		for (auto x : range(im.width()))
		{
			Vec2f pointAfter;
			pointAfter.fill(0.f);
			sumW = 0;

			for (auto i : range(segsBefore.size()))
			{
				floatX = (float) x;
				floatY = (float) y;

				dist = segsAfter[i].dist(floatX, floatY);
				u = segsAfter[i].getU(floatX, floatY);
				v = segsAfter[i].getV(floatX, floatY);
				newPoint = segsBefore[i].UVtoX(u, v);

				//w = segsAfter[i].weight(newPoint[0], newPoint[1], a, b, p);
				//w = segsAfter[i].weight(floatX, floatY, a, b, p);

				w = segsBefore[i].weight(newPoint[0], newPoint[1], a, b, p);
				

				pointAfter[0] += newPoint[0] * w;
				pointAfter[1] += newPoint[1] * w;
				sumW += w;
			}

			pointAfter[0] /= sumW;
			pointAfter[1] /= sumW;

			for (auto z : range(im.channels()))
				imOut(x, y, z) = interpolateLin(im, pointAfter[0], pointAfter[1], z, true);
		}
	}
	return imOut;

}

// return a vector of N images in addition to the two inputs that morphs going between im1 and im2 for the corresponding segments
vector<FloatImage> morph(const FloatImage &im1, const FloatImage &im2, vector<Segment> &segsBefore,
                         vector<Segment> &segsAfter, int N, float a, float b, float p)
{

	float alpha, beta;
	vector<FloatImage> imMorph;
	int j;

	if ((im1.width() != im2.width()) ||
		(im1.height() != im2.height()) ||
		(im1.channels() != im2.channels()))
		throw MismatchedDimensionsException();

	imMorph.push_back(im1);

	N = N + 2;
	for (int i = 1; i < N - 1; i++)
	{
		alpha = ((float) i) / ((float) N - 1.0);
		beta = 1.0 - alpha;
		j = 0;

		vector<Segment> segsIntermediate;
		for (auto s : range(segsBefore.size()))
		{
			Segment sb = segsBefore[s];
			Segment sa = segsAfter[s];

			Segment si(alpha * sa.P[0] + beta * sb.P[0], alpha * sa.P[1] + beta * sb.P[1],
				       alpha * sa.Q[0] + beta * sb.Q[0], alpha * sa.Q[1] + beta * sb.Q[1]);

			segsIntermediate.push_back(si);

			j++;
		}

		FloatImage bwarp = warp(im1, segsBefore, segsIntermediate);
		FloatImage awarp = warp(im2, segsAfter, segsIntermediate);

		FloatImage bwarp2 = brightness(bwarp, beta);
		FloatImage awarp2 = brightness(awarp, alpha);
		FloatImage tmpMorph = awarp2 + bwarp2;

		imMorph.push_back(tmpMorph);

		segsIntermediate.clear();
	}

	imMorph.push_back(im2);
	return imMorph;

}

/**************************************************************
 //                 SEGMENT CLASS FUNCTIONS                  //
 *************************************************************/

// Segment constructor takes in 2 points (x1,y1) and (x2,y2) corresponding to the ends of a segment and computes:
// P - 2-element vector to point (x1, y1)
// Q - 2-element vector to pont (x2, y2)
// PQ - 2-element vector from P to Q
// PQ2 - squared distance between P and Q
// PQlength - distance between P and Q
// PQDivByPQlength2 - 2-element vector PQ normalized by PQ2
// perpPQ - 2-element vector perpendicular to PQ
// perpPQDivByPQlength - 2-element vector perpPQ normalized by PQlength
Segment::Segment(float x1, float y1, float x2, float y2)
{
	P.fill(0.f);
	Q.fill(0.f);
	perpPQ.fill(0.f);


	P = {{x1, y1}};
	Q = {{x2, y2}};

	// subtract vector P from Q to get the vector PQ
	PQ = subtract(P, Q);

	// get the dot product between PQ and itself and then take the square root to get the length of the vector PQ
	PQ2 = dot(PQ, PQ);
	PQlength = sqrt(PQ2);

	// normalize the PQ vector by the square of its distance (why square?)
	float normalizer1 = 1 / PQ2;
	PQDivByPQlength2 = scalarMult(PQ, normalizer1);

	//get the perpendicual vector and normalize it
	perpPQ[0] = -PQ[1];
	perpPQ[1] = PQ[0];
	float normalizer2 = 1 / PQlength;
	perpPQDivByPQlength = scalarMult(perpPQ, normalizer2);

}

// Implement the computation of the u coordinate of a point (x,y) with respect to a segment
float Segment::getU(float x, float y)
{

	return dot(subtract(P, Vec2f{{x, y}}), PQDivByPQlength2);

}

// Implement the computation of the v coordinate of a point (x,y) with respect to a segment
float Segment::getV(float x, float y)
{

	return dot(subtract(P, Vec2f{{x, y}}), perpPQDivByPQlength);

}

// compute the new (x, y) position of a point given by the (u,v) location relative to another segment.
// return the point (x,y) in a 2-element vector
Vec2f Segment::UVtoX(float u, float v)
{
	// takes the u,v values and returns a coordinate - to be called from target segment


	Vec2f PQmultU = scalarMult(PQ, u);
	Vec2f perpNormPQmultV = scalarMult(perpPQDivByPQlength, v);
	return add(P, add(PQmultU, perpNormPQmultV));

}

// Implement distance from a point (x,y) to the segment. Remember the 3 cases from class
float Segment::dist(float x, float y)
{

	float u = getU(x, y);
	float v = getV(x, y);
	float dist2;

	if (u < 0)
		dist2 = v * v + u * u * PQ2;
	else if (u > 1)
		dist2 = v * v + (u - 1) * (u - 1) * PQ2;
	else
		dist2 = v * v;

	return sqrt(dist2);

}

// compute the weight of a segment to a point (x,y) given the weight parameters a,b, and p
float Segment::weight(float x, float y, float a, float b, float p)
{

	return pow((pow(PQlength, p) / (a + dist(x, y))), b);

}

/**************************************************************
 //               DON'T EDIT BELOW THIS LINE                //
 *************************************************************/

// add 2 vectors of the same length.
Vec2f Segment::add(const Vec2f &vec1, const Vec2f &vec2)
{
	return Vec2f{{vec2[0]+vec1[0], vec2[1]+vec1[1]}};
}

// subtracts 2 vectors of the same length.
Vec2f Segment::subtract(const Vec2f &vec1, const Vec2f &vec2)
{
	return Vec2f{{vec2[0]-vec1[0], vec2[1]-vec1[1]}};
}

// takes the dot product between 2 vectors of the same length
float Segment::dot(const Vec2f &vec1, const Vec2f &vec2)
{
	return vec2[0]*vec1[0] + vec2[1]*vec1[1];
}

// mutliplies an entire vector by a scalor value
Vec2f Segment::scalarMult(const Vec2f &vec, float factor)
{
	return Vec2f{{vec[0]*factor, vec[1]*factor}};
}

// destructor
Segment::~Segment() {} // Nothing to clean up
