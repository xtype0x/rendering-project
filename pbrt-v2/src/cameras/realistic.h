#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_CAMERAS_REALISTIC_H
#define PBRT_CAMERAS_REALISTIC_H

#include "pbrt.h"
#include "camera.h"
#include "film.h"

// Example representation of an autofocus zone.
class AfZone {
	public:
	  // from the data file
	  float left, right;
	  float top, bottom;
	  int xres,yres;
};
//	add lens class
//	Tao Du
//	taodu@stanford.edu
//	May 17, 2014
class Lens
{
	public:
	float radius;
	float thickness;
	float refraction;
	float aperture;
	float zPos;	//	the position in the camera frame
};

class ThickLens
{
public:
	float ImgPos(float z);
	float F, P;
	float Fprime, Pprime;
	float f, fprime;
};

class RealisticCamera : public Camera {
public:
   RealisticCamera(const AnimatedTransform &cam2world,
      float hither, float yon, float sopen,
      float sclose, float filmdistance, float aperture_diameter,
      const string &specfile,
	  const string &autofocusfile,
      float filmdiag,
	  Film *film);
   ~RealisticCamera();
   float GenerateRay(const CameraSample &sample, Ray *) const;
   void  AutoFocus(Renderer * renderer, const Scene * scene, Sample * origSample);
   void  ParseAfZones(const string& filename);
	//	parse data file
	void ParseLens(const string& filename);
	//	thick lens approximation
	void CompThickLens();

	//	estimate the filmPos in autofocus
	Point EstimateAutoFocusPos(AfZone &zone, Renderer * renderer, const Scene * scene, Sample * origSample);

	//	evaluate the auto focus results
	float EvalAutoFocus(Renderer * renderer, const Scene * scene, Sample * origSample, float fPos, AfZone &zone);
private:
   bool  autofocus;
   vector<AfZone> afZones;
   float ShutterOpen;
   float ShutterClose;
   Film * film;
	//	store the size of the file
	float filmDiag;
	float filmLengthInX;
	float filmLengthInY;
	//	film position in the camera space
	float filmDist;
	//	absolute film distance in the camera space
	float filmPos;
	//	real aperture
	float aperture;

	//	add a vector of lens
	vector<Lens> lenses;

	//	thick lens approximation
	ThickLens tLens;

};

//	n1 represents the refraction in the side of Rin
//	n2 is the refraction on the other side
bool RefractFromLens(Lens lens, Ray Rin, Ray &Rout, float n1, float n2);

RealisticCamera *CreateRealisticCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film);

//	f measure we use
float varMeasurement(float *rgb, int height, int width);
float smlMeasurement(float *rgb, int height, int width);

#endif
