// cameras/realistic.cpp*
#include "stdafx.h"
#include "cameras/realistic.h"
#include "paramset.h"
#include "sampler.h"
#include "montecarlo.h"
#include "filters/box.h"
#include "film/image.h"
#include "samplers/stratified.h"
#include "intersection.h"
#include "renderer.h"
#include "reflection.h"

#include <fstream>
#include <algorithm>
#include <vector>
#include <string>
#include <iostream>

using namespace std;

//	Tao Du
//	thick lens approximation
float ThickLens::ImgPos(float z)
{
	float zprime = 1 / (1 / fprime + 1 / (z - P));
	return Pprime + zprime;
}

RealisticCamera *CreateRealisticCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film) {
	   // Extract common camera parameters from \use{ParamSet}
	   float hither = params.FindOneFloat("hither", -1);
	   float yon = params.FindOneFloat("yon", -1);
	   float shutteropen = params.FindOneFloat("shutteropen", -1);
	   float shutterclose = params.FindOneFloat("shutterclose", -1);

	   // Realistic camera-specific parameters
	   //string specfile = params.FindOneString("specfile", "");
	   string specfile = params.FindOneFilename("specfile", "");
	   float filmdistance = params.FindOneFloat("filmdistance", 70.0); // about 70 mm default to film
	   float fstop = params.FindOneFloat("aperture_diameter", 1.0);
	   float filmdiag = params.FindOneFloat("filmdiag", 35.0);
	   string autofocusfile = params.FindOneFilename("af_zones", "");
	   //string autofocusfile = params.FindOneString("af_zones", "");
	   assert(hither != -1 && yon != -1 && shutteropen != -1 &&
	      shutterclose != -1 && filmdistance!= -1);
	   if (specfile == "") {
	       Severe( "No lens spec file supplied!\n" );
	   }
	   return new RealisticCamera(cam2world, hither, yon,
	      shutteropen, shutterclose, filmdistance, fstop,
	      specfile, autofocusfile, filmdiag, film);
}

RealisticCamera::RealisticCamera(const AnimatedTransform &cam2world,
                                 float hither, float yon,
                                 float sopen, float sclose,
                                 float filmdistance, float aperture_diameter,
                                 const string &specfile,
								 const string &autofocusfile,
                                 float filmdiag,
								 Film *f)
                                 : Camera(cam2world, sopen, sclose, f),
								   ShutterOpen(sopen),
								   ShutterClose(sclose),
								   film(f)
{

	// YOUR CODE HERE -- build and store datastructures representing the given lens
	// and film placement.

	//	compute the film size
	filmDiag = filmdiag;
	float xRes = float(f->xResolution);
	float yRes = float(f->yResolution);
	float diagRes = sqrt(xRes * xRes + yRes * yRes);
	filmLengthInX = filmdiag / diagRes * xRes;
	filmLengthInY = filmdiag / diagRes * yRes;
	//	store the film position
	filmDist = filmdistance;
	//	store the aperture
	aperture = aperture_diameter;
	//	parse the specfile
	ParseLens(specfile);

	//	Tao Du
	//	thick lens approximation
	CompThickLens();

	// If 'autofocusfile' is the empty string, then you should do
	// nothing in any subsequent call to AutoFocus()
	autofocus = false;

	if (autofocusfile.compare("") != 0)  {
		std::cout << "autoocusfile is not empty!" << std::endl;
		ParseAfZones(autofocusfile);
		autofocus = true;
	}
}

void RealisticCamera::CompThickLens()
{
	//	fill out the value in tLens: F P Fprime Pprime f fprime
	//	trace the ray to compute P and F
		//	Tao Du
	//	thick lens approximation
	//	we first compute 
	int nLens = (int)lenses.size();
	Point Pstart(lenses[nLens - 1].aperture / 4, 0.f, filmPos);
	//	check for Pstart
	Ray Rin(Pstart, Vector(0.f, 0.f, 1.f), 0.f, INFINITY);
	//	start to iterate all the lenses
	Ray Rout;
	float n2;
	for (int i = nLens - 1; i >= 0; i--)
	{
		if (i > 0)
			n2 = lenses[i - 1].refraction;
		else
			n2 = 1.0f;
		bool succeed = RefractFromLens(lenses[i], Rin, Rout, lenses[i].refraction, n2);
		//	if Rin and lens won't intersect
		//	the function will return false		
		if (!succeed)
		{
			//	we don't want to see that actually
			//	otherwise we have to try another ray parallel
			//	and closer to z-axis...
			break;
		}		
		//	update Rin
		Rin = Rout;
	}
	//	compute the intersection point in z axis
	float t = -Rin.o.x / Rin.d.x;
	Point p = Rin(t);
	tLens.F = p.z;
	//	now, compute tLens.P
	t = (Pstart.x - Rin.o.x) / Rin.d.x;
	tLens.P = Rin.o.z + t * Rin.d.z;
	//	compute f
	tLens.f = tLens.F - tLens.P;

	//	now, shoot another ray from the object space
	Rin.o.z = 5.f;
	Rin.d = Vector(0.f, 0.f, -1.f);
	Ray Rstart = Rin;
	float n1;
	for (int i = 0; i < nLens; i++)
	{
		if (i == 0)
			n1 = 1.f;
		else
			n1 = lenses[i - 1].refraction;
		bool succeed = RefractFromLens(lenses[i], Rin, Rout, n1, lenses[i].refraction);
		//	if Rin and lens won't intersect
		//	the function will return false		
		if (!succeed)
		{
			//	we don't want to see that actually
			//	otherwise we have to try another ray parallel
			//	and closer to z-axis...
			break;
		}
		//	update Rin
		Rin = Rout;
	}
	//	compute the intersection point in z axis
	t = -Rin.o.x / Rin.d.x;
	p = Rin(t);
	tLens.Fprime = p.z;
	//	now, compute tLens.Pprime
	t = (Rstart.o.x - Rin.o.x) / Rin.d.x;
	tLens.Pprime = Rin.o.z + t * Rin.d.z;
	//	compute fprime
	tLens.fprime = tLens.Fprime - tLens.Pprime;
	//	print the results to double check
}

void RealisticCamera::ParseLens(const string& filename)
{
	ifstream specfile(filename.c_str());
	if (!specfile)
	{
		fprintf(stderr, "Cannot open file %s\n", filename.c_str());
		exit(-1);
	}
	char line[512];
	filmPos = 0.0;
	while (!specfile.eof()) 
	{
		specfile.getline(line, 512);
		if (line[0] != '\0' && line[0] != '#' && 
		line[0] != ' ' && line[0] != '\t' && line[0] != '\n')
		{
			lenses.resize(lenses.size() + 1);
			Lens& lens = lenses[lenses.size() - 1];
			sscanf(line, "%f %f %f %f\n", &lens.radius, &lens.thickness, &lens.refraction, &lens.aperture);
			//	if we find the aperture stop
			if (lens.radius == 0.0f)
			{
				lens.aperture = min(lens.aperture, aperture);
				//	aperture stop should be in air
				lens.refraction = 1.f;			
			}
			lens.zPos = filmPos;
			filmPos -= lens.thickness;
		}
	}
	filmPos -= filmDist;
	//printf("Read in %zu lens from %s\n", lenses.size(), filename.c_str());
	specfile.close();
}

// parses the AF zone file
void RealisticCamera::ParseAfZones(const string& filename)
{
  ifstream specfile(filename.c_str());
   if (!specfile) {
      fprintf(stderr, "Cannot open file %s\n", filename.c_str());
      exit (-1);
   }

   char line[512];

   while (!specfile.eof()) {
      specfile.getline(line, 512);
      if (line[0] != '\0' && line[0] != '#' &&
         line[0] != ' ' && line[0] != '\t' && line[0] != '\n')
      {
		afZones.resize(afZones.size()+1);
		AfZone& zone = afZones[afZones.size()-1];
		sscanf(line, "%f %f %f %f\n", &zone.left, &zone.right, &zone.top, &zone.bottom);
      }
   }
   specfile.close();
	//printf("Read in %zu AF zones from %s\n", afZones.size(), filename.c_str());
}

RealisticCamera::~RealisticCamera()
{
	
}

bool RefractFromLens(Lens lens, Ray Rin, Ray &Rout, float n1, float n2)
{
	//	compute the intersection of lens and Rin
	//	if they intersect
	//	return true and build a refracted ray in Rout
	//	otherwise return false
	//	compute the center of the circle

	//	if we encounter the aperture stop
	if (lens.radius == 0.0f)
	{
		//	the aperture shot
		//	Rin.o + t * Rin.d
		float zPos = lens.zPos;
		float t = (zPos - Rin.o.z) / Rin.d.z;
		Point Pinter = Rin(t);
		//	decide whether Pinter is inside the aperture
		float aper = lens.aperture;
		if (Pinter.x * Pinter.x + Pinter.y * Pinter.y > aper * aper / 4)
			return false;
		else
		{
			Rout = Rin;
			return true;
		}
	}
	//	else
	float lensAper = lens.aperture;
	float radius = lens.radius;
	float radiusSq = radius * radius;
	float zCenter = lens.zPos - radius;
	//	if zCenter and ray.o.z are in the same side of lens.zPos
	//	we should move the ray closer in case it will intersect with
	//	the other half
	//	compute the distance between center and the ray
	Point center(0.0, 0.0, zCenter);
	//	(Rin.o + t * Rin.d - center) * Rin.d = 0
	//	(Rin.o - center) * Rin.d + t * Rin.d * Rin.d = 0;
	float dSq = Rin.d.LengthSquared();
	float t = Dot(Rin.o - center, Rin.d) / -dSq;
	Point midPoint = Rin(t);
	float midSq = (midPoint - center).LengthSquared();
	if (midSq > radiusSq)	//	no intersection
		return false;
	//	find the intersection point
	float dLen = Rin.d.Length();
	float deltaT = sqrt((radiusSq - midSq)) / dLen;
	float t0 = t - deltaT;
	float t1 = t + deltaT;
	float tInt;
	//	three cases:
	//	t1 > t0 > 0.f
	//	t1 > 0.f > t0
	//	0.f > t1 > t0
	if (t1 <= 0.f)
		return false;
	else if (t0 < 0.f)
		tInt = t1;
	else	//	t0 > 0.f && t1 > 0.f
	{
		//	find the correct intersection point
		//	test whether t0 is a valid intersection point
		Point p = Rin(t0);
		if ((p.z - lens.zPos) * (p.z - zCenter) < 0.f
			&& p.x * p.x + p.y * p.y < lensAper * lensAper / 4)
			tInt = t0;
		else
		{
			//	t0 is invalid, test t1
			p = Rin(t1);
			if ((p.z -lens.zPos) * (p.z - zCenter) < 0.f)
				tInt = t1;
			else
				return false;
		}
	}		
	Point Pinter = Rin(tInt);
	float apSq = Pinter.x * Pinter.x + Pinter.y * Pinter.y;
	if (apSq > lensAper * lensAper / 4)
		return false;		//	no intersection within the aperture
	Rout.o = Pinter;
	//	now decide the direction of the refracted light
	Vector normal = Normalize(Pinter - center);
	if (Rin.d.z * normal.z < 0.f)
		normal = -normal;
	float cos1 = Dot(normal, Rin.d) / dLen;
	float sin1 = sqrt(1 - cos1 * cos1);
	float sin2 = n1 * sin1 / n2;
	if (sin2 > 1.0f)
		return false;
	//	now let's build the dir vector
	float cos2 = sqrt(1 - sin2 * sin2);
	Vector v1 = Cross(Rin.d, normal);
	//	in some weird case, Rin.d // normal
	if (v1.Length() == 0.f)
	{
		Rout = Rin;
		return true;
	}
	Vector v2 = Normalize(Cross(normal, v1));
	Rout.d = cos2 * normal + sin2 * v2;
	return true;
}

float RealisticCamera::GenerateRay(const CameraSample &sample, Ray *ray) const
{
  // YOUR CODE HERE -- make that ray!

  // use sample->imageX and sample->imageY to get raster-space coordinates
  // of the sample point on the film.
  // use sample->lensU and sample->lensV to get a sample position on the lens
	
	//	the range of imageX: 0 to xResolution
	//	the range of imageY: 0 to yResolution
	float xCamera, yCamera;
	xCamera = -(sample.imageX * 1.0 / film->xResolution - 0.5) * filmLengthInX;
	yCamera = (sample.imageY * 1.0 / film->yResolution - 0.5) * filmLengthInY;
	Point Pcamera(xCamera, yCamera, filmPos);

	//	generate lens samples
	float lensU, lensV;
	ConcentricSampleDisk(sample.lensU, sample.lensV, &lensU, &lensV);
	//	the last lens
	int nLens = (int)lenses.size();
	float LensRadius = lenses[nLens - 1].aperture / 2;
	lensU *= LensRadius;
	lensV *= LensRadius;
	float r = lenses[nLens - 1].radius;
	float delta = fabs(r) - sqrt(r * r - LensRadius * LensRadius);
	float z = lenses[nLens - 1].zPos;	
	if (r > 0.f)
	{
		z -= delta;
	}
	else
	{	
		z += delta;
	}
	Point Phit(lensU, lensV, z);
	Ray Rin(Pcamera, Phit - Pcamera, 0.f, INFINITY);
	//	start to iterate all the lenses
	Ray Rout;
	float n2;
	for (int i = nLens - 1; i >= 0; i--)
	{
		if (i > 0)
			n2 = lenses[i - 1].refraction;
		else
			n2 = 1.0f;
		bool succeed = RefractFromLens(lenses[i], Rin, Rout, lenses[i].refraction, n2);
		//	if Rin and lens won't intersect
		//	the function will return false		
		if (!succeed)
			return 0.0f;
		//	update Rin
		Rin = Rout;
	}
	//	return Rin
	*ray = Ray(Rin.o, Normalize(Vector(Rin.d)), 0.f, INFINITY);
	ray->time = sample.time;
	CameraToWorld(*ray, ray);
	ray->d = Normalize(ray->d);
	//	compute the weight
	Vector XX = Phit - Pcamera;
	float XXlenSq = XX.LengthSquared();
	float w = XX.z * LensRadius / XXlenSq;
	return M_PI * w * w;
}

Point RealisticCamera::EstimateAutoFocusPos(AfZone &zone, Renderer * renderer, const Scene * scene, Sample * origSample)
{
	RNG rng;
	MemoryArena arena;
	Filter * filter = new BoxFilter(.5f,.5f);
	//	consider only the central part of the zone
	float xcenter = (zone.left + zone.right) / 2;
	float ycenter = (zone.top + zone.bottom) / 2;
	float xsize = zone.left - xcenter;
	float ysize = zone.top - ycenter;
	float scale = 1.f;
	const float crop[] = {xcenter + scale * xsize, 
		xcenter - scale * xsize, 
		ycenter + scale * ysize,
		ycenter - scale * ysize};

	ImageFilm sensor(film->xResolution, film->yResolution, filter, crop, "foo.exr", false);
	int xstart,xend,ystart,yend;
	sensor.GetSampleExtent(&xstart,&xend,&ystart,&yend);

	StratifiedSampler sampler(xstart, xend, ystart, yend,
	                          16, 16, true, ShutterOpen, ShutterClose);

	// Allocate space for samples and intersections
	int maxSamples = sampler.MaximumSampleCount();
	Sample *samples = origSample->Duplicate(maxSamples);
	RayDifferential *rays = new RayDifferential[maxSamples];
	Spectrum *Ls = new Spectrum[maxSamples];
	Spectrum *Ts = new Spectrum[maxSamples];
	Intersection *isects = new Intersection[maxSamples];

	// Get samples from _Sampler_ and update image
	int sampleCount;
	Point pVal = Point(0.f, 0.f, 0.f);
	int pCnt = 0;
	Transform transform;
	while ((sampleCount = sampler.GetMoreSamples(samples, rng)) > 0) {
		// Generate camera rays and compute radiance along rays
		for (int i = 0; i < sampleCount; ++i) {
			// Find camera ray for _sample[i]_

			float rayWeight = this->GenerateRayDifferential(samples[i], &rays[i]);
			rays[i].ScaleDifferentials(1.f / sqrtf(sampler.samplesPerPixel));


			// Evaluate radiance along camera ray

			if (rayWeight > 0.f)
			{
				Ls[i] = rayWeight * renderer->Li(scene, rays[i], &samples[i], rng,
												 arena, &isects[i], &Ts[i]);
				//	extract the z value
				//	if we don't have intersection
				//	skip this sample
				if (!isects[i].primitive)
					continue;
				BSDF *bsdf = isects[i].GetBSDF(rays[i], arena);
    			const Point &p = bsdf->dgShading.p;		
				CameraToWorld.Interpolate(samples[i].time, &transform);
				transform = Inverse(transform);
				Point Pcamera = transform(p);
				pVal += Pcamera;
				pCnt++;			
			}			
			else {
				Ls[i] = 0.f;
				Ts[i] = 1.f;
			}

			// Issue warning if unexpected radiance value returned
			if (Ls[i].HasNaNs()) {
				Error("Not-a-number radiance value returned "
					  "for image sample.  Setting to black.");
				Ls[i] = Spectrum(0.f);
			}
			else if (Ls[i].y() < -1e-5) {
				Error("Negative luminance value, %f, returned"
					  "for image sample.  Setting to black.", Ls[i].y());
				Ls[i] = Spectrum(0.f);
			}
			else if (isinf(Ls[i].y())) {
				Error("Infinite luminance value returned"
					  "for image sample.  Setting to black.");
				Ls[i] = Spectrum(0.f);
			}

		}

		// Report sample results to _Sampler_, add contributions to image
		if (sampler.ReportResults(samples, rays, Ls, isects, sampleCount))
		{
			for (int i = 0; i < sampleCount; ++i)
			{
				sensor.AddSample(samples[i], Ls[i]);
			}
		}

		// Free _MemoryArena_ memory from computing image sample values
		arena.FreeAll();
	}
	//	report the average zVal
	pVal /= pCnt;
	return pVal;
}

void  RealisticCamera::AutoFocus(Renderer * renderer, const Scene * scene, Sample * origSample) {
	// YOUR CODE HERE:
	// The current code shows how to create a new Sampler, and Film cropped to the size of the auto focus zone.
	// It then renders the film, producing rgb values.  You need to:
	//
	// 1. Modify this code so that it can adjust film plane of the camera
	// 2. Use the results of raytracing to evaluate whether the image is in focus
	// 3. Search over the space of film planes to find the best-focused plane.

	//	now we have a ray from Pcamera, and it points towards
	if(!autofocus)
		return;
	std::cout << "autofocus is true !" << std::endl;
	int nZones = (int)afZones.size();
	Point *estPos = new Point[nZones];
	for (int i = 0; i < nZones; i++)
	{
		estPos[i] = EstimateAutoFocusPos(afZones[i], renderer, scene, origSample);
	}
	//	select the proper filmPos here
	int afOption;
	if (nZones == 1)
		afOption = 1;
	else
	{
		printf("select the auto focus mode:\n");
		printf("1: closest object\n");
		printf("2: object in the middle of the image\n");
		printf("3: farthest object\n");
		int retVal = scanf("%d", &afOption);
		if (retVal == EOF)
		{	
			//	do nothing
			printf("fail to input!\n");
			return;
		}
		//	default: 1: closest object
		if (afOption < 1 || afOption > 3)
			afOption = 1;
		
	}
	float objPos = 0.f;
	int zoneId = -1;
	switch (afOption)
	{
	case 2:
		{
			//	scan all the afZones, find the one closest to (0.5, 0.5)
			float minDist = 0.f;
			for (int i = 0; i < nZones; i++)
			{
				AfZone zone = afZones[i];
				float xmid = (zone.left + zone.right) / 2;
				float ymid = (zone.top + zone.bottom) / 2;
				float dist = (xmid - 0.5f) * (xmid - 0.5f) 
					+ (ymid - 0.5f) * (ymid - 0.5f);
				if (i == 0 || dist < minDist)
				{
					minDist = dist;
					zoneId = i;
				}
			}
			objPos = estPos[zoneId].z;
			AfZone zone = afZones[zoneId];
			printf("select zones: %f %f %f %f\n", 
				zone.left, zone.right, zone.top, zone.bottom);
		}
		break;
	case 3:
		for (int i = 0; i < nZones; i++)
		{
			if (i == 0 || objPos < estPos[i].z)
			{
				objPos = estPos[i].z;
				zoneId = i;		
			}	
		}
		break;	
	default:	//	case 1
		for (int i = 0; i < nZones; i++)
		{
			if (i == 0 || objPos > estPos[i].z)
			{
				objPos = estPos[i].z;
				zoneId = i;		
			}	
		}
	}
	filmPos = tLens.ImgPos(objPos);
	delete []estPos;
	//	now, do golden search in a small range of filmPos
	float tau = (sqrt(5) - 1) / 2;
	float a = filmPos - 5.f;
	float b = filmPos + 5.f;
	float x0 = a + (1 - tau) * (b - a);
	float x1 = a + tau * (b - a);
	AfZone zone = afZones[zoneId];
	float f0 = -EvalAutoFocus(renderer, scene, origSample, x0, zone);
	float f1 = -EvalAutoFocus(renderer, scene, origSample, x1, zone);
	while (b - a >= 0.5)
	{
		if (f0 >= f1)
		{
			a = x0;
			x0 = x1;
			f0 = f1;
			x1 = a + tau * (b - a);
			f1 = -EvalAutoFocus(renderer, scene, origSample, x1, zone);
		}
		else
		{
			b = x1;
			x1 = x0;
			f1 = f0;
			x0 = a + (1 - tau) * (b - a);
			f0 = -EvalAutoFocus(renderer, scene, origSample, x0, zone);
		}
	}
	filmPos = x1;
	printf("filmPos = %f\n", filmPos);
}

float RealisticCamera::EvalAutoFocus(Renderer * renderer, const Scene * scene, Sample * origSample, float fPos, AfZone &zone)
{
	filmPos = fPos;
	RNG rng;
	MemoryArena arena;
	Filter * filter = new BoxFilter(.5f,.5f);
	const float crop[] = {zone.left,zone.right,zone.top,zone.bottom};
	ImageFilm sensor(film->xResolution, film->yResolution, filter, crop,"foo.exr",false);
	int xstart,xend,ystart,yend;
	sensor.GetSampleExtent(&xstart,&xend,&ystart,&yend);
	StratifiedSampler sampler(xstart, xend, ystart, yend,
	                          16, 16, true, ShutterOpen, ShutterClose);
	// Allocate space for samples and intersections
	int maxSamples = sampler.MaximumSampleCount();
	Sample *samples = origSample->Duplicate(maxSamples);
	RayDifferential *rays = new RayDifferential[maxSamples];
	Spectrum *Ls = new Spectrum[maxSamples];
	Spectrum *Ts = new Spectrum[maxSamples];
	Intersection *isects = new Intersection[maxSamples];
	// Get samples from _Sampler_ and update image
	int sampleCount;
	while ((sampleCount = sampler.GetMoreSamples(samples, rng)) > 0) {
		// Generate camera rays and compute radiance along rays
		for (int i = 0; i < sampleCount; ++i) {
			// Find camera ray for _sample[i]_

			float rayWeight = this->GenerateRayDifferential(samples[i], &rays[i]);
			rays[i].ScaleDifferentials(1.f / sqrtf(sampler.samplesPerPixel));


			// Evaluate radiance along camera ray

			if (rayWeight > 0.f)
				Ls[i] = rayWeight * renderer->Li(scene, rays[i], &samples[i], rng,
												 arena, &isects[i], &Ts[i]);
			else {
				Ls[i] = 0.f;
				Ts[i] = 1.f;
			}

			// Issue warning if unexpected radiance value returned
			if (Ls[i].HasNaNs()) {
				Error("Not-a-number radiance value returned "
					  "for image sample.  Setting to black.");
				Ls[i] = Spectrum(0.f);
			}
			else if (Ls[i].y() < -1e-5) {
				Error("Negative luminance value, %f, returned"
					  "for image sample.  Setting to black.", Ls[i].y());
				Ls[i] = Spectrum(0.f);
			}
			else if (isinf(Ls[i].y())) {
				Error("Infinite luminance value returned"
					  "for image sample.  Setting to black.");
				Ls[i] = Spectrum(0.f);
			}

		}

		// Report sample results to _Sampler_, add contributions to image
		if (sampler.ReportResults(samples, rays, Ls, isects, sampleCount))
		{
			for (int i = 0; i < sampleCount; ++i)
			{

				sensor.AddSample(samples[i], Ls[i]);

			}
		}

		// Free _MemoryArena_ memory from computing image sample values
		arena.FreeAll();
	}

	float * rgb;
	int width;
	int height;
	sensor.WriteRGB(&rgb,&width,&height,1.f);
	// YOUR CODE HERE! The rbg contents of the image for this zone
	// are now stored in the array 'rgb'.  You can now do whatever
	// processing you wish
	
	//	first, scale the color to [0, 1]
	float maxColor = 0.0;
	for (int c = 0; c < width * height * 3; c++)
		if (rgb[c] > maxColor)
			maxColor = rgb[c];
	if (maxColor == 0.f)
		return -1.f;
	for (int c = 0; c < width * height * 3; c++)
		rgb[c] /= maxColor;
	
	//	the layout of rgb:
	//	width, width, width, ... width
	//	rgbrgbrgbrgb ... rgbrgbrgb
	//	compute the focus measure
	//float fMeasure = varMeasurement(rgb, height, width);
	float fMeasure = smlMeasurement(rgb, height, width);
	printf("filmPos = %f fMeasure = %f\n", filmPos, fMeasure);
	//you own rgb  now so make sure to delete it:
	delete [] rgb;
	//if you want to see the output rendered from your sensor, uncomment this line (it will write a file called foo.exr)
	//sensor.WriteImage(1.f);
	delete[] samples;
	delete[] rays;
	delete[] Ls;
	delete[] Ts;
	delete[] isects;
	return fMeasure;
}

float varMeasurement(float *rgb, int height, int width)
{
	//	compute the variance in three channels
	float fMeasure = 0.f;
	float mean[3] = {0.f, 0.f, 0.f};
	float var[3] = {0.f, 0.f, 0.f};
	for (int c = 0; c < height * width; c++)
	{
		mean[0] += rgb[3 * c];
		mean[1] += rgb[3 * c + 1];
		mean[2] += rgb[3 * c + 2];
	}
	mean[0] /= (height * width);
	mean[1] /= (height * width);
	mean[2] /= (height * width);
	for (int c = 0; c < height * width; c++)
	{
		var[0] += ((rgb[3 * c] - mean[0]) * (rgb[3 * c] - mean[0]));	
		var[1] += ((rgb[3 * c + 1] - mean[1]) * (rgb[3 * c + 1] - mean[1]));	
		var[2] += ((rgb[3 * c + 2] - mean[2]) * (rgb[3 * c + 2] - mean[2]));	
	}
	var[0] /= (height * width);
	var[1] /= (height * width);
	var[2] /= (height * width);
	var[0] = sqrt(var[0]);
	var[1] = sqrt(var[1]);
	var[2] = sqrt(var[2]);
	fMeasure = var[0] + var[1] + var[2];
	return fMeasure;
}

float smlMeasurement(float *rgb, int height, int width)
{
	float fMeasure = 0.f;
	//	let's use a 3x3 window, i.e., step = 1
	int step = 1;
	for (int w = step; w < width - step; w++)
	{
		for (int h = step; h < height - step; h++)
		{
			//	extract r g and b
			int id = h * width + w;
			int left = id - step;
			int right = id + step;
			int up = id - width * step;
			int down = id + width * step;
			for (int channel = 0; channel < 3; channel++)
			{
				float color = rgb[3 * id + channel];
				fMeasure += (fabs(2 * color - rgb[3 * left + channel] - rgb[3 * right + channel]) 
					+ fabs(2 * color - rgb[3 * up + channel] - rgb[3 * down + channel]));
			}
		}	
	}
	return fMeasure;
}
