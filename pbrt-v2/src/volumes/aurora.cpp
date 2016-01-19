/*
Tao Du
taodu@stanford.edu
Jun 4, 2014
*/

// volumes/aurora.cpp*
#include "stdafx.h"
#include "volumes/aurora.h"
#include "paramset.h"
#include "core/perlin.h"
#include "core/perlin2d.h"
#include <iostream>
using namespace std;

AuroraGrid::AuroraGrid(const BBox &e, int x, int y, int z, float r)
    :noise2d(0.5,0.2,3,5)
{
	Vector vox = e.pMax - e.pMin;
	/*
	float dx = vox.x / x;
	float dy = vox.y / y;
	float dz = vox.z / z;
	step = max(dx, max(dy, dz));
	step = max(step, 2 * radius);
	*/
	step = 2 * r;
	nx = int(vox.x / step) + 1;
	ny = int(vox.y / step) + 1;
	nz = int(vox.z / step) + 1;
	vox.x = step * nx;
	vox.y = step * ny;
	vox.z = step * nz;
	extent.pMin = e.pMin;
	extent.pMax = e.pMin + vox;
	grids = new AuroraVoxel[nx * ny * nz];
	radius = r;
}


void AuroraGrid::AddPhoton(const AuroraPhoton &photon)
{
	Point p = photon.p;
	if (!extent.Inside(p))
		return;
	Vector vox = photon.p - extent.pMin;
    int vx = int(vox.x / step);
	int vy = int(vox.y / step);
	int vz = int(vox.z / step);
	if (vx == nx || vy == ny || vz == nz)
		return;
	grids[vz*nx*ny + vy*nx + vx].photons.push_back(photon);
}

void AuroraGrid::SearchInGrid(const Point &p, float &r, float &g, float &b) const
{
	r = g = b = 0.f;
	if (!extent.Inside(p))
		return;
	Vector vox = p - extent.pMin;
	float fx = vox.x / step;
	float fy = vox.y / step;
	float fz = vox.z / step;
	int vx = int(fx); int vy = int(fy); int vz = int(fz);
	int xmin = fx - vx > 0.5f ? vx : vx - 1;
	int ymin = fy - vy > 0.5f ? vy : vy - 1;
	int zmin = fz - vz > 0.5f ? vz : vz - 1;
	xmin = Clamp(xmin, 0, nx - 1);
	ymin = Clamp(ymin, 0, ny - 1);
	zmin = Clamp(zmin, 0, nz - 1);
	int xmax = xmin + 1; int ymax = ymin + 1; int zmax = zmin + 1;
	xmax = Clamp(xmax, 0, nx - 1);
	ymax = Clamp(ymax, 0, ny - 1);
	zmax = Clamp(zmax, 0, nz - 1);

	Vector offset = extent.Offset(p);

    //std::cout<<offset[0]<<" "<<offset[1]<<" "<<noise2d.Evaluate(offset[1],offset[2])<<std::endl;
    float p_radius=0.5f*(noise2d.Evaluate(offset[2],offset[1])/1.6f)-0.4f;
    float t_radius=radius+p_radius;
    //std::cout<<"Using radius "<<t_radius<<std::endl;
	//	scan all the grids within the range above
	float radiusSq = t_radius * t_radius;
	//	sigam in gaussian kernel
	float sigma = t_radius;
	float halfInvSigmaSq = .5f / sigma / sigma;
	float sum = 0.f;
	for (int i = xmin; i <= xmax; i++)
		for (int j = ymin; j <= ymax; j++)
			for (int k = zmin; k <= zmax; k++)
			{
				AuroraVoxel voxel = grids[k*nx*ny + j*nx + i];
				int photonNum = voxel.PhotonNum();
				for (int s = 0; s < photonNum; s++)
				{
					//	get the position of the photon
					AuroraPhoton photon = voxel.photons[s];
					Vector v = photon.p - p;
					float lenSq = v.LengthSquared();
					if (lenSq >= radiusSq)
						continue;
					//	compute the gaussian weight
					float weight = expf(-lenSq * halfInvSigmaSq);
					r += (weight * photon.r);
					g += (weight * photon.g);
					b += (weight * photon.b);
					sum += weight;
				}
			}
	if (sum > 0.f)
	{
		r /= sum;
		g /= sum;
		b /= sum;
	}
	//	otherwise r = g = b = 0.f
}

float AuroraGrid::LoadFactor()
{
	int count = 0;
	for (int i = 0; i < nx * ny * nz; i++){
		if (grids[i].PhotonNum() > 0)
			count++;
	}
  return count * 1.f / nx / ny / nz;
}

AuroraDensity::AuroraDensity(const Spectrum &sa, const Spectrum &ss, float gg, const BBox &e,
  const Transform &VolumeToWorld, float aa, float bb,
  const Vector &up, int x, int y, int z, float radius, const float *d,
  float cr, float cg, float cb, float is,
  const string rcolor, const string gcolor, const string bcolor, const string intensity,
  const Vector &b, int bn, float a, float l, float dd, float dA, float ee,
  float _persistence, float _frequency, int _level, int _seed, float _magnitude, float _perturb )
    :sig_a(sa), sig_s(ss), g(gg), extent(e),
      WorldToVolume(Inverse(VolumeToWorld)), a(aa), b(bb), 
      nx(x), ny(y), nz(z), grid(e, nx, ny, nz, radius),
      B(b), beamNum(bn), alphaD(a), L(l), dt(dd), deltaAlpha(dA), eleThreshold(ee),
      persistence(_persistence), frequency(_frequency), level(_level), seed(_seed), magnitude(_magnitude), perturb(_perturb)
{
  upDir = Normalize(up);
  eleDensity = new float[nx*ny*nz];
      memcpy(eleDensity, d, nx * ny * nz * sizeof(float));
  //  color and intensity
  //  read aurora color information
  std::ifstream fin;
  float height, value;
  //  r color
  fin.open(rcolor.c_str());
  while (fin >> height)
  {
    fin >> value;
    auroraColor[0].Add_Sample(height, value*cr);
  }
  auroraColor[0].Build();
  fin.close();
  //  g color
  fin.open(gcolor.c_str());
  while (fin >> height)
  {
    fin >> value;
    auroraColor[1].Add_Sample(height, value*cg);
  }
  auroraColor[1].Build();
  fin.close();
  //  b color
  fin.open(bcolor.c_str());
  while (fin >> height)
  {
    fin >> value;
    auroraColor[2].Add_Sample(height, value*cb);
  }
  auroraColor[2].Build();
  fin.close();
  //  intensity
  fin.open(intensity.c_str());
  while (fin >> height)
  {
    fin >> value;
    auroraIntensity.Add_Sample(height, value*is);
  }
  auroraIntensity.Build();
  fin.close();
  //  generate photons
  GeneratePhotons();
}

AuroraDensity::~AuroraDensity()
{
	delete []eleDensity;
}

float AuroraDensity::EleDensity(const Point &Pobj) const
{
	if (!extent.Inside(Pobj)) return 0;
    // Compute voxel coordinates and offsets for _Pobj_
    Vector vox = extent.Offset(Pobj);
    vox.x *= nx;
    vox.y *= ny;
    vox.z *= nz;
    int vx = (int)vox.x;
	int vy = (int)vox.y;
	int vz = (int)vox.z;
	if (vx == nx || vy == ny || vz == nz)
		return 0;
    float dx = vox.x - vx, dy = vox.y - vy, dz = vox.z - vz;
    // Trilinearly interpolate density values to compute local density
    float d00 = Lerp(dx, D(vx, vy, vz),     D(vx+1, vy, vz));
    float d10 = Lerp(dx, D(vx, vy+1, vz),   D(vx+1, vy+1, vz));
    float d01 = Lerp(dx, D(vx, vy, vz+1),   D(vx+1, vy, vz+1));
    float d11 = Lerp(dx, D(vx, vy+1, vz+1), D(vx+1, vy+1, vz+1));
    float d0 = Lerp(dy, d00, d10);
    float d1 = Lerp(dy, d01, d11);
    return Lerp(dz, d0, d1);
}

float AuroraDensity::Density(const Point &Pobj) const
{
    if (!extent.Inside(Pobj)) return 0;
    float height = Dot(Vector(Pobj), upDir);
    return a * expf(-b * height);
}

Spectrum AuroraDensity::Lve(const Point &p, const Vector &, float) const
{
	const Point Pobj = WorldToVolume(p);
	float den = Density(Pobj);
	float rgb[3];
	grid.SearchInGrid(Pobj, rgb[0], rgb[1], rgb[2]);
	
	rgb[0] *= den;
	rgb[1] *= den;
	rgb[2] *= den;
	
	return Spectrum::FromRGB(rgb);
}

Spectrum AuroraDensity::tau(const Ray &r, float stepSize, float u) const 
{
    float t0, t1;
    float length = r.d.Length();
    if (length == 0.f) return 0.f;
    Ray rn(r.o, r.d / length, r.mint * length, r.maxt * length, r.time);
    if (!IntersectP(rn, &t0, &t1)) return 0.;
    Spectrum tau(0.);
    t0 += u * stepSize;
    while (t0 < t1) {
        tau += sigma_t(rn(t0), -rn.d, r.time);
        t0 += stepSize;
    }
    return tau * stepSize;
}

void AuroraDensity::GeneratePhotons()
{
    Perlin start_point_noise(0.5f,0.25,3,1);
    Perlin noise(persistence,frequency,level,seed);
	//	generate photons
	int count = 0;
	Vector vox = extent.pMax - extent.pMin;
	//	build a local geomagnetic coordinate
	Vector u = B / B.Length();
	Vector v = Cross(B, Vector(1, 0, 0));
	v /= v.Length();
	Vector w = Cross(u, v);
	w /= w.Length();
	int maxEleNum = int(1.f / dt);
	//	count the photon number
	int photonNum = 0;
	int percent = 0;
	while (count < beamNum)
	{
		float dx = rand() * 1.f / RAND_MAX;
		float dy = rand() * 1.f / RAND_MAX;
		float dz = rand() * 1.f / RAND_MAX;
		Vector offset = vox;
		offset.x *= dx;
		offset.y *= dy;
		offset.z *= dz;
		Point start = extent.pMin + offset;
        start+=upDir*(30.0f*(start_point_noise.Evaluate(dz)-0.0f));
		float density = EleDensity(start);
		if (density > eleThreshold)
		{
			count++;
			if (count * 100.f / beamNum >= percent)
			{
				std::cout <<"\r"<< percent << "%\t" << std::flush;
				percent++;
			}
			//	start to simulate a new beam
			Point p = start;
			for (int i = 0; i < maxEleNum; i++)
			{
				//	generate a deflection point
				float t = rand() * 1.f / RAND_MAX;
				float alpha = Radians(alphaD - t * deltaAlpha);
				float beta = 2 * M_PI * rand() * 1.f / RAND_MAX;
				float len = L * dt * t;
				p += (len * u);
				p += (tanf(alpha) * len * (cosf(beta) * v + sinf(beta) * w));
				if (EleDensity(p) <= eleThreshold)
					break;
				else
				{
					//	add a new photon
					float h = Dot(Vector(p), upDir);
                    //h+=80.0f*(noise.Evaluate(dz)-0.5f);
					float h0 = Dot(Vector(extent.pMin), upDir);
					float intensity = auroraIntensity.Evaluate(h - h0);
					float r = auroraColor[0].Evaluate(h+magnitude*(noise.Evaluate(dz)+perturb)) * intensity;
					float g = auroraColor[1].Evaluate(h+magnitude*(noise.Evaluate(dz)+perturb)) * intensity;
					float b = auroraColor[2].Evaluate(h+magnitude*(noise.Evaluate(dz)+perturb)) * intensity;
					AuroraPhoton photon(p, r, g, b);
					grid.AddPhoton(photon);
					photonNum++;
				}
			}
		}
	}
	std::cout << std::endl << "photon number: " << photonNum << std::endl;
	//	analyse the grids
	std::cout << "load factor: " << grid.LoadFactor() << std::endl;
}

// AuroraDensity Method Definitions
AuroraDensity *CreateAuroraVolumeRegion(const Transform &volume2world,
        const ParamSet &params) {
    // Initialize common volume region parameters
    Spectrum sigma_a = params.FindOneSpectrum("sigma_a", 0.);
    Spectrum sigma_s = params.FindOneSpectrum("sigma_s", 0.);
    float g = params.FindOneFloat("g", 0.);
    Point p0 = params.FindOnePoint("p0", Point(0,0,0));
    Point p1 = params.FindOnePoint("p1", Point(1,1,1));
	float a = params.FindOneFloat("a", 1.);
    float b = params.FindOneFloat("b", 1.);
    Vector up = params.FindOneVector("updir", Vector(0,1,0));
	float radius = params.FindOneFloat("radius", 1.);
	int nitems;
    const float *data = params.FindFloat("density", &nitems);
    if (!data) {
        Error("No \"density\" values provided for volume grid?");
        return NULL;
    }
    int nx = params.FindOneInt("nx", 1);
    int ny = params.FindOneInt("ny", 1);
    int nz = params.FindOneInt("nz", 1);
    if (nitems != nx*ny*nz) {
        Error("VolumeGridDensity has %d density values but nx*ny*nz = %d",
            nitems, nx*ny*nz);
        return NULL;
    }
	float cr = params.FindOneFloat("rscale", 1.f);
	float cg = params.FindOneFloat("gscale", 1.f);
	float cb = params.FindOneFloat("bscale", 1.f);
	float is = params.FindOneFloat("iscale", 1.f);
	string rcolor = params.FindOneFilename("aurora_r", "aurora_r.txt.even");
	string gcolor = params.FindOneFilename("aurora_g", "aurora_g.txt.even");
	string bcolor = params.FindOneFilename("aurora_b", "aurora_b.txt.even");
	string intensity = params.FindOneFilename("aurora_intensity", "rays.txt.even");
	Vector B = params.FindOneVector("B", Vector(1, -1, 0));
	int bn = params.FindOneInt("beamNumber", 1500);
	float alphaD = params.FindOneFloat("alphaD", 3);
	float L = params.FindOneFloat("length", 175);
	float dt = params.FindOneFloat("dt", 1.f / 300);
	float dA = params.FindOneFloat("dA", 0.86f);
	float ee = params.FindOneFloat("threshold", 0.03f);
    float persistence = params.FindOneFloat("persistence", 0.25);
    float frequency = params.FindOneFloat("frequency", 0.2);
    int level = params.FindOneInt("level", 2);
    int seed = params.FindOneInt("seed", 0);
    float magnitude = params.FindOneFloat("magnitude", 60.0f);
    float perturb = params.FindOneFloat("perturb", 0.0f);
	//	TODO: what is our parameter for AuroraDensity?
	return new AuroraDensity(sigma_a, sigma_s, g, BBox(p0, p1),
        volume2world, a, b, up, nx, ny, nz, radius, data, cr, cg, cb, is, rcolor, gcolor, bcolor, intensity,
		B, bn, alphaD, L, dt, dA, ee,
        persistence, frequency, level, seed, magnitude, perturb);
}
