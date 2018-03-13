#include "bsdf.h"

#include <algorithm>
#include <iostream>
#include <utility>


using std::min;
using std::max;
using std::swap;

namespace CMU462 {

void make_coord_space(Matrix3x3& o2w, const Vector3D& n) {
  Vector3D z = Vector3D(n.x, n.y, n.z);
  Vector3D h = z;
  if (fabs(h.x) <= fabs(h.y) && fabs(h.x) <= fabs(h.z))
    h.x = 1.0;
  else if (fabs(h.y) <= fabs(h.x) && fabs(h.y) <= fabs(h.z))
    h.y = 1.0;
  else
    h.z = 1.0;

  z.normalize();
  Vector3D y = cross(h, z);
  y.normalize();
  Vector3D x = cross(z, y);
  x.normalize();

  o2w[0] = x;
  o2w[1] = y;
  o2w[2] = z;
}

// Diffuse BSDF //

Spectrum DiffuseBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return albedo * (1.0 / PI);
}

Spectrum DiffuseBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  // TODO (PathTracer):
  // Implement DiffuseBSDF

	*wi = sampler.get_sample(pdf);

	return f(wo, *wi);

  //return Spectrum();
}

// Mirror BSDF //

Spectrum MirrorBSDF::f(const Vector3D& wo, const Vector3D& wi) {

	return reflectance*(1/fabs(wo.z));

}

Spectrum MirrorBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  // TODO (PathTracer):
  // Implement MirrorBSDF
	//std::cout << "mirror" << std::endl;
	//std::cout << "wo: " << wo << std::endl;
	reflect(wo, wi);
	*pdf = 1;

	Spectrum result = reflectance*(1 / fabs(wo.z));

	return result;
}

// Glossy BSDF //

/*
Spectrum GlossyBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum GlossyBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  *pdf = 1.0f;
  return reflect(wo, wi, reflectance);
}
*/

// Refraction BSDF //

Spectrum RefractionBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum RefractionBSDF::sample_f(const Vector3D& wo, Vector3D* wi,
                                  float* pdf) {
  // TODO (PathTracer):
  // Implement RefractionBSDF
    //std::cout << "glass" << std::endl;
	refract(wo, wi, ior);
	*pdf = 1;
	return transmittance*(1 / fabs(wo.z));
}

// Glass BSDF //

Spectrum GlassBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum GlassBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  // TODO (PathTracer):
  // Compute Fresnel coefficient and either reflect or refract based on it.
	//std::cout << "ior: " << ior << std::endl;
	//case 1 total internal reflection
	if (!refract(wo, wi, ior)) {
		reflect(wo,wi);
		*pdf = 1;
		return reflectance*(1 / fabs(wo.z));
	}

	//case 2 reflect or refract a ray
	double r0, r1;
	//std::cout << "cos in: " << cos_theta(wo) << std::endl;
	//std::cout << "cos out: " << cos_theta(*wi) << std::endl;
	r0 = ((double)ior*abs_cos_theta(wo) - abs_cos_theta(*wi)) / ((double)ior*abs_cos_theta(wo) + abs_cos_theta(*wi));
	r1 = (abs_cos_theta(wo) - (double)ior*abs_cos_theta(*wi)) / (abs_cos_theta(wo) + (double)ior*abs_cos_theta(*wi));
	//std::cout << "ro: " << r0 << std::endl;
	//std::cout << "r1: " << r1 << std::endl;

	double f = (pow(r0, 2) + pow(r1, 2)) / 2;
	//std::cout << "f: " << f << std::endl;

	f = clamp(f, 0., 1.);

	//std::cout << "f: " << f << std::endl;

	double randomNum = (double)(std::rand()) / RAND_MAX;

	if (randomNum < f) {
		
		*pdf = (float)f;
		//std::cout << "****refrected*****" << std::endl;
		return *pdf * transmittance*(1 / fabs(wo.z));
	}
	else {
		*pdf = (double)(1 - f);
		reflect(wo, wi);
		//std::cout << "reflected" << std::endl;
		return *pdf * reflectance*(1 / fabs(wo.z));
	}
    
}

void BSDF::reflect(const Vector3D& wo, Vector3D* wi) {
  // TODO (PathTracer):
  // Implement reflection of wo about normal (0,0,1) and store result in wi.
	//Vector3D n = Vector3D(0., 0., 1.);

	//*wi = (-1 * wo + 2 * dot(wo, n) * n).norm();

	//wi->unit();
	//std::cout << "wi: " << *wi << std::endl;
	*wi = Vector3D(-wo.x, -wo.y, wo.z);
	//wi->unit();
}

bool BSDF::refract(const Vector3D& wo, Vector3D* wi, float ior) {
  // TODO (PathTracer):
  // Use Snell's Law to refract wo surface and store result ray in wi.
  // Return false if refraction does not occur due to total internal reflection
  // and true otherwise. When dot(wo,n) is positive, then wo corresponds to a
  // ray entering the surface through vacuum.
	Vector3D n = Vector3D(0., 0., 1.);

	//float a = 1. / ior;

	if (dot(wo, n) >= 0.) {
		//from vacumm
		//if ( sin_theta(wo) > ior) { //theta larger than the critical angle
		//	return false;
		//}
		//*wi = Vector3D(-sin_theta(wo)*cos_phi(wo),
		//	           -sin_theta(wo)*sin_phi(wo),
		//	           -sqrt(1 - sin_theta(wo)*sin_theta(wo) / (ior*ior)));	
		float b = 1. / ior;
		*wi = Vector3D(-wo.x/ior, -wo.y/ior,
			           -sqrt(1 - b * b + b * b * wo.z * wo.z));
	}
	else {
		//enter vaccum
		float a = ior;
		if ((1 - a * a + a * a * wo.z * wo.z) < 0) { // theta larger than the critical angle
			printf("internal\n");
			return false;
		}
		//*wi = Vector3D(-sin_theta(wo)*cos_phi(wo),
		//	           -sin_theta(wo)*sin_phi(wo),
		//	           -sqrt(1 - sin_theta(wo)*sin_theta(wo)*(ior*ior)));
		*wi = Vector3D(-wo.x*ior, -wo.y*ior,
			           sqrt(1 - a * a + a * a * wo.z * wo.z));
		//printf("here\n");
		//std::cout << "wo: " << wo << std::endl;
		//std::cout << "wi: " << 1 - sin_theta(wo)*sin_theta(wo) * (ior*ior) << std::endl;
	}

    return true;
}

// Emission BSDF //

Spectrum EmissionBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum EmissionBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  *wi = sampler.get_sample(pdf);
  return Spectrum();
}

}  // namespace CMU462
