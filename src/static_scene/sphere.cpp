#include "sphere.h"

#include <cmath>

#include "../bsdf.h"
#include "../misc/sphere_drawing.h"

namespace CMU462 {
namespace StaticScene {

bool Sphere::test(const Ray& r, double& t1, double& t2) const {
  // TODO (PathTracer):
  // Implement ray - sphere intersection test.
  // Return true if there are intersections and writing the
  // smaller of the two intersection times in t1 and the larger in t2.
  
	double a = r.d.norm2();
	double b = 2 * dot(r.o - o, r.d);
	double c = (r.o - o).norm2() - r2;
	double delta = b*b - 4 * a * c;

	if (delta < 0) {
		return false;
	}
	else {
		double min = (-b - sqrt(delta)) / (2 * a);
		if (min >= r.min_t && min <= r.max_t) {
			t1 = min;
		}
		else {
			t1 = (-b + sqrt(delta)) / (2 * a);
		}
		return true;
	}
}

bool Sphere::intersect(const Ray& r) const {
  // TODO (PathTracer):
  // Implement ray - sphere intersection.
  // Note that you might want to use the the Sphere::test helper here.

	//double t1 = 0;
	//double t2 = 0;
	//if (test(r, t1, t2) && t1 >= r.min_t && t1 <= r.max_t) {
	//	r.max_t = t1;
	//	return true;
	//}
	Intersection i;

	return intersect(r, &i);
}

bool Sphere::intersect(const Ray& r, Intersection* isect) const {
  // TODO (PathTracer):
  // Implement ray - sphere intersection.
  // Note again that you might want to use the the Sphere::test helper here.
  // When an intersection takes place, the Intersection data should be updated
  // correspondingly.

	double t1 = 0;
	double t2 = 0;
	if (test(r, t1, t2) && t1 >= r.min_t && t1 <= r.max_t) {
		r.max_t = t1;

		isect->t = t1;
		isect->primitive = this;
		//check this
		isect->n = normal(r.o + t1 * r.d);
		isect->bsdf = get_bsdf();
		return true;
	}
	return false;
}

void Sphere::draw(const Color& c) const { Misc::draw_sphere_opengl(o, r, c); }

void Sphere::drawOutline(const Color& c) const {
  // Misc::draw_sphere_opengl(o, r, c);
}

}  // namespace StaticScene
}  // namespace CMU462
