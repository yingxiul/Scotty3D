#include "triangle.h"

#include "CMU462/CMU462.h"
#include "GL/glew.h"

namespace CMU462 {
namespace StaticScene {

Triangle::Triangle(const Mesh* mesh, vector<size_t>& v) : mesh(mesh), v(v) {}
Triangle::Triangle(const Mesh* mesh, size_t v1, size_t v2, size_t v3)
    : mesh(mesh), v1(v1), v2(v2), v3(v3) {}

BBox Triangle::get_bbox() const {
  // TODO (PathTracer):
  // compute the bounding box of the triangle

	Vector3D p0 = mesh->positions[v1];
	Vector3D p1 = mesh->positions[v2];
	Vector3D p2 = mesh->positions[v3];
	BBox a(p0);
	a.expand(p1);
	a.expand(p2);
	return a;
}

bool Triangle::intersect(const Ray& r) const {
  // TODO (PathTracer): implement ray-triangle intersection

  //get 3 vertices position

	Vector3D p0 = mesh->positions[v1];
	Vector3D p1 = mesh->positions[v2];
	Vector3D p2 = mesh->positions[v3];

	Vector3D e1 = p1 - p0;
	Vector3D e2 = p2 - p0;
	Vector3D s = r.o - p0;

	double scalar = 1. / dot(cross(e1, r.d), e2);

	Vector3D vec = scalar * Vector3D(-dot(cross(s, e2), r.d), dot(cross(e1, r.d), s), -dot(cross(s, e2), e1));

	if (vec.z >= r.min_t && vec.z <= r.max_t &&
		vec.x >= 0 && vec.x <= 1 &&
		vec.y >= 0 && vec.y <= 1 &&
		1 - vec.x - vec.y >= 0 && 1 - vec.x - vec.y <= 1) {

		r.max_t = vec.z;
		return true;
	}

	return false;
}

bool Triangle::intersect(const Ray& r, Intersection* isect) const {
  // TODO (PathTracer):
  // implement ray-triangle intersection. When an intersection takes
  // place, the Intersection data should be updated accordingly

	Vector3D p0 = mesh->positions[v1];
	Vector3D p1 = mesh->positions[v2];
	Vector3D p2 = mesh->positions[v3];

	Vector3D e1 = p1 - p0;
	Vector3D e2 = p2 - p0;
	Vector3D s = r.o - p0;

	double scalar = 1. / dot(cross(e1, r.d), e2);

	Vector3D vec = scalar * Vector3D(-dot(cross(s, e2), r.d), dot(cross(e1, r.d), s), -dot(cross(s, e2), e1));

	if (vec.z >= r.min_t && vec.z <= r.max_t &&
		vec.x >= 0 && vec.x <= 1 &&
		vec.y >= 0 && vec.y <= 1 &&
		(1 - vec.x - vec.y) >= 0 && (1 - vec.x - vec.y) <= 1) {

		r.max_t = vec.z;

		isect->t = vec.z;
		isect->primitive = this;
		//check this
		isect->n = (1 - vec.x - vec.y) * mesh->normals[v1] + vec.x * mesh->normals[v2] + vec.y * mesh->normals[v3];
		if (dot(isect->n,r.d)>0) {
			isect->n = -1 * isect->n;
		}
		isect->bsdf = get_bsdf();

		return true;
	}

	return false;
}

void Triangle::draw(const Color& c) const {
  glColor4f(c.r, c.g, c.b, c.a);
  glBegin(GL_TRIANGLES);
  glVertex3d(mesh->positions[v1].x, mesh->positions[v1].y,
             mesh->positions[v1].z);
  glVertex3d(mesh->positions[v2].x, mesh->positions[v2].y,
             mesh->positions[v2].z);
  glVertex3d(mesh->positions[v3].x, mesh->positions[v3].y,
             mesh->positions[v3].z);
  glEnd();
}

void Triangle::drawOutline(const Color& c) const {
  glColor4f(c.r, c.g, c.b, c.a);
  glBegin(GL_LINE_LOOP);
  glVertex3d(mesh->positions[v1].x, mesh->positions[v1].y,
             mesh->positions[v1].z);
  glVertex3d(mesh->positions[v2].x, mesh->positions[v2].y,
             mesh->positions[v2].z);
  glVertex3d(mesh->positions[v3].x, mesh->positions[v3].y,
             mesh->positions[v3].z);
  glEnd();
}

}  // namespace StaticScene
}  // namespace CMU462
