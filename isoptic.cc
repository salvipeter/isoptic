#include <cmath>
#include <unordered_set>
#include <dc.hh>
#include <geometry.hh>

using namespace Geometry;

double area(const TriMesh &mesh, const Point3D &p) {
  // As in: G. Csima, J. Szirmai: Isoptic surfaces of polyhedra.
  //        CAGD 47, pp. 55-60, 2016.
  double result = 0.0;
  for (const auto &tri : mesh.triangles()) {
    double Omega = 2 * M_PI;
    for (size_t j = 0; j < 3; ++j) {
      auto v = mesh[tri[j]] - p;
      auto v_1 = mesh[tri[(j+2)%3]] - p;
      auto v1 = mesh[tri[(j+1)%3]] - p;
      auto c1 = (v_1 ^ v).normalize();
      auto c2 = (v ^ v1).normalize();
      Omega -= std::acos(std::min(std::max(c1 * c2, -1.0), 1.0));
    }
    result += Omega;
  }
  return result / 2;
}

using Normals = std::vector<Vector3D>;

std::vector<Normals> findNormals(const TriMesh &mesh) {
  std::vector<Normals> result(mesh.points().size());
  for (const auto &tri : mesh.triangles()) {
    auto n = ((mesh[tri[1]] - mesh[tri[0]]) ^ (mesh[tri[2]] - mesh[tri[0]])).normalize();
    for (size_t i = 0; i < 3; ++i)
      result[tri[i]].push_back(n);
  }
  return result;
}

bool isSilhouette(const Point3D &q, const Normals &ns, const Point3D &p) {
  if (ns.size() == 1)
    return true;
  auto u = q - p;
  char c = 0;
  for (const auto &n : ns) {
    if (n * u < 0) {
      if (c == 2)
        return true;
      c = 1;
    } else {
      if (c == 1)
        return true;
      c = 2;
    }
  }
  return false;
}

double area_concave(const TriMesh &mesh, const std::vector<Normals> &normals, const Point3D &p) {
  // Alternative definition, also for concave meshes
  std::unordered_set<size_t> silhouette;
  size_t n = mesh.points().size();
  for (size_t i = 0; i < n; ++i)
    if (isSilhouette(mesh[i], normals[i], p))
      silhouette.insert(i);

  double result = 0.0;
  for (auto i = silhouette.begin(); i != silhouette.end(); ++i)
    for (auto j = i; j != silhouette.end(); ++j) {
      const auto &q1 = mesh[*i], &q2 = mesh[*j];
      auto d = std::acos((q1 - p).normalize() * (q2 - p).normalize());
      if (d > result)
        result = d;
    } 

  return result;
}

std::array<Point3D, 2> boundingBox(const TriMesh &mesh, double scaling) {
  Point3D boxmin, boxmax;
  const auto &points = mesh.points();
  boxmin = boxmax = points[0];
  for (const auto &p : points)
    for (int i = 0; i < 3; ++i) {
      boxmin[i] = std::min(boxmin[i], p[i]);
      boxmax[i] = std::max(boxmax[i], p[i]);
    }
  auto mean = (boxmin + boxmax) / 2;
  boxmin = mean + (boxmin - mean) * scaling;
  boxmax = mean + (boxmax - mean) * scaling;
  return { boxmin, boxmax };
}

int main(int argc, char **argv) {
  if (argc < 2 || argc > 5) {
    std::cerr << "Usage: " << argv[0] << " <input.obj> [scaling] [res] [alpha]" << std::endl;
    return 1;
  }

  auto mesh = TriMesh::readOBJ(argv[1]);
  double scaling = 2.5;
  if (argc >= 3)
    scaling = std::strtod(argv[2], nullptr);
  size_t res = 30;
  if (argc >= 4)
    res = std::atoi(argv[3]);
  double alpha = M_PI / 2;
  if (argc >= 5)
    alpha = std::strtod(argv[4], nullptr);

  auto normals = findNormals(mesh);
  auto f = [&](const DualContouring::Point3D &p) {
    // return area(mesh, { p[0], p[1], p[2] });
    return area_concave(mesh, normals, { p[0], p[1], p[2] });
  };
  auto bbox = boundingBox(mesh, scaling);
  std::array<DualContouring::Point3D, 2> dc_bbox = { {
      { bbox[0][0], bbox[0][1], bbox[0][2] },
      { bbox[1][0], bbox[1][1], bbox[1][2] }
    } };

  auto quad = DualContouring::isosurface(f, alpha, dc_bbox, { res, res, res });
  quad.writeOBJ("output.obj");
}
