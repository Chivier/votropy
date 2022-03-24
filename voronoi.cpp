#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <vector>

// CGAL
#include <CGAL/Convex_hull_3/dual/halfspace_intersection_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/point_generators_3.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_3.h>
#include <fstream>
#include <list>
#include <vector>

#include <CGAL/Convex_hull_3/dual/halfspace_intersection_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/point_generators_3.h>
#include <list>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Plane_3 Plane;
typedef Kernel::Vector_3 Vector;
typedef Kernel::Point_3 Point;
typedef CGAL::Surface_mesh<Point> Surface_mesh;

// compute the tangent plane of a point
Plane tangent_plane(Point const &p1, Point const &p2) {
    Point p{(p1.x() + p2.x()) / 2, (p1.y() + p2.y()) / 2, (p1.z() + p2.z()) / 2};
    Vector v(p.x() - p1.x(), p.y() - p1.y(), p.z() - p1.z());
    v = v / sqrt(v.squared_length());
    Plane plane(p, v);
    return plane;
}

int main(int argc, char **argv) {
    std::vector<Point> points;

    std::fstream input;
    if (argc == 1) {
        input.open("test.txt");
    } else {
        input.open(argv[1]);
    }
    if (!input.is_open()) {
        std::cerr << "Error: cannot open file " << argv[1] << std::endl;
        return 1;
    }

    double xmin = std::numeric_limits<double>::infinity();
    double ymin = std::numeric_limits<double>::infinity();
    double zmin = std::numeric_limits<double>::infinity();
    double xmax = -std::numeric_limits<double>::infinity();
    double ymax = -std::numeric_limits<double>::infinity();
    double zmax = -std::numeric_limits<double>::infinity();

    int point_nums = 0;
    input >> point_nums;

    for (int point_index = 0; point_index < point_nums; point_index++) {
        double x, y, z;
        input >> x >> y >> z;
        points.push_back(Point{x, y, z});
        if (x > xmax)
            xmax = x;
        if (x < xmin)
            xmin = x;
        if (y > ymax)
            ymax = y;
        if (y < ymin)
            ymin = y;
        if (z > zmax)
            zmax = z;
        if (z < zmin)
            zmin = z;
    }
    double epsilon = 10;
    xmax += epsilon;
    xmin -= epsilon;
    ymax += epsilon;
    ymin -= epsilon;
    zmax += epsilon;
    zmin -= epsilon;

    std::map<int, int> face_count;
    std::map<int, int> segment_count;
    std::map<int, int> point_count;

    for (size_t index = 0; index < point_nums; ++index) {
        std::list<Plane> planes;
        // Intermediate planes
        for (size_t subindex = 0; subindex < point_nums; ++subindex) {
            if (index != subindex) {
                planes.push_back(tangent_plane(points[index], points[subindex]));
            }
        }
        // Boundary planes
        planes.push_back(Plane(Point(xmin, 0, 0), Vector(-1.0, 0, 0)));
        planes.push_back(Plane(Point(xmax, 0, 0), Vector(1.0, 0, 0)));
        planes.push_back(Plane(Point(0, ymin, 0), Vector(0, -1.0, 0)));
        planes.push_back(Plane(Point(0, ymax, 0), Vector(0, 1.0, 0)));
        planes.push_back(Plane(Point(0, 0, zmin), Vector(0, 0, -1.0)));
        planes.push_back(Plane(Point(0, 0, zmax), Vector(0, 0, 1.0)));

        // define polyhedron to hold the intersection
        Surface_mesh chull;
        // compute the intersection
        // if no point inside the intersection is provided, one
        // will be automatically found using linear programming
        CGAL::halfspace_intersection_3(planes.begin(), planes.end(), chull, points[index]);

        if (face_count[num_faces(chull)] == 0) {
            face_count[num_faces(chull)] = 1;
        } else {
            face_count[num_faces(chull)]++;
        }

        if (segment_count[num_edges(chull)] == 0) {
            segment_count[num_edges(chull)] = 1;
        } else {
            segment_count[num_edges(chull)]++;
        }

        if (point_count[num_vertices(chull)] == 0) {
            point_count[num_vertices(chull)] = 1;
        } else {
            point_count[num_vertices(chull)]++;
        }
    }

    double voronoi_entropy1 = 0.0;
    double voronoi_entropy2 = 0.0;
    double voronoi_entropy3 = 0.0;
    int convex_hull_count = point_nums;

    for (auto item : face_count) {
        double p = 1.0 * item.second / convex_hull_count;
        voronoi_entropy1 += p * std::log(p);
    }

    for (auto item : segment_count) {
        double p = 1.0 * item.second / convex_hull_count;
        voronoi_entropy2 += p * std::log(p);
    }

    for (auto item : point_count) {
        double p = 1.0 * item.second / convex_hull_count;
        voronoi_entropy3 += p * std::log(p);
    }

    printf("Voronoi entropy(face) : %lf\n", -voronoi_entropy1);
    printf("Voronoi entropy(segs) : %lf\n", -voronoi_entropy2);
    printf("Voronoi entropy(point): %lf\n", -voronoi_entropy3);
    return 0;
}
