#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <limits>
#include <string>
#include "gif.h"  // taken from https://github.com/charlietangora/gif-h 

#define M_PI 3.14159265358979323846

using namespace std;

// 3D vector class that handles color values.
class Vector3D {
public:
    double x, y, z;
    Vector3D() : x(0), y(0), z(0) {}
    Vector3D(double x, double y, double z) : x(x), y(y), z(z) {}

    // Allow negation of a vector.
    Vector3D operator-() const { return Vector3D(-x, -y, -z); }
    Vector3D operator+(const Vector3D& v) const { return Vector3D(x+v.x, y+v.y, z+v.z); }
    Vector3D operator-(const Vector3D& v) const { return Vector3D(x-v.x, y-v.y, z-v.z); }
    Vector3D operator*(double s) const { return Vector3D(x*s, y*s, z*s); }
    Vector3D operator/(double s) const { return Vector3D(x/s, y/s, z/s); }
    double dot(const Vector3D& v) const { return x*v.x + y*v.y + z*v.z; }
    double norm() const { return sqrt(x*x + y*y + z*z); }
    Vector3D normalize() const { return *this / norm(); }
};


class Sphere {
public:
    Vector3D center, color;
    double radius, specular, reflective;
    Sphere(Vector3D c, double r, Vector3D col, double spec, double refl)
        : center(c), radius(r), color(col), specular(spec), reflective(refl) {}
};


class Light {
    public: 
    Vector3D direction;
    Vector3D position;
    double intensity;
    string type;
};


vector<Sphere> spheres = {
    Sphere(Vector3D(0, -1, 3), 1, Vector3D(255, 0, 150), 500, 0.6),
    Sphere(Vector3D(2, 0, 4), 1, Vector3D(0, 255, 150), 500, 0.3),
    Sphere(Vector3D(-2, 0, 4), 1, Vector3D(150, 100, 255), 10, 0.4)
};

vector<Light> lights = {
    {Vector3D(0, 0, 0), Vector3D(0, 0, 0), 0.2, "ambient"},
    {Vector3D(-1, -1, -1).normalize(), Vector3D(0, 0, 0), 0.6, "directional"},
    {Vector3D(0, 0, 0), Vector3D(2, 1, 0), 0.3, "point"}
};

// constant and camera origin.
const int width = 500, height = 500;
const double viewport_size = 1.0, projection_plane_z = 1.0;
const Vector3D background_color(0, 0, 0);
Vector3D O(0, 0, 0); //camara


Vector3D canvasToViewport(int x, int y) {
    return Vector3D(x * viewport_size / width, y * viewport_size / height, projection_plane_z);
}


double intersectRaySphere(Vector3D O, Vector3D D, Sphere sphere) {
    Vector3D CO = O - sphere.center;
    double a = D.dot(D);
    double b = 2 * CO.dot(D);
    double c = CO.dot(CO) - sphere.radius * sphere.radius;
    double discriminant = b*b - 4*a*c;
    if (discriminant < 0)
        return numeric_limits<double>::infinity();
    double t1 = (-b + sqrt(discriminant)) / (2*a);
    double t2 = (-b - sqrt(discriminant)) / (2*a);
    return (t2 > 0) ? t2 : t1;
}


double intersectRayPlane(Vector3D O, Vector3D D, double planeY) {
    if (fabs(D.y) < 1e-6)
        return numeric_limits<double>::infinity();
    double t = (planeY - O.y) / D.y;
    return (t > 0) ? t : numeric_limits<double>::infinity();
}

Vector3D reflectRay(Vector3D incident, Vector3D N) {
    return incident - N * (2 * incident.dot(N));
}


Vector3D computeLighting(Vector3D P, Vector3D N, Vector3D V, double specular) {
    double intensity = 0.0;
    for (auto &light : lights) {
        Vector3D L;
        double t_max;
        if (light.type == "ambient") {
            intensity += light.intensity;
        } else {
            if (light.type == "point") {
                L = (light.position - P).normalize();
                t_max = 1.0;
            } else { 
                L = light.direction;
                t_max = numeric_limits<double>::infinity();
            }

            double t_shadow = intersectRaySphere(P, L, spheres[0]);
            if (t_shadow < numeric_limits<double>::infinity())
                continue;
            double dotNL = max(0.0, N.dot(L));
            intensity += light.intensity * dotNL;
            if (specular != -1) {
                Vector3D R = reflectRay(L, N);
                double r_dot_v = max(0.0, R.dot(V));
                if (r_dot_v > 0)
                    intensity += light.intensity * pow(r_dot_v / (R.norm() * V.norm()), specular);
            }
        }
    }
    return Vector3D(intensity, intensity, intensity);
}

// trace ray and returns color
Vector3D traceRay(Vector3D O, Vector3D D, int depth) {
    if (depth > 3)
        return background_color;
    double min_t = numeric_limits<double>::infinity();
    bool hitPlane = false;
    Sphere* closest_sphere = nullptr;
    for (auto &sphere : spheres) {
        double t = intersectRaySphere(O, D, sphere);
        if (t < min_t) {
            min_t = t;
            closest_sphere = &sphere;
            hitPlane = false;
        }
    }
    double t_plane = intersectRayPlane(O, D, -1);
    if (t_plane < min_t) {
        min_t = t_plane;
        hitPlane = true;
        closest_sphere = nullptr;
    }
    if (min_t == numeric_limits<double>::infinity())
        return background_color;
    Vector3D P = O + D * min_t;
    if (hitPlane) {
        Vector3D N(0, 1, 0);
        Vector3D V = -D;
        Vector3D lighting = computeLighting(P, N, V, 1000);
        Vector3D plane_color(255, 255, 0); // Yellow floor
        return plane_color * lighting.x;
    }
    
    Vector3D N = (P - closest_sphere->center).normalize();
    Vector3D V = -D;
    Vector3D lighting = computeLighting(P, N, V, closest_sphere->specular);
    Vector3D local_color = closest_sphere->color * lighting.x;
    double r = closest_sphere->reflective;
    if (depth <= 0 || r <= 0)
        return local_color;
    Vector3D R = reflectRay(D, N);
    Vector3D reflected_color = traceRay(P, R, depth - 1);
    return local_color * (1 - r) + reflected_color * r;
}

// Renders the scene into a frame buffer.
void renderFrameToBuffer(unsigned char* frameBuffer) {
    for (int j = 0; j < height; j++) {
        int canvas_y = height/2 - j;
        for (int i = 0; i < width; i++) {
            int canvas_x = i - width/2;
            Vector3D color = traceRay(O, canvasToViewport(canvas_x, canvas_y), 0);
            int r = (int)min(255.0, color.x);
            int g = (int)min(255.0, color.y);
            int b = (int)min(255.0, color.z);
            int index = 4 * (j * width + i);
            frameBuffer[index + 0] = (unsigned char)r;
            frameBuffer[index + 1] = (unsigned char)g;
            frameBuffer[index + 2] = (unsigned char)b;
            frameBuffer[index + 3] = 255;
        }
    }
}

int main() {
    const int frame_count = 30;
    const float amplitude = 0.3;
    float base_y[3] = {-1, 0, 0};


    GifWriter writer;
    const int delay = 10;
    GifBegin(&writer, "Raytracing animation.gif", width, height, delay);

    unsigned char* frameBuffer = new unsigned char[width * height * 4];

    // Animate and render each frame.
    for (int frame = 0; frame < frame_count; frame++) {
        for (int i = 0; i < 3; i++) {
            spheres[i].center.y = base_y[i] + amplitude * sin(frame * 2 * M_PI / frame_count);
        }
        renderFrameToBuffer(frameBuffer);
        GifWriteFrame(&writer, frameBuffer, width, height, delay);

        // Save frame as a simple PPM image.
        string filename = "frame" + to_string(frame) + ".ppm";
        ofstream ofs(filename);
        ofs << "P3" << width << " " << height << "255";
        for (int j = 0; j < height; j++) {
            for (int i = 0; i < width; i++) {
                int index = 4 * (j * width + i);
                ofs << (int)frameBuffer[index] << " "
                    << (int)frameBuffer[index+1] << " "
                    << (int)frameBuffer[index+2] << " ";
            }
        }
        ofs.close();
        cout << "Frame " << frame + 1<< " generated" << endl;
    }

    GifEnd(&writer);
    delete[] frameBuffer;
    return 0;
}
