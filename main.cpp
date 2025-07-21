#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include <random>
#include <cmath>
#include <limits>

// Vector3 class for 3D operations
class Vec3 {
public:
    double x, y, z;
    
    Vec3() : x(0), y(0), z(0) {}
    Vec3(double x, double y, double z) : x(x), y(y), z(z) {}
    
    Vec3 operator-() const { return Vec3(-x, -y, -z); }
    Vec3 operator+(const Vec3& v) const { return Vec3(x + v.x, y + v.y, z + v.z); }
    Vec3 operator-(const Vec3& v) const { return Vec3(x - v.x, y - v.y, z - v.z); }
    Vec3 operator*(double t) const { return Vec3(x * t, y * t, z * t); }
    Vec3 operator*(const Vec3& v) const { return Vec3(x * v.x, y * v.y, z * v.z); }
    Vec3 operator/(double t) const { return *this * (1.0 / t); }
    
    Vec3& operator+=(const Vec3& v) { x += v.x; y += v.y; z += v.z; return *this; }
    Vec3& operator*=(double t) { x *= t; y *= t; z *= t; return *this; }
    Vec3& operator/=(double t) { return *this *= 1.0 / t; }
    
    double length() const { return std::sqrt(length_squared()); }
    double length_squared() const { return x*x + y*y + z*z; }
    
    Vec3 unit_vector() const { return *this / length(); }
    
    static Vec3 random() {
        static std::random_device rd;
        static std::mt19937 gen(rd());
        static std::uniform_real_distribution<> dis(0.0, 1.0);
        return Vec3(dis(gen), dis(gen), dis(gen));
    }
    
    static Vec3 random(double min, double max) {
        static std::random_device rd;
        static std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(min, max);
        return Vec3(dis(gen), dis(gen), dis(gen));
    }
    
    static Vec3 random_in_unit_sphere() {
        Vec3 p;
        do {
            p = Vec3::random(-1, 1);
        } while (p.length_squared() >= 1.0);
        return p;
    }
    
    static Vec3 random_unit_vector() {
        return random_in_unit_sphere().unit_vector();
    }
    
    bool near_zero() const {
        const auto s = 1e-8;
        return (std::fabs(x) < s) && (std::fabs(y) < s) && (std::fabs(z) < s);
    }
};

using Point3 = Vec3;
using Color = Vec3;

// Utility functions
inline Vec3 operator*(double t, const Vec3& v) { return v * t; }

inline double dot(const Vec3& u, const Vec3& v) {
    return u.x * v.x + u.y * v.y + u.z * v.z;
}

inline Vec3 cross(const Vec3& u, const Vec3& v) {
    return Vec3(u.y * v.z - u.z * v.y,
                u.z * v.x - u.x * v.z,
                u.x * v.y - u.y * v.x);
}

inline Vec3 reflect(const Vec3& v, const Vec3& n) {
    return v - 2 * dot(v, n) * n;
}

inline Vec3 refract(const Vec3& uv, const Vec3& n, double etai_over_etat) {
    auto cos_theta = std::fmin(dot(-uv, n), 1.0);
    Vec3 r_out_perp = etai_over_etat * (uv + cos_theta * n);
    Vec3 r_out_parallel = -std::sqrt(std::fabs(1.0 - r_out_perp.length_squared())) * n;
    return r_out_perp + r_out_parallel;
}

// Ray class
class Ray {
public:
    Point3 origin;
    Vec3 direction;
    
    Ray() {}
    Ray(const Point3& origin, const Vec3& direction) : origin(origin), direction(direction) {}
    
    Point3 at(double t) const {
        return origin + t * direction;
    }
};

// Hit record structure
struct HitRecord {
    Point3 p;
    Vec3 normal;
    double t;
    bool front_face;
    std::shared_ptr<class Material> mat_ptr;
    
    void set_face_normal(const Ray& r, const Vec3& outward_normal) {
        front_face = dot(r.direction, outward_normal) < 0;
        normal = front_face ? outward_normal : -outward_normal;
    }
};

// Abstract hittable class
class Hittable {
public:
    virtual ~Hittable() = default;
    virtual bool hit(const Ray& r, double t_min, double t_max, HitRecord& rec) const = 0;
};

// Material classes
class Material {
public:
    virtual ~Material() = default;
    virtual bool scatter(const Ray& r_in, const HitRecord& rec, Color& attenuation, Ray& scattered) const = 0;
};

class Lambertian : public Material {
public:
    Color albedo;
    
    Lambertian(const Color& a) : albedo(a) {}
    
    virtual bool scatter(const Ray& r_in, const HitRecord& rec, Color& attenuation, Ray& scattered) const override {
        (void)r_in; // Unused parameter
        auto scatter_direction = rec.normal + Vec3::random_unit_vector();
        
        if (scatter_direction.near_zero())
            scatter_direction = rec.normal;
        
        scattered = Ray(rec.p, scatter_direction);
        attenuation = albedo;
        return true;
    }
};

class Metal : public Material {
public:
    Color albedo;
    double fuzz;
    
    Metal(const Color& a, double f) : albedo(a), fuzz(f < 1 ? f : 1) {}
    
    virtual bool scatter(const Ray& r_in, const HitRecord& rec, Color& attenuation, Ray& scattered) const override {
        Vec3 reflected = reflect(r_in.direction.unit_vector(), rec.normal);
        scattered = Ray(rec.p, reflected + fuzz * Vec3::random_in_unit_sphere());
        attenuation = albedo;
        return (dot(scattered.direction, rec.normal) > 0);
    }
};

class Dielectric : public Material {
public:
    double ir; // Index of Refraction
    
    Dielectric(double index_of_refraction) : ir(index_of_refraction) {}
    
    virtual bool scatter(const Ray& r_in, const HitRecord& rec, Color& attenuation, Ray& scattered) const override {
        attenuation = Color(1.0, 1.0, 1.0);
        double refraction_ratio = rec.front_face ? (1.0/ir) : ir;
        
        Vec3 unit_direction = r_in.direction.unit_vector();
        double cos_theta = std::fmin(dot(-unit_direction, rec.normal), 1.0);
        double sin_theta = std::sqrt(1.0 - cos_theta*cos_theta);
        
        bool cannot_refract = refraction_ratio * sin_theta > 1.0;
        Vec3 direction;
        
        if (cannot_refract || reflectance(cos_theta, refraction_ratio) > random_double())
            direction = reflect(unit_direction, rec.normal);
        else
            direction = refract(unit_direction, rec.normal, refraction_ratio);
        
        scattered = Ray(rec.p, direction);
        return true;
    }
    
private:
    static double reflectance(double cosine, double ref_idx) {
        // Use Schlick's approximation for reflectance.
        auto r0 = (1-ref_idx) / (1+ref_idx);
        r0 = r0*r0;
        return r0 + (1-r0)*std::pow((1 - cosine),5);
    }
    
    static double random_double() {
        static std::random_device rd;
        static std::mt19937 gen(rd());
        static std::uniform_real_distribution<> dis(0.0, 1.0);
        return dis(gen);
    }
};

// Sphere class
class Sphere : public Hittable {
public:
    Point3 center;
    double radius;
    std::shared_ptr<Material> mat_ptr;
    
    Sphere() {}
    Sphere(Point3 cen, double r, std::shared_ptr<Material> m) : center(cen), radius(r), mat_ptr(m) {};
    
    virtual bool hit(const Ray& r, double t_min, double t_max, HitRecord& rec) const override {
        Vec3 oc = r.origin - center;
        auto a = r.direction.length_squared();
        auto half_b = dot(oc, r.direction);
        auto c = oc.length_squared() - radius*radius;
        
        auto discriminant = half_b*half_b - a*c;
        if (discriminant < 0) return false;
        auto sqrtd = std::sqrt(discriminant);
        
        // Find the nearest root that lies in the acceptable range.
        auto root = (-half_b - sqrtd) / a;
        if (root < t_min || t_max < root) {
            root = (-half_b + sqrtd) / a;
            if (root < t_min || t_max < root)
                return false;
        }
        
        rec.t = root;
        rec.p = r.at(rec.t);
        Vec3 outward_normal = (rec.p - center) / radius;
        rec.set_face_normal(r, outward_normal);
        rec.mat_ptr = mat_ptr;
        
        return true;
    }
};

// Hittable list class
class HittableList : public Hittable {
public:
    std::vector<std::shared_ptr<Hittable>> objects;
    
    HittableList() {}
    HittableList(std::shared_ptr<Hittable> object) { add(object); }
    
    void clear() { objects.clear(); }
    void add(std::shared_ptr<Hittable> object) { objects.push_back(object); }
    
    virtual bool hit(const Ray& r, double t_min, double t_max, HitRecord& rec) const override {
        HitRecord temp_rec;
        bool hit_anything = false;
        auto closest_so_far = t_max;
        
        for (const auto& object : objects) {
            if (object->hit(r, t_min, closest_so_far, temp_rec)) {
                hit_anything = true;
                closest_so_far = temp_rec.t;
                rec = temp_rec;
            }
        }
        
        return hit_anything;
    }
};

// Camera class
class Camera {
public:
    Point3 origin;
    Point3 lower_left_corner;
    Vec3 horizontal;
    Vec3 vertical;
    
    Camera(Point3 lookfrom, Point3 lookat, Vec3 vup, double vfov, double aspect_ratio) {
        auto theta = vfov * M_PI / 180.0;
        auto h = std::tan(theta/2);
        auto viewport_height = 2.0 * h;
        auto viewport_width = aspect_ratio * viewport_height;
        
        auto w = (lookfrom - lookat).unit_vector();
        auto u = cross(vup, w).unit_vector();
        auto v = cross(w, u);
        
        origin = lookfrom;
        horizontal = viewport_width * u;
        vertical = viewport_height * v;
        lower_left_corner = origin - horizontal/2 - vertical/2 - w;
    }
    
    Ray get_ray(double s, double t) const {
        return Ray(origin, lower_left_corner + s*horizontal + t*vertical - origin);
    }
};

// Utility function for random double
double random_double() {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<> dis(0.0, 1.0);
    return dis(gen);
}

// Ray color function
Color ray_color(const Ray& r, const Hittable& world, int depth) {
    HitRecord rec;
    
    // If we've exceeded the ray bounce limit, no more light is gathered.
    if (depth <= 0)
        return Color(0,0,0);
    
    if (world.hit(r, 0.001, std::numeric_limits<double>::infinity(), rec)) {
        Ray scattered;
        Color attenuation;
        if (rec.mat_ptr->scatter(r, rec, attenuation, scattered))
            return attenuation * ray_color(scattered, world, depth-1);
        return Color(0,0,0);
    }
    
    Vec3 unit_direction = r.direction.unit_vector();
    auto t = 0.5*(unit_direction.y + 1.0);
    return (1.0-t)*Color(1.0, 1.0, 1.0) + t*Color(0.5, 0.7, 1.0);
}

// Write color to output
void write_color(std::ostream &out, Color pixel_color, int samples_per_pixel) {
    auto r = pixel_color.x;
    auto g = pixel_color.y;
    auto b = pixel_color.z;
    
    // Divide the color by the number of samples and gamma-correct for gamma=2.0.
    auto scale = 1.0 / samples_per_pixel;
    r = std::sqrt(scale * r);
    g = std::sqrt(scale * g);
    b = std::sqrt(scale * b);
    
    // Write the translated [0,255] value of each color component.
    auto clamp = [](double x, double min, double max) { return x < min ? min : (x > max ? max : x); };
    out << static_cast<int>(256 * clamp(r, 0.0, 0.999)) << ' '
        << static_cast<int>(256 * clamp(g, 0.0, 0.999)) << ' '
        << static_cast<int>(256 * clamp(b, 0.0, 0.999)) << '\n';
}

// Create custom scene
HittableList create_scene() {
    HittableList world;
    
    // Ground
    auto ground_material = std::make_shared<Lambertian>(Color(0.5, 0.5, 0.5));
    world.add(std::make_shared<Sphere>(Point3(0,-1000,0), 1000, ground_material));
    
    // Large central glass sphere
    auto glass_material = std::make_shared<Dielectric>(1.5);
    world.add(std::make_shared<Sphere>(Point3(0, 1, 0), 1.0, glass_material));
    
    // Large left diffuse sphere
    auto diffuse_material = std::make_shared<Lambertian>(Color(0.4, 0.2, 0.1));
    world.add(std::make_shared<Sphere>(Point3(-4, 1, 0), 1.0, diffuse_material));
    
    // Large right metal sphere
    auto metal_material = std::make_shared<Metal>(Color(0.7, 0.6, 0.5), 0.0);
    world.add(std::make_shared<Sphere>(Point3(4, 1, 0), 1.0, metal_material));
    
    // Additional smaller spheres for visual interest
    world.add(std::make_shared<Sphere>(Point3(2, 0.5, -1), 0.5, 
        std::make_shared<Lambertian>(Color(0.7, 0.3, 0.3))));
    
    world.add(std::make_shared<Sphere>(Point3(-2, 0.5, -1), 0.5, 
        std::make_shared<Metal>(Color(0.8, 0.8, 0.9), 0.3)));
    
    world.add(std::make_shared<Sphere>(Point3(0, 0.5, -2), 0.5, 
        std::make_shared<Dielectric>(1.5)));
    
    world.add(std::make_shared<Sphere>(Point3(1, 0.3, 1), 0.3, 
        std::make_shared<Lambertian>(Color(0.2, 0.7, 0.2))));
    
    world.add(std::make_shared<Sphere>(Point3(-1, 0.3, 1), 0.3, 
        std::make_shared<Metal>(Color(0.9, 0.7, 0.2), 0.1)));
    
    return world;
}

int main() {
    // Image
    const auto aspect_ratio = 16.0 / 9.0;
    const int image_width = 800;
    const int image_height = static_cast<int>(image_width / aspect_ratio);
    const int samples_per_pixel = 100;
    const int max_depth = 50;
    
    // World
    HittableList world = create_scene();
    
    // Camera - positioned for an interesting perspective
    Point3 lookfrom(13, 2, 3);
    Point3 lookat(0, 0, 0);
    Vec3 vup(0, 1, 0);
    
    Camera cam(lookfrom, lookat, vup, 20, aspect_ratio);
    
    // Render
    std::ofstream file("image.ppm");
    file << "P3\n" << image_width << " " << image_height << "\n255\n";
    
    std::cout << "Rendering " << image_width << "x" << image_height << " image with " 
              << samples_per_pixel << " samples per pixel..." << std::endl;
    
    for (int j = image_height-1; j >= 0; --j) {
        std::cout << "\rScanlines remaining: " << j << ' ' << std::flush;
        for (int i = 0; i < image_width; ++i) {
            Color pixel_color(0, 0, 0);
            for (int s = 0; s < samples_per_pixel; ++s) {
                auto u = (i + random_double()) / (image_width-1);
                auto v = (j + random_double()) / (image_height-1);
                Ray r = cam.get_ray(u, v);
                pixel_color += ray_color(r, world, max_depth);
            }
            write_color(file, pixel_color, samples_per_pixel);
        }
    }
    
    std::cout << "\nDone! Image saved as image.ppm" << std::endl;
    return 0;
}