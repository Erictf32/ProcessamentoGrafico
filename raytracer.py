import math
import random
from typing import List, Optional


# Utility functions and classes
class Vec3:
    def __init__(self, x: float, y: float, z: float):
        self.x = x
        self.y = y
        self.z = z

    def __add__(self, other):
        return Vec3(self.x + other.x, self.y + other.y, self.z + other.z)

    def __sub__(self, other):
        return Vec3(self.x - other.x, self.y - other.y, self.z - other.z)

    def __mul__(self, other):
        if isinstance(other, Vec3):
            return Vec3(self.x * other.x, self.y * other.y, self.z * other.z)
        else:
            return Vec3(self.x * other, self.y * other, self.z * other)

    __rmul__ = __mul__

    def __truediv__(self, t: float):
        return self * (1.0 / t)

    def length(self):
        return math.sqrt(self.length_squared())

    def length_squared(self):
        return self.x * self.x + self.y * self.y + self.z * self.z

    def unit(self):
        return self / self.length()

    def dot(self, other):
        return self.x * other.x + self.y * other.y + self.z * other.z

    def cross(self, other):
        return Vec3(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x,
        )

    def near_zero(self):
        s = 1e-8
        return abs(self.x) < s and abs(self.y) < s and abs(self.z) < s

    def __repr__(self):
        return f"Vec3({self.x}, {self.y}, {self.z})"

    def __neg__(self):
        return Vec3(-self.x, -self.y, -self.z)


def clamp(x: float, min_val: float, max_val: float):
    if x < min_val:
        return min_val
    if x > max_val:
        return max_val
    return x


def random_in_unit_sphere():
    while True:
        p = Vec3(random.uniform(-1, 1), random.uniform(-1, 1), random.uniform(-1, 1))
        if p.length_squared() >= 1:
            continue
        return p


def random_unit_vector():
    return random_in_unit_sphere().unit()


def reflect(v: Vec3, n: Vec3):
    return v - n * (2 * v.dot(n))


def refract(uv: Vec3, n: Vec3, etai_over_etat: float):
    cos_theta = min((-uv).dot(n), 1.0)
    r_out_perp = (uv + n * cos_theta) * etai_over_etat
    r_out_parallel = n * (-math.sqrt(abs(1.0 - r_out_perp.length_squared())))
    return r_out_perp + r_out_parallel


def schlick(cosine: float, ref_idx: float):
    r0 = (1 - ref_idx) / (1 + ref_idx)
    r0 = r0 * r0
    return r0 + (1 - r0) * pow((1 - cosine), 5)


class Ray:
    def __init__(self, origin: Vec3, direction: Vec3):
        self.origin = origin
        self.direction = direction

    def at(self, t: float):
        return self.origin + self.direction * t


# Hittable and sphere
class HitRecord:
    def __init__(self):
        self.p: Optional[Vec3] = None
        self.normal: Optional[Vec3] = None
        self.t: float = 0.0
        self.front_face: bool = True
        self.material = None

    def set_face_normal(self, r: Ray, outward_normal: Vec3):
        self.front_face = r.direction.dot(outward_normal) < 0
        self.normal = outward_normal if self.front_face else outward_normal * -1


class Hittable:
    def hit(self, r: Ray, t_min: float, t_max: float, rec: HitRecord) -> bool:
        raise NotImplementedError()


class Sphere(Hittable):
    def __init__(self, center: Vec3, radius: float, material):
        self.center = center
        self.radius = radius
        self.material = material

    def hit(self, r: Ray, t_min: float, t_max: float, rec: HitRecord) -> bool:
        oc = r.origin - self.center
        a = r.direction.length_squared()
        half_b = oc.dot(r.direction)
        c = oc.length_squared() - self.radius * self.radius
        discriminant = half_b * half_b - a * c
        if discriminant < 0:
            return False
        sqrtd = math.sqrt(discriminant)

        root = (-half_b - sqrtd) / a
        if root < t_min or root > t_max:
            root = (-half_b + sqrtd) / a
            if root < t_min or root > t_max:
                return False

        rec.t = root
        rec.p = r.at(rec.t)
        outward_normal = (rec.p - self.center) / self.radius
        rec.set_face_normal(r, outward_normal)
        rec.material = self.material
        return True


class HittableList(Hittable):
    def __init__(self, objects: Optional[List[Hittable]] = None):
        self.objects: List[Hittable] = objects or []

    def add(self, obj: Hittable):
        self.objects.append(obj)

    def hit(self, r: Ray, t_min: float, t_max: float, rec: HitRecord) -> bool:
        temp_rec = HitRecord()
        hit_anything = False
        closest_so_far = t_max

        for obj in self.objects:
            if obj.hit(r, t_min, closest_so_far, temp_rec):
                hit_anything = True
                closest_so_far = temp_rec.t
                rec.t = temp_rec.t
                rec.p = temp_rec.p
                rec.normal = temp_rec.normal
                rec.front_face = temp_rec.front_face
                rec.material = temp_rec.material

        return hit_anything


# Materials
class Material:
    def scatter(self, r_in: Ray, rec: HitRecord):
        raise NotImplementedError()


class Lambertian(Material):
    def __init__(self, albedo: Vec3):
        self.albedo = albedo

    def scatter(self, r_in: Ray, rec: HitRecord):
        scatter_direction = rec.normal + random_unit_vector()
        if scatter_direction.near_zero():
            scatter_direction = rec.normal
        scattered = Ray(rec.p, scatter_direction)
        attenuation = self.albedo
        return True, attenuation, scattered


class Metal(Material):
    def __init__(self, albedo: Vec3, fuzz: float):
        self.albedo = albedo
        self.fuzz = min(fuzz, 1)

    def scatter(self, r_in: Ray, rec: HitRecord):
        reflected = reflect(r_in.direction.unit(), rec.normal)
        scattered = Ray(rec.p, reflected + random_in_unit_sphere() * self.fuzz)
        attenuation = self.albedo
        return (scattered.direction.dot(rec.normal) > 0), attenuation, scattered


class Dielectric(Material):
    def __init__(self, ref_idx: float):
        self.ref_idx = ref_idx

    def scatter(self, r_in: Ray, rec: HitRecord):
        attenuation = Vec3(1.0, 1.0, 1.0)
        etai_over_etat = 1 / self.ref_idx if rec.front_face else self.ref_idx

        unit_direction = r_in.direction.unit()
        cos_theta = min((-unit_direction).dot(rec.normal), 1.0)
        sin_theta = math.sqrt(1.0 - cos_theta * cos_theta)

        cannot_refract = etai_over_etat * sin_theta > 1.0
        if cannot_refract or schlick(cos_theta, etai_over_etat) > random.random():
            direction = reflect(unit_direction, rec.normal)
        else:
            direction = refract(unit_direction, rec.normal, etai_over_etat)

        scattered = Ray(rec.p, direction)
        return True, attenuation, scattered


# Camera
class Camera:
    def __init__(
        self,
        lookfrom: Vec3,
        lookat: Vec3,
        vup: Vec3,
        vfov: float,
        aspect_ratio: float,
        aperture: float,
        focus_dist: float,
    ):
        self.origin = lookfrom
        theta = math.radians(vfov)
        h = math.tan(theta / 2)
        viewport_height = 2.0 * h
        viewport_width = aspect_ratio * viewport_height

        w = (lookfrom - lookat).unit()
        u = vup.cross(w).unit()
        v = w.cross(u)

        self.horizontal = u * viewport_width * focus_dist
        self.vertical = v * viewport_height * focus_dist
        self.lower_left_corner = (
            self.origin - self.horizontal / 2 - self.vertical / 2 - w * focus_dist
        )

        self.lens_radius = aperture / 2
        self.u = u
        self.v = v

    def get_ray(self, s: float, t: float):
        rd = random_in_unit_sphere() * self.lens_radius
        offset = self.u * rd.x + self.v * rd.y
        return Ray(
            self.origin + offset,
            (self.lower_left_corner + self.horizontal * s + self.vertical * t)
            - self.origin
            - offset,
        )


def ray_color(r: Ray, world: Hittable, depth: int):
    if depth <= 0:
        return Vec3(0, 0, 0)

    rec = HitRecord()
    if world.hit(r, 0.001, math.inf, rec):
        scatter_result = rec.material.scatter(r, rec)
        if scatter_result[0]:
            attenuation, scattered = scatter_result[1], scatter_result[2]
            return attenuation * ray_color(scattered, world, depth - 1)
        return Vec3(0, 0, 0)

    unit_direction = r.direction.unit()
    t = 0.5 * (unit_direction.y + 1.0)
    return Vec3(1.0, 1.0, 1.0) * (1.0 - t) + Vec3(0.5, 0.7, 1.0) * t


# Render
if __name__ == "__main__":
    # Image
    aspect_ratio = 16.0 / 9.0
    image_width = 400
    image_height = int(image_width / aspect_ratio)
    samples_per_pixel = 100
    max_depth = 50

    # World
    world = HittableList()

    ground_material = Lambertian(Vec3(0.8, 0.8, 0.0))
    world.add(Sphere(Vec3(0, -1000.5, -1), 1000, ground_material))

    glass = Dielectric(1.5)
    center_glass = Sphere(Vec3(-0.7, 0.2, -1), 0.3, glass)
    world.add(center_glass)

    diffuse = Lambertian(Vec3(0.1, 0.2, 0.9))
    world.add(Sphere(Vec3(0.0, 0.2, -1.0), 0.2, diffuse))

    metal = Metal(Vec3(0.8, 0.6, 0.2), 0.0)
    world.add(Sphere(Vec3(0.7, 0.3, -1.5), 0.3, metal))

    # Camera
    lookfrom = Vec3(1.5, 1.0, 2.0)
    lookat = Vec3(0, 0.2, -1)
    vup = Vec3(0, 1, 0)
    dist_to_focus = (lookfrom - lookat).length()
    aperture = 0.1

    cam = Camera(lookfrom, lookat, vup, 40, aspect_ratio, aperture, dist_to_focus)

    # Render
    with open("final_scene.ppm", "w") as f:
        f.write(f"P3\n{image_width} {image_height}\n255\n")
        for j in range(image_height - 1, -1, -1):
            print(f"Scanlines remaining: {j}", end="\r")
            for i in range(image_width):
                color = Vec3(0, 0, 0)
                for s in range(samples_per_pixel):
                    u = (i + random.random()) / (image_width - 1)
                    v = (j + random.random()) / (image_height - 1)
                    r = cam.get_ray(u, v)
                    color += ray_color(r, world, max_depth)
                scale = 1.0 / samples_per_pixel
                r = math.sqrt(scale * color.x)
                g = math.sqrt(scale * color.y)
                b = math.sqrt(scale * color.z)
                ir = int(256 * clamp(r, 0.0, 0.999))
                ig = int(256 * clamp(g, 0.0, 0.999))
                ib = int(256 * clamp(b, 0.0, 0.999))
                f.write(f"{ir} {ig} {ib} ")
            f.write("\n")
    print("\nDone.")