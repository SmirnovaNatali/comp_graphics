#include "Geometry.h"


bool Plane::Intersect(const Ray &ray, float t_min, float t_max, SurfHit &surf) const
{
    //находим время
  surf.t = dot((point - ray.o), normal) / dot(ray.d, normal);

  if(surf.t > t_min && surf.t < t_max)
  {
    surf.hit      = true;
    surf.hitPoint = ray.o + surf.t * ray.d;//формула луча
    surf.normal   = normal;
    surf.m_ptr    = m_ptr;// материал объекта пересечения
    return true;
  }

  return false;
}

/////////////////////////////////////////
bool Sphere::Intersect(const Ray& ray, float t_min, float t_max, SurfHit& surf) const
{
    //считаем по формуле
    float3 k = ray.o - center;
    //пересчение луча с базовой сферой
    //квадрат вектора скорости
    float a = dot(ray.d, ray.d);
    //точка испускания луча на вектор скорости на 2
    float b = dot(2 * k, ray.d);
    //квалрат (исходная точка-центр сферы)- квадрат радиуса
    float c = dot(k, k) - r_sq;
    //дискриминант
    float diskr = b * b - 4 * a * c;
    //смотрим на корни
    if (diskr < 0) 
        return false;
    //D < 0 – пересечения нет
    // D = 0 – касание в одной точке
    //D > 0 – пересечение в двух точках
    //1 корень
    surf.t = (-b - sqrt(diskr)) / 2 * a;
    if (surf.t > t_min && surf.t < t_max)
    {
        surf.hit = true;
        surf.hitPoint = ray.o + surf.t * ray.d;
        surf.normal = normalize(surf.hitPoint - center);
        surf.m_ptr = m_ptr;
        return true;
    }

    //2 корень
    surf.t = (-b + sqrt(diskr)) / 2 * a;
    if (surf.t > t_min && surf.t < t_max)
    {
        surf.hit = true;
        surf.hitPoint = ray.o + surf.t * ray.d;
        surf.normal = normalize(surf.hitPoint - center);
        surf.m_ptr = m_ptr;
        return true;
    }

    return false;
}


bool Triangle::Intersect(const Ray& ray, float tmin, float tmax, SurfHit& surf) const
{
    //из Алгоритма Моллера-Трумбора
    float3 T = ray.o - a;//пересечение луча
    float3 E1 = b - a;//вектор 1 стороны
    float3 E2 = c - a;//2 стороны
    float3 D = ray.d;
    float3 P = cross(D, E2);
    float3 Q = cross(T, E1);
    float det = dot(E1, P);//время



    if (det < tmin && det > tmax) {
        return false;
    }
    //u и v должны удовлетворять условиям, меньше 1, больше 0
    float invDet = 1 / det;
    float u = dot(T, P) * invDet;

    if (u < 0 || u > 1) {
        return false;
    }

    float v = dot(ray.d, Q) * invDet;

    if (v < 0 || u + v > 1) {
        return false;
    }


    surf.t = dot(E2, Q) * invDet;
    if (surf.t > tmin && surf.t < tmax) {
        surf.hit = true;
        surf.hitPoint = float3(surf.t, u, v);
        surf.normal = cross(b - a, c - a);//нормаль через умнжение
        surf.m_ptr = m_ptr;
        return true;
    }


    return false;
}

//выровняем по осям координат,строим по 2 точкам
bool Parallel::Intersect(const Ray& ray, float tmin, float tmax, SurfHit& surf) const
{
    //противоположные грани прям.парал. лежат в плоскостях, которые парал координатным плоскостям
   //находим для координаты t
    //координаты вершин параллепипида- это коорд.источника луча,делен на координаты вектора времени

    float x1 = (t_min.x - ray.o.x) / ray.d.x;
    float x2 =(t_max.x - ray.o.x) / ray.d.x;
    float y1 = (t_min.y - ray.o.y) / ray.d.y;
    float y2 = (t_max.y - ray.o.y) / ray.d.y;
    float z1 = (t_min.z - ray.o.z) / ray.d.z;
    float z2 = (t_max.z - ray.o.z) /ray.d.z;

    //оба парметра отрицатлеьны, тогда луч не пересекаетэту пару плосксотей, не пересекает параллеп
    //необходимо взять максимум по осям для нашего ближайшего расстояние tmin
    float tMin = max(max(min(x1, x2), min(y1, y2)), min(z1, z2));
    float tMax = min(min(max(x1, x2), max(y1, y2)), max(z1, z2));


    surf.t = tMin;
    //если tmin<=tmax и tmax>0, пересечение с кубоидом
    if (tMin < tMax && tMax > 0 && surf.t > tmin && surf.t < tmax) {
        surf.hit = true;
        surf.hitPoint = ray.o + surf.t * ray.d;
        surf.normal = normalize(surf.hitPoint);
        surf.m_ptr = m_ptr;
        return true;
    }
    return false;
}


bool Square::Intersect(const Ray& ray, float tmin, float tmax, SurfHit& surf) const
{
   
    float3 d = float3((c.x+a.x-b.x), (c.y + a.y - b.y), (c.z + a.z - b.z)); // находим d по коорлдинатам, задаем по уравнению плоскостей
    
    if (Triangle(a, b, c, new IdealMirror(float3(1.00f, 1.0f, 0.00f))).Intersect(ray, tmin, tmax, surf))
        return true;

    if (Triangle(a, c, d, new IdealMirror(float3(1.00f, 1.0f, 0.00f))).Intersect(ray, tmin, tmax, surf))
        return true;


    return false;
}

