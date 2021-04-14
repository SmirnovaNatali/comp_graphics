#include "Material.h"

SurfHit& SurfHit::operator= (const SurfHit &rhs)
{
  if (this == &rhs)
    return (*this);

  hit = rhs.hit;
  hitPoint = rhs.hitPoint;
  normal = rhs.normal;
  m_ptr = rhs.m_ptr;
  t = rhs.t;

  return (*this);
}

//////////////////

float3 reflect(const float3& v, const float3& normal)
{
    //если поверхность обладает отражающими свойствами, то строится  вторичный луч отражения.
    //направление луча определяется по закону отражения 
    //ray=i-2*n(n-i)
    return v - 2 * dot(v, normal) * normal;
}


//идеальное зеркало
bool IdealMirror::Scatter(const Ray &ray_in, const SurfHit &surf, float3 &attenuation, Ray &ray_out)
{
    //генерирует луч
    float3 reflection_dir = reflect(normalize(ray_in.d), surf.normal); 
    //смещаем исхождную точку для избежания перекрытия объекта
    ray_out = Ray(surf.hitPoint + surf.normal * 1e-4, reflection_dir); 
    attenuation = color;

    return (dot(ray_out.d, surf.normal) > 0);
}
//диффузное отражение

bool  Diffuse::Scatter(const Ray &ray_in, const SurfHit &surf, float3 &attenuation, Ray &ray_out)
{
    //косинус высчитывается через перемножение единичных векторов
    //смотрим на угол между нормалью точки перечесения и вектором, описывающим направление света.
    float kd = dot(normalize(surf.normal), ray_in.d);
    attenuation = color * max(kd, 0.0f);
    return false;
}

bool Light::Scatter(const Ray &ray_in, const SurfHit &surf, float3 &attenuation, Ray &ray_out)
{
    attenuation = attenuation * color;
    return true;
}
