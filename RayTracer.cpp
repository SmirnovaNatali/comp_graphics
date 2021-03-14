#include <limits>
#include "RayTracer.h"
#include "Material.h"
#include "Geometry.h"


//Базовый алгоритм трассировки луча
float3 SimpleRT::TraceRay(const Ray & ray, const std::vector<std::shared_ptr<GeoObject>> &geo, const int &depth)
{
  float tnear = std::numeric_limits<float>::max();
  int   geoIndex = -1;

  SurfHit surf;

  for (int i = 0; i < geo.size(); ++i)
  {
    SurfHit temp;
    //
    if (geo.at(i)->Intersect(ray, 0.001, tnear, temp))
    {
      if (temp.t < tnear)
      {
        tnear = temp.t;
        geoIndex = i;//перебирает вектор
       // если удовлетворяет , запоминаем его
        surf = temp;
      }
    }
  }


  if (geoIndex == -1)
  {
    float3 unit_direction = normalize(ray.d);
    float t = 0.5f * (unit_direction.y + 1.0f);

    return (1.0f - t) * float3(1.0f, 1.0f, 1.0f) + t * bg_color;
  }

  float3 surfColor(0.0f, 0.0f, 0.0f);
  if (dot(ray.d, surf.normal) > 0)
  {
    surf.normal = -surf.normal;
  }

  Ray scattered;// новый луч
  if(depth < max_ray_depth && surf.m_ptr->Scatter(ray, surf, surfColor, scattered))
  {
    return surfColor * TraceRay(scattered, geo, depth + 1);//если глубина не достигнута предела и есть пересечения, пускаем луч
  }
  else
  {
    return float3(0.0f, 0.0f, 0.0f);
  }
}

float3 WhittRT::TraceRay(const Ray& ray, const std::vector<std::shared_ptr<GeoObject>>&geo, const std::vector<std::shared_ptr<Lights>>&lig, int depth)
{
  float3 color = float3(1.0f, 1.0f, 1.0f);// белый цвет
  float3 timeсolor = float3(1.0f, 1.0f, 1.0f);//белый цвет
  //первичный луч 
  Ray timeofray = ray;
  while (depth < max_depth) {
	color *= timeсolor;
	float tnear = std::numeric_limits<float>::max();
    int   geoIndex = -1;
	SurfHit surf;
	//найти ближайшее пересечение между лучом и объектом
	for (int i = 0; i < geo.size(); ++i)
	{
		SurfHit temp;
		//пересечение найдено
		if (geo.at(i)->Intersect(timeofray, 0.001, tnear, temp))
		{
			if (temp.t < tnear)
			{
				tnear = temp.t;
				geoIndex = i;
				//перебирает вектор
	            // если удовлетворяет , запоминаем его
				surf = temp;
			}
		}
	}
		// Если пересечение не найдено
	if (geoIndex == -1)
	{
	  float3 unit_direction = normalize(timeofray.d);
	  float t = 0.5f * (unit_direction.y + 1.0f);
      timeсolor = (1.0f - t) * float3(1.0f, 1.0f, 1.0f) + t * bg_color;
	  depth++;//идем по глубине
	  break;
	}

	if (dot(timeofray.d, surf.normal) > 0)
	{
	  surf.normal = -surf.normal;
	}
	Ray scattered;//новый луч

	// Если не источник света
	if (typeid(*surf.m_ptr) != typeid(Light))
	{
		//диффузный материал
		if (typeid(*surf.m_ptr).name() == typeid(Defuse).name()) 
		{
			timeсolor = float3(0.0f, 0.0f, 0.0f);
		    float3 time;
			int countoflig = 0;
			for (int i = 0; i < lig.size(); i++) {
				Ray rayIn;
				rayIn.o = lig.at(i)->position;
				rayIn.d = normalize(rayIn.o - surf.hitPoint);
				Ray shadow(surf.hitPoint + normalize(surf.normal) * 1e-4, rayIn.d);
				// на пути луча нет проблем
				if (!ShadowRay(shadow, geo))
				{
					surf.m_ptr->Scatter(rayIn, surf, time, scattered);
					timeсolor += time * lig.at(i)->color;
					countoflig++;
				}
			}
			break;
		}
	// Иначе находим пересечение с зеркалом
	else if (surf.m_ptr->Scatter(timeofray, surf, timeсolor, scattered))
	{
	   timeofray = scattered;
	   depth++;
	}
	else
	{
	   depth++;
	   timeсolor = float3(0.0f, 0.0f, 0.0f);
	}
	}
    // Если источник света
	else 
    {
		surf.m_ptr->Scatter(timeofray, surf, timeсolor, scattered);
		break;
	}
  }
  color *= timeсolor;
  return color;

}

// теневой луч показывает,есть ли объекты на пути луча от точки пересечения, которую нашел
bool WhittRT::ShadowRay(const Ray& ray, const std::vector<std::shared_ptr<GeoObject>>& geo) {
  Ray timeRay = ray;
  float tnear = std::numeric_limits<float>::max();
  int   geoIndex = -1;
  SurfHit surf;
  for (int i = 0; i < geo.size(); ++i)
  {
	SurfHit temp;
	if (geo.at(i)->Intersect(timeRay, 0.001, tnear, temp))
	{
		if (temp.t < tnear)
		{
		   tnear = temp.t;
			geoIndex = i;
			//перебирает вектор
			 // если удовлетворяет , запоминаем его
			surf = temp;
		}
		return true;
	}

  }
  return false;
}
